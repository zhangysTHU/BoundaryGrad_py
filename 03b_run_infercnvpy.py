from __future__ import annotations

# 03b：纯 Python 的 infercnvpy 路线。
# 它按基因组位置对表达做滑窗平滑，估计每个 spot 的 CNV profile，
# 再用 KMeans 把 CNV profile 分成若干 CNVLabel，并输出统一格式的 03_cnv_calls.tsv。
# 输入：intermediate/02_TumorST_clustered.h5ad 和 resources/gencode_v38_gene_pos.txt。
# 输出：
# - intermediate/03_cnv_calls.tsv：cell_ID、CNVLabel、cnv_score，供 04 读取。
# - intermediate/03_infercnvpy_result.h5ad：包含 infercnvpy 结果 obsm["X_cnv"] 的 AnnData，供复核。
# - output/03_infercnv/infercnvpy/cnv_calls.tsv：03_cnv_calls.tsv 的输出目录副本。
# 备注：COTTRAZM_FAST_CNV=1 只生成结构正确的 smoke-test 文件，不代表真实 CNV。
import sys
import os
import threading
import time
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.cluster import KMeans

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.config_loader import load_config
from lib.io_utils import save_cnv_calls

cfg = load_config()
paths = cfg.PATHS
out_dir = paths["output"] / "03_infercnv" / "infercnvpy"
out_dir.mkdir(parents=True, exist_ok=True)


def log_step(message: str) -> None:
    # 带时间戳的日志，便于长时间运行时确认停在哪个大步骤。
    print(f"{datetime.now():%Y-%m-%d %H:%M:%S} | {message}", flush=True)


@contextmanager
def heartbeat(label: str, interval_seconds: int = 120):
    # infercnvpy.tl.infercnv 是一个独立耗时块，内部不持续吐进度；
    # heartbeat 至少能证明进程仍在该块中运行。
    stop = threading.Event()
    started = time.time()

    def beat() -> None:
        while not stop.wait(interval_seconds):
            elapsed = int(time.time() - started)
            log_step(f"{label} still running; elapsed {elapsed // 60} min {elapsed % 60} sec")

    thread = threading.Thread(target=beat, daemon=True)
    thread.start()
    try:
        yield
    finally:
        stop.set()
        thread.join(timeout=1)


def env_int(name: str, default: int) -> int:
    # 用环境变量控制并行数/进度间隔，便于在不同机器上调低负载。
    value = os.getenv(name)
    if value in (None, ""):
        return default
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"{name} must be an integer, got {value!r}") from exc
    if parsed < 1:
        raise ValueError(f"{name} must be >= 1, got {parsed}")
    return parsed


try:
    import infercnvpy as cnv
except ImportError as exc:
    raise ImportError("03b requires infercnvpy. Install it or run 03a_run_infercnv_external_r.py instead.") from exc


def main() -> None:
    log_step("loading clustered AnnData")
    adata = sc.read_h5ad(paths["intermediate"] / "02_TumorST_clustered.h5ad")
    log_step(f"loaded object: {adata.n_obs} spots x {adata.n_vars} genes")
    if os.getenv("COTTRAZM_FAST_CNV") == "1":
        # smoke-test 模式：跳过真实 CNV 推断，只生成结构正确的 CNV calls。
        # 适合测试 04-10 是否能跑通，但不能用于真实生物学结论。
        log_step("COTTRAZM_FAST_CNV=1 detected; writing smoke-test CNV calls")
        normal_cluster = adata.obs.groupby("seurat_clusters")["NormalScore"].mean().sort_values(ascending=False).index[0]
        scores = adata.obs["nCount_Spatial"].astype(float)
        scores = scores / scores.max()
        calls = pd.DataFrame({"cell_ID": adata.obs_names, "CNVLabel": adata.obs["seurat_clusters"].astype(str).values, "cnv_score": scores.values})
        calls.loc[calls["CNVLabel"] == str(normal_cluster), ["CNVLabel", "cnv_score"]] = ["Normal", 0.0]
        save_cnv_calls(paths["intermediate"] / "03_cnv_calls.tsv", calls)
        calls.to_csv(out_dir / "cnv_calls.tsv", sep="\t", index=False)
        log_step(f"wrote {paths['intermediate'] / '03_cnv_calls.tsv'}")
        raise SystemExit(0)

    # infercnvpy 需要每个基因的染色体坐标。这里和 gencode_v38_gene_pos.txt 取交集。
    log_step("loading gene position table")
    gene_pos = pd.read_csv(paths["resources"] / "gencode_v38_gene_pos.txt", sep="\t", header=None)
    gene_pos = gene_pos.iloc[:, :4]
    gene_pos.columns = ["gene", "chromosome", "start", "end"]
    gene_pos = gene_pos.drop_duplicates("gene").set_index("gene")
    common = adata.var_names.intersection(gene_pos.index)
    log_step(f"gene position overlap: {len(common)} genes")
    log_step("subsetting AnnData and joining gene coordinates")
    adata = adata[:, common].copy()
    adata.var = adata.var.join(gene_pos.loc[common])
    adata.var["start"] = pd.to_numeric(adata.var["start"], errors="coerce")
    adata.var["end"] = pd.to_numeric(adata.var["end"], errors="coerce")

    # 参考组选择沿用 Cottrazm 思路：NormalScore 最高的空间 cluster 作为 normal/reference。
    normal_cluster = adata.obs.groupby("seurat_clusters")["NormalScore"].mean().sort_values(ascending=False).index[0]
    adata.obs["infercnv_reference"] = np.where(adata.obs["seurat_clusters"].astype(str) == str(normal_cluster), "reference", "observation")
    log_step(f"reference cluster selected: {normal_cluster}")

    infercnv_n_jobs = env_int("COTTRAZM_INFERCNVPY_N_JOBS", 4)
    infercnv_chunksize = env_int("COTTRAZM_INFERCNVPY_CHUNKSIZE", 1000)
    heartbeat_seconds = env_int("COTTRAZM_PROGRESS_SECONDS", 120)
    log_step(
        "starting infercnvpy.tl.infercnv "
        f"(window_size=100, step=10, n_jobs={infercnv_n_jobs}, chunksize={infercnv_chunksize})"
    )
    with heartbeat("infercnvpy.tl.infercnv", heartbeat_seconds):
        # 核心 CNV 推断：按基因组顺序滑窗聚合表达，和 reference 比较得到 CNV-like signal。
        cnv.tl.infercnv(
            adata,
            reference_key="infercnv_reference",
            reference_cat=["reference"],
            window_size=100,
            n_jobs=infercnv_n_jobs,
            chunksize=infercnv_chunksize,
        )
    log_step("infercnvpy.tl.infercnv completed")

    # infercnvpy 把结果放在 obsm["X_cnv"]，每行对应 spot，每列对应 CNV window/feature。
    cnv_matrix = adata.obsm.get("X_cnv")
    if cnv_matrix is None:
        raise RuntimeError("infercnvpy finished but did not produce adata.obsm['X_cnv'].")
    if sp.issparse(cnv_matrix):
        cnv_matrix = cnv_matrix.toarray()
    else:
        cnv_matrix = np.asarray(cnv_matrix)
    if cnv_matrix.ndim != 2:
        raise RuntimeError(f"Expected a 2D CNV matrix, got shape {cnv_matrix.shape}.")
    log_step(f"CNV matrix ready: {cnv_matrix.shape[0]} spots x {cnv_matrix.shape[1]} windows/features")
    scores = np.abs(cnv_matrix).sum(axis=1)

    # 用 KMeans 把连续 CNV profile 离散成 cnv_k 个标签；reference spot 强制标为 Normal。
    log_step("running KMeans on CNV matrix")
    labels = KMeans(n_clusters=cfg.PARAMS["cnv_k"], random_state=0, n_init=10).fit_predict(cnv_matrix) + 1
    calls = pd.DataFrame({"cell_ID": adata.obs_names, "CNVLabel": labels.astype(str), "cnv_score": scores})
    calls.loc[adata.obs["infercnv_reference"].values == "reference", ["CNVLabel", "cnv_score"]] = ["Normal", 0.0]

    log_step("saving CNV calls and AnnData result")
    save_cnv_calls(paths["intermediate"] / "03_cnv_calls.tsv", calls)
    adata.write_h5ad(paths["intermediate"] / "03_infercnvpy_result.h5ad")
    calls.to_csv(out_dir / "cnv_calls.tsv", sep="\t", index=False)
    log_step("step 03b completed")


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
