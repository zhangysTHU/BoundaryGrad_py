from __future__ import annotations

# 07：空间 spot 反卷积，估计每个 spot 中不同参考细胞类型的比例。
# 原理分两层：
# 1) 用 marker enrichment 判断每个 Location_cluster topic 可能有哪些细胞类型；
# 2) 对候选细胞类型做非负且和为 1 的最小二乘拟合，得到比例矩阵 DeconData。
# 输入：
# - intermediate/05_TumorST_boundary_defined.h5ad：需要 obs["Location"] 和 seurat_clusters。
# - intermediate/06_sig_exp.csv：gene x celltype signature。
# - intermediate/06_clustermarkers_list.json：celltype 到 marker genes 的映射。
# 输出：
# - intermediate/07_DeconData.tsv：cell_ID + 各 celltype 比例；08 重构和 10 pie/bar 图直接使用。
# - intermediate/07_TumorST_for_decon.h5ad：保存归一化后、带 Decon_topics 的对象，供复核。
# - output/07_spatial_deconvolution/DeconData.xlsx：同一反卷积结果的 Excel 版本，方便手动查看。
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.config_loader import load_config
from lib.decon_helpers import deconvolve_topic, enrich_analysis, get_enrich_matrix

cfg = load_config()
paths = cfg.PATHS
out_dir = paths["output"] / "07_spatial_deconvolution"
out_dir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(paths["intermediate"] / "05_TumorST_boundary_defined.h5ad")
sig_exp = pd.read_csv(paths["intermediate"] / "06_sig_exp.csv", index_col=0)
markers = json.loads((paths["intermediate"] / "06_clustermarkers_list.json").read_text(encoding="utf-8"))

# 反卷积使用 log-normalized 表达；nCount_Spatial 仍来自原始 QC，用于辅助判断低测序深度 spot。
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Decon_topics 把空间区域 Location 和形态/表达 cluster 合并，减少逐 spot 独立判断的噪声。
adata.obs["Decon_topics"] = adata.obs["Location"].astype(str) + "_" + adata.obs["seurat_clusters"].astype(str)
expr = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
log_expr = pd.DataFrame(expr.T, index=adata.var_names, columns=adata.obs_names)
nolog_expr = np.power(2, log_expr) - 1

# meta 中每个参考细胞类型的列，是该细胞类型 top marker 在 spot 中的平均表达，用于 topic 级筛选。
meta = adata.obs[["nCount_Spatial", "Decon_topics", "Location"]].copy()
for cluster, genes in markers.items():
    present = [g for g in genes[:25] if g in log_expr.index]
    meta[cluster] = log_expr.loc[present].mean(axis=0) if present else 0.0

# 只保留空间表达和单细胞 signature 都存在的基因。
genes = sig_exp.index.intersection(nolog_expr.index)
filter_sig = sig_exp.loc[genes]
filter_expr = nolog_expr.loc[genes]
filter_log_expr = log_expr.loc[genes]
enrich_matrix = get_enrich_matrix(filter_sig, markers)
enrich_result = enrich_analysis(filter_log_expr, enrich_matrix)

# decon 行为 spot，列为参考细胞类型；每一行最终近似和为 1。
decon = pd.DataFrame(0.0, index=adata.obs_names, columns=filter_sig.columns)
for topic, topic_obs in meta.groupby("Decon_topics"):
    topic_cells = topic_obs.index.tolist()
    # topic 内先按 marker enrichment 选候选细胞类型，避免所有类型一起拟合导致不稳定。
    avg_enrich = enrich_result[topic_cells].max(axis=1).sort_values(ascending=False)
    selected = list(avg_enrich.head(3).index)
    if topic.startswith("Mal_") and cfg.PARAMS["decon_malignant_cluster"] in filter_sig.columns:
        selected.append(cfg.PARAMS["decon_malignant_cluster"])
    if topic_obs["nCount_Spatial"].median() < 5000 and cfg.PARAMS["decon_stromal_cluster"] in filter_sig.columns:
        selected.append(cfg.PARAMS["decon_stromal_cluster"])
    selected = [x for x in dict.fromkeys(selected) if x in filter_sig.columns]
    if len(selected) == 1:
        decon.loc[topic_cells, selected[0]] = 1
        continue
    # 用候选类型 marker 基因做 constrained least squares，估计每个 spot 的组成。
    topic_markers = sorted(set(g for ct in selected for g in markers.get(ct, []) if g in filter_expr.index))
    if not topic_markers:
        topic_markers = filter_expr.index.tolist()
    res = deconvolve_topic(filter_expr.loc[topic_markers, topic_cells], filter_sig.loc[topic_markers, selected])
    decon.loc[topic_cells, selected] = res.T.loc[topic_cells, selected]

# 输出宽表：cell_ID + 各参考细胞类型比例。
decon.insert(0, "cell_ID", decon.index)
decon.to_csv(paths["intermediate"] / "07_DeconData.tsv", sep="\t", index=False)
decon.to_excel(out_dir / "DeconData.xlsx", index=False)
adata.write_h5ad(paths["intermediate"] / "07_TumorST_for_decon.h5ad")
