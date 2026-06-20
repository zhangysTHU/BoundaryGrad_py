from __future__ import annotations

# 08：把 07 的 spot 级细胞比例和单细胞 signature 结合，重构指定区域的“伪单细胞/亚型表达矩阵”。
# 默认只重构 Bdy 边界区：每个 spot 会按反卷积比例拆成若干 Subtype_spot 列。
# 输入：
# - intermediate/05_TumorST_boundary_defined.h5ad：提供原始 counts、Location、spot 名。
# - intermediate/06_sig_exp.csv 和 06_clustermarkers_list.json：提供 signature 和重构基因集合。
# - intermediate/07_DeconData.tsv：提供每个 spot 的细胞类型比例。
# 输出：
# - intermediate/08_reconstructed_matrix.tsv：gene x reconstructed_cell，列名形如 celltype_spotbarcode。
# - intermediate/08_TumorST_reconstructed.h5ad：上述矩阵的 AnnData，obs 含 Subtypes、orig.ident、Location。
# 下游：当前 09/10 不直接依赖 08；08 主要用于边界区表达重构的单独分析。
import json
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.config_loader import load_config
from lib.recon_helpers import get_recon_matrix

cfg = load_config()
paths = cfg.PATHS
out_dir = paths["output"] / "08_spatial_reconstruction"
out_dir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(paths["intermediate"] / "05_TumorST_boundary_defined.h5ad")
sig_exp = pd.read_csv(paths["intermediate"] / "06_sig_exp.csv", index_col=0)
markers = json.loads((paths["intermediate"] / "06_clustermarkers_list.json").read_text(encoding="utf-8"))
decon = pd.read_csv(paths["intermediate"] / "07_DeconData.tsv", sep="\t").set_index("cell_ID")

# 使用原始 counts 做表达量分配，避免在 log 空间分解 UMI。
counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
counts = counts.toarray() if sparse.issparse(counts) else counts
counts = pd.DataFrame(counts.T, index=adata.var_names, columns=adata.obs_names)

# 核心重构：对每个目标 spot 和其中非零细胞类型，按 signature 权重分配该 spot 的 gene count。
mtx = get_recon_matrix(counts, adata.obs, sig_exp, markers, decon, cfg.PARAMS["recon_locations"])

# AnnData 行是重构出来的 Subtype_spot，列是基因。
recon = AnnData(X=mtx.T.values)
recon.obs_names = mtx.columns
recon.var_names = mtx.index
recon.obs["Subtypes"] = [x.split("_")[0].replace(".", "_") for x in recon.obs_names]
recon.obs["orig.ident"] = [x.split("_", 1)[1] if "_" in x else x for x in recon.obs_names]
recon.obs["Location"] = adata.obs.reindex(recon.obs["orig.ident"])["Location"].values
recon.write_h5ad(paths["intermediate"] / "08_TumorST_reconstructed.h5ad")
mtx.to_csv(paths["intermediate"] / "08_reconstructed_matrix.tsv", sep="\t")
