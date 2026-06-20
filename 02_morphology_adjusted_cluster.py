from __future__ import annotations

# 02：用 stLearn 的 H&E 形态特征校正空间表达，再做 PCA/邻居图/UMAP/Leiden 聚类。
# 输入：intermediate/01_TumorST_preprocessed.h5ad，里面 X/counts 是 spot x gene 表达矩阵。
# 输出：
# - intermediate/02_TumorST_clustered.h5ad：加入 Morph 层、X_umap、seurat_clusters、NormalScore，供 03/04/05 使用。
# - intermediate/InferCNV/raw_counts_matrix.mtx：gene x cell 稀疏矩阵，供 03a 外部 R inferCNV 使用。
# - intermediate/InferCNV/genes.tsv / cells.tsv：矩阵行列名。
# - intermediate/InferCNV/CellAnnotation.txt：两列无表头，cell_ID 和 cluster，用于 inferCNV 分组。
# - intermediate/InferCNV/reference_cluster.txt：NormalScore 最高的 cluster，作为 inferCNV reference。
# - output/02_morphology_cluster/*：聚类图、NormalScore 图和 stLearn tile 图像。
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.config_loader import load_config
from lib.io_utils import ensure_dir, write_infercnv_inputs
from lib.plot_utils import save_embedding_scatter, save_spatial_scatter

cfg = load_config()
paths = cfg.PATHS
out_dir = ensure_dir(paths["output"] / "02_morphology_cluster")

adata = sc.read_h5ad(paths["intermediate"] / "01_TumorST_preprocessed.h5ad")

# stLearn 会切 tile、提取图像特征，并用 SME_normalize 融合空间邻域、形态和表达。
# 若 stLearn 或图像处理失败，则退回 raw counts，保证流水线仍可继续。
try:
    import stlearn as st
    from pathlib import Path as _Path

    st_data = st.Read10X(path=str(paths["spaceranger"]))
    st_data.var_names_make_unique()
    st_data.layers["raw_count"] = st_data.X
    tile_path = _Path(out_dir / f"{cfg.SAMPLE_NAME}_tile")
    tile_path.mkdir(parents=True, exist_ok=True)
    st.pp.tiling(st_data, tile_path, crop_size=40)
    st.pp.extract_feature(st_data)
    st.pp.normalize_total(st_data)
    st.pp.log1p(st_data)
    st.em.run_pca(st_data, n_comps=50, random_state=0)
    st.spatial.SME.SME_normalize(st_data, use_data="raw", weights="weights_matrix_gd_md")
    morph = sparse.csr_matrix(np.asarray(st_data.obsm["raw_SME_normalized"]))
except Exception as exc:
    print(f"stlearn morphology normalization failed, falling back to raw counts: {exc}")
    morph = adata.layers["counts"].copy() if "counts" in adata.layers else adata.X.copy()

# Morph 层保存形态校正后的矩阵；X 暂时替换为 Morph，用于聚类。
adata.layers["Morph"] = morph
adata.X = morph.copy()

# 标准单细胞/空间转录组降维聚类流程：
# normalize_total -> log1p -> 高变基因 -> scale -> PCA -> 邻居图 -> UMAP -> Leiden。
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_pcs=50)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=cfg.PARAMS["cluster_resolution"], key_added="seurat_clusters")

# Cottrazm 原流程用免疫/B 细胞等正常细胞 marker 给 cluster 打 NormalScore，
# 后面 inferCNV 会选 NormalScore 最高的 cluster 作为参考。
normal_features = ["PTPRC", "CD2", "CD3D", "CD3E", "CD3G", "CD5", "CD7", "CD79A", "MS4A1", "CD19"]
present = [g for g in normal_features if g in adata.var_names]
if present:
    vals = adata[:, present].X
    adata.obs["NormalScore"] = np.asarray(vals.mean(axis=1)).ravel()
else:
    adata.obs["NormalScore"] = 0.0

save_spatial_scatter(adata, "seurat_clusters", out_dir / f"{cfg.SAMPLE_NAME}_Spatial_SeuratCluster.pdf", cfg.CLUSTER_COLS)
save_embedding_scatter(adata, "X_umap", "seurat_clusters", out_dir / f"{cfg.SAMPLE_NAME}_UMAP_SeuratCluster.pdf", cfg.CLUSTER_COLS)

import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(7, 4))
adata.obs.boxplot(column="NormalScore", by="seurat_clusters", ax=ax)
ax.set_title("NormalScore")
fig.suptitle("")
fig.tight_layout()
fig.savefig(out_dir / f"{cfg.SAMPLE_NAME}_NormalScore.pdf")
plt.close(fig)

# 为 03a 外部 R inferCNV 路线写入 matrix/cells/genes/CellAnnotation。
infer_dir = ensure_dir(paths["intermediate"] / "InferCNV")
write_infercnv_inputs(adata, infer_dir, group_key="seurat_clusters")

# reference_cluster.txt 记录 inferCNV 参考 cluster。
normal_cluster = adata.obs.groupby("seurat_clusters")["NormalScore"].mean().sort_values(ascending=False).index[0]
(infer_dir / "reference_cluster.txt").write_text(str(normal_cluster), encoding="utf-8")
adata.write_h5ad(paths["intermediate"] / "02_TumorST_clustered.h5ad")
