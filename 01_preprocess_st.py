from __future__ import annotations

# 01：读取 10x Visium / Space Ranger 输出，补齐空间坐标，计算基础 QC 指标。
# 输入：input/spaceranger_outs/，应包含 filtered_feature_bc_matrix(.h5 或目录) 和 spatial/ 图像坐标文件。
# 输出：
# - intermediate/01_TumorST_preprocessed.h5ad：AnnData，保留 layers["counts"]、obs QC 指标、obsm["spatial"]。
# - output/01_preprocess/QC/QCData.xlsx：每个 spot 的 nCount_Spatial、nFeature_Spatial、Mito.percent。
# - output/01_preprocess/QC/Vlnplot.pdf：三个 QC 指标的小提琴图。
# - output/01_preprocess/QC/*_spatial.pdf：QC 指标分箱后的空间分布图，用于检查组织区域/技术偏差。
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.config_loader import load_config
from lib.io_utils import read_scale_factors, read_visium_positions
from lib.plot_utils import save_spatial_scatter

cfg = load_config()
paths = cfg.PATHS
out_dir = paths["output"] / "01_preprocess"
(out_dir / "QC").mkdir(parents=True, exist_ok=True)

# scanpy.read_visium 会读取 filtered_feature_bc_matrix 和 spatial 图像/坐标。
# counts 层保留原始 UMI 计数，后续 normalize/log1p 会改 X，但 counts 仍可供 CNV/差异分析使用。
adata = sc.read_visium(paths["spaceranger"])
adata.var_names_make_unique()
adata.layers["counts"] = adata.X.copy()

# 不同 Space Ranger 版本的 tissue_positions 文件名和列名略有差异；
# 这里统一整理为 row/col/imagerow/imagecol，并保证 obs 和 obsm["spatial"] 都有坐标。
if "spatial" not in adata.obsm:
    pos = read_visium_positions(paths["spaceranger"] / "spatial")
    pos = pos.loc[adata.obs_names]
    adata.obsm["spatial"] = pos[["imagecol", "imagerow"]].to_numpy()
    adata.obs = adata.obs.drop(columns=adata.obs.columns.intersection(pos.columns), errors="ignore").join(pos, how="left")
else:
    pos = read_visium_positions(paths["spaceranger"] / "spatial")
    pos = pos.loc[adata.obs_names]
    adata.obs = adata.obs.drop(columns=adata.obs.columns.intersection(pos.columns), errors="ignore").join(pos, how="left")

adata.uns["spatial_scale_factors"] = read_scale_factors(paths["spaceranger"] / "spatial")
def _ravel(x):
    return np.asarray(x).ravel() if not hasattr(x, "A1") else x.A1


import numpy as np

# 计算常见 Visium QC：线粒体比例、总 UMI 数、检测到的基因数。
mito_sum = _ravel(adata[:, adata.var_names.str.startswith("MT-")].X.sum(axis=1))
total_sum = _ravel(adata.X.sum(axis=1))
adata.obs["Mito.percent"] = mito_sum / np.maximum(total_sum, 1) * 100
adata.obs["nCount_Spatial"] = total_sum
adata.obs["nFeature_Spatial"] = _ravel((adata.X > 0).sum(axis=1))

qc = adata.obs[["nCount_Spatial", "nFeature_Spatial", "Mito.percent"]]
qc.to_excel(out_dir / "QC" / "QCData.xlsx")

# 生成类似 Seurat VlnPlot 的 QC 小提琴图。
fig, axes = plt.subplots(1, 3, figsize=(9, 3))
for ax, key in zip(axes, ["nFeature_Spatial", "nCount_Spatial", "Mito.percent"]):
    ax.violinplot(qc[key].dropna())
    ax.set_title(key)
    ax.set_xticks([])
fig.tight_layout()
fig.savefig(out_dir / "QC" / "Vlnplot.pdf")
plt.close(fig)

# 把连续 QC 指标分成 5 个等级，在空间坐标上查看是否有明显组织/技术偏差。
for key in ["nFeature_Spatial", "nCount_Spatial", "Mito.percent"]:
    adata.obs[f"{key}_bin"] = pd.qcut(adata.obs[key].rank(method="first"), q=5, labels=False).astype(str)
    save_spatial_scatter(adata, f"{key}_bin", out_dir / "QC" / f"{key}_spatial.pdf", title=key)
    del adata.obs[f"{key}_bin"]

# h5ad 是后续 Python 复现流水线的主要对象格式。
adata.write_h5ad(paths["intermediate"] / "01_TumorST_preprocessed.h5ad")
