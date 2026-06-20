from __future__ import annotations

# 04：把 03 得到的 CNV calls 合并回 AnnData，并生成 CNV 标签/分数可视化。
# 输入统一为 intermediate/03_cnv_calls.tsv，因此 03a 和 03b 两条路线都可接到这里。
# 输入：
# - intermediate/02_TumorST_clustered.h5ad：spot 聚类、UMAP、空间坐标。
# - intermediate/03_cnv_calls.tsv：cell_ID、CNVLabel、cnv_score。
# 输出：
# - intermediate/04_TumorST_cnv_scored.h5ad：在 obs 中加入 CNVLabel/cnv_score，供 05 边界划分使用。
# - output/04_cnv_score/*.pdf：CNVLabel 空间图、H&E 叠加图、UMAP 图和 score 分布图。
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.config_loader import load_config
from lib.plot_utils import save_embedding_scatter, save_spatial_image_overlay, save_spatial_scatter

cfg = load_config()
paths = cfg.PATHS
out_dir = paths["output"] / "04_cnv_score"
out_dir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(paths["intermediate"] / "02_TumorST_clustered.h5ad")
calls = pd.read_csv(paths["intermediate"] / "03_cnv_calls.tsv", sep="\t")

# 按 cell_ID 对齐，缺失的 spot 默认 Normal / 0 分，避免下游因为少量缺失直接失败。
calls = calls.set_index("cell_ID").reindex(adata.obs_names)
adata.obs["CNVLabel"] = calls["CNVLabel"].fillna("Normal").astype(str)
adata.obs["cnv_score"] = calls["cnv_score"].fillna(0).astype(float)

# 生成空间图、H&E 叠加图、UMAP 图和 CNV score 小提琴/箱线图。
save_spatial_scatter(adata, "CNVLabel", out_dir / f"{cfg.SAMPLE_NAME}_cnv_label.pdf", cfg.CLUSTER_COLS)
save_spatial_image_overlay(adata, "CNVLabel", out_dir / f"{cfg.SAMPLE_NAME}_cnv_label_HE_overlay.pdf", cfg.CLUSTER_COLS, title="CNVLabel on H&E")
save_spatial_image_overlay(adata, "cnv_score", out_dir / f"{cfg.SAMPLE_NAME}_cnv_score_HE_overlay.pdf", title="CNV score on H&E")
save_embedding_scatter(adata, "X_umap", "CNVLabel", out_dir / f"{cfg.SAMPLE_NAME}_reduction_cnvlabel.pdf", cfg.CLUSTER_COLS)

fig, ax = plt.subplots(figsize=(7, 4))
plot_obs = adata.obs.copy()
label_order = plot_obs.groupby("CNVLabel", observed=False)["cnv_score"].median().sort_values(ascending=False).index.tolist()
sns.violinplot(data=plot_obs, x="CNVLabel", y="cnv_score", order=label_order, ax=ax, inner=None, alpha=0.5)
sns.boxplot(
    data=plot_obs,
    x="CNVLabel",
    y="cnv_score",
    order=label_order,
    ax=ax,
    width=0.35,
    showcaps=True,
    boxprops={"facecolor": "white", "zorder": 2},
    fliersize=0.5,
)
ax.set_title("CNV Scores")
ax.set_xlabel("CNVLabel")
ax.set_ylabel("CNV_scores")
ax.tick_params(axis="x", rotation=45)
fig.tight_layout()
fig.savefig(out_dir / f"{cfg.SAMPLE_NAME}_cnv_observation_vlnplot.pdf")
plt.close(fig)

# 带 CNVLabel/cnv_score 的对象进入 05 边界定义。
adata.write_h5ad(paths["intermediate"] / "04_TumorST_cnv_scored.h5ad")
