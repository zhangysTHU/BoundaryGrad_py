from __future__ import annotations

# 10：汇总可视化。
# 输出反卷积细胞比例柱状图、空间 pie 图，以及每个 Location 的差异火山图。
# 输入：
# - intermediate/05_TumorST_boundary_defined.h5ad：Location 和空间坐标。
# - intermediate/07_DeconData.tsv：反卷积比例，用于 bar/pie 图。
# - intermediate/09_DiffGenes.pkl：差异表，用于火山图。
# 输出：
# - output/10_plots/DeconBarplot.pdf：各 Location 的细胞组分百分比。
# - output/10_plots/DeconPieplot.pdf：每个 spot 的细胞组分 pie，叠加低分辨率 H&E。
# - output/10_plots/DiffVolcano_<Location>.pdf：各 Location 差异基因火山图。
import pickle
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib.patches import Wedge

try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.config_loader import load_config

cfg = load_config()
paths = cfg.PATHS
out_dir = paths["output"] / "10_plots"
out_dir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(paths["intermediate"] / "05_TumorST_boundary_defined.h5ad")
decon = pd.read_csv(paths["intermediate"] / "07_DeconData.tsv", sep="\t").set_index("cell_ID")
with open(paths["intermediate"] / "09_DiffGenes.pkl", "rb") as fh:
    diff_tables = pickle.load(fh)

# 1) 各空间区域中细胞类型比例的堆叠柱状图。
plot_cols = decon.columns.tolist()
meta = decon.join(adata.obs["Location"], how="left")
bar = meta.groupby("Location")[plot_cols].sum()
bar = bar.div(bar.sum(axis=1), axis=0) * 100
fig, ax = plt.subplots(figsize=(7, 5))
bottom = np.zeros(len(bar))
for ct in plot_cols:
    ax.bar(bar.index.astype(str), bar[ct], bottom=bottom, label=ct)
    bottom += bar[ct].values
ax.set_ylabel("Percent")
ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
fig.tight_layout()
fig.savefig(out_dir / "DeconBarplot.pdf")
plt.close(fig)

# 2) 在低分辨率 H&E 图上，用 pie chart 显示每个 spot 的反卷积组成。
img_path = paths["spaceranger"] / "spatial" / "tissue_lowres_image.png"
if img_path.exists():
    img = mpimg.imread(img_path)
    scale = adata.uns.get("spatial_scale_factors", {}).get("tissue_lowres_scalef", 1.0)
    coords = pd.DataFrame(adata.obsm["spatial"], index=adata.obs_names, columns=["x", "y"]) * scale
    pie_data = coords.join(decon, how="inner")
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(img)
    colors = sns.color_palette("tab20", n_colors=len(plot_cols))
    color_map = dict(zip(plot_cols, colors))
    radius = cfg.PARAMS["pie_scale"] * 10
    for _, row in pie_data.iterrows():
        # 每个 Wedge 是一个细胞类型的扇区，所有扇区合起来就是一个 spot 的组成。
        start = 0.0
        total = max(row[plot_cols].sum(), 1e-12)
        for ct in plot_cols:
            frac = row[ct] / total
            if frac <= 0:
                continue
            wedge = Wedge((row["x"], row["y"]), radius, start * 360, (start + frac) * 360, facecolor=color_map[ct], alpha=cfg.PARAMS["scatterpie_alpha"], edgecolor=cfg.PARAMS["pie_border_color"], linewidth=0.1)
            ax.add_patch(wedge)
            start += frac
    ax.set_xlim(0, img.shape[1])
    ax.set_ylim(img.shape[0], 0)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_dir / "DeconPieplot.pdf")
    plt.close(fig)

# 3) 对 09 的每个差异表画火山图，显著上调区域为红色，其它显著下调为蓝色。
for loc, tab in diff_tables.items():
    tab = tab.loc[~tab["Symbol"].str.match(r"^IG[HJKL]|^RNA|^MT-|^RPS|^RPL")].copy()
    tab["color"] = np.where(
        (-np.log10(tab["FDR"].clip(lower=1e-300)) < cfg.PARAMS["volcano_p_cutoff_log10"]) | (tab["Diff"].abs() <= cfg.PARAMS["diff_logfc_cutoff"]),
        "grey",
        np.where(tab["Diff"] > 0, loc, "Other"),
    )
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.scatterplot(data=tab, x="Diff", y=-np.log10(tab["pvalue"].clip(lower=1e-300)), hue="color", palette={loc: "red", "Other": "blue", "grey": "grey"}, s=10, ax=ax, linewidth=0)
    selected = tab[tab["color"] != "grey"].assign(absdiff=lambda x: x["Diff"].abs()).nlargest(cfg.PARAMS["volcano_label_n"], "absdiff")
    texts = [ax.text(row["Diff"], -np.log10(max(row["pvalue"], 1e-300)), row["Symbol"], fontsize=8) for _, row in selected.iterrows()]
    if adjust_text is not None:
        adjust_text(texts, ax=ax)
    ax.axvline(cfg.PARAMS["diff_logfc_cutoff"], ls="--", c="black", lw=0.8)
    ax.axvline(-cfg.PARAMS["diff_logfc_cutoff"], ls="--", c="black", lw=0.8)
    ax.axhline(cfg.PARAMS["volcano_p_cutoff_log10"], ls="--", c="black", lw=0.8)
    ax.set_xlabel("log2(FoldChange)")
    ax.set_ylabel("-log10(pvalue)")
    fig.tight_layout()
    fig.savefig(out_dir / f"DiffVolcano_{loc}.pdf")
    plt.close(fig)
