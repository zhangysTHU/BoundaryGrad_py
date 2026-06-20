from __future__ import annotations

# helper：通用绘图函数。
# 主脚本只负责决定画什么变量，这里负责空间散点、H&E 叠加和 UMAP/embedding 散点。
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _palette_for_categories(cats: list[str], palette: list[str] | None = None) -> dict[str, str]:
    # 给离散类别循环分配颜色；类别数超过 palette 长度时会重复使用。
    colors = palette or plt.rcParams["axes.prop_cycle"].by_key()["color"]
    return {cat: colors[i % len(colors)] for i, cat in enumerate(cats)}


def save_spatial_scatter(adata, color_key: str, path: Path, palette: list[str] | None = None, title: str | None = None) -> None:
    # 不叠加组织图片，只用 adata.obsm["spatial"] 画 spot 空间分布。
    coords = adata.obsm["spatial"]
    values = adata.obs[color_key].astype(str)
    cats = list(pd.Categorical(values).categories)
    color_map = _palette_for_categories(cats, palette)
    fig, ax = plt.subplots(figsize=(7, 7))
    for cat in cats:
        mask = values == cat
        ax.scatter(coords[mask, 0], coords[mask, 1], s=12, c=color_map[cat], label=cat, alpha=0.8, linewidths=0)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_title(title or color_key)
    ax.axis("off")
    ax.legend(markerscale=2, frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def save_spatial_image_overlay(
    adata,
    color_key: str,
    path: Path,
    palette: list[str] | None = None,
    title: str | None = None,
    image_key: str = "hires",
    spot_size: float | None = None,
    alpha: float = 0.75,
) -> None:
    # 如果 AnnData 中保留了 scanpy.read_visium 读入的 H&E 图，则把 spot 画到图上；
    # 否则退回普通空间散点图。
    spatial = adata.uns.get("spatial")
    if not spatial:
        save_spatial_scatter(adata, color_key, path, palette, title)
        return

    library_id = next(iter(spatial))
    library = spatial[library_id]
    images = library.get("images", {})
    if image_key not in images:
        image_key = next(iter(images))
    image = images[image_key]
    scalefactors = library.get("scalefactors", {})
    scale = scalefactors.get(f"tissue_{image_key}_scalef", 1.0)
    coords = adata.obsm["spatial"] * scale

    values = adata.obs[color_key]
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(image)
    if pd.api.types.is_numeric_dtype(values):
        scatter = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            s=spot_size or 8,
            c=values.astype(float),
            cmap="viridis",
            alpha=alpha,
            linewidths=0,
        )
        fig.colorbar(scatter, ax=ax, fraction=0.035, pad=0.02, label=color_key)
    else:
        values = values.astype(str)
        cats = list(pd.Categorical(values).categories)
        color_map = _palette_for_categories(cats, palette)
        for cat in cats:
            mask = values == cat
            ax.scatter(
                coords[mask, 0],
                coords[mask, 1],
                s=spot_size or 8,
                c=color_map[cat],
                label=cat,
                alpha=alpha,
                linewidths=0,
            )
        ax.legend(markerscale=2, frameon=False, loc="best")
    ax.set_title(title or color_key)
    ax.set_xlim(0, image.shape[1])
    ax.set_ylim(image.shape[0], 0)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def save_embedding_scatter(adata, basis: str, color_key: str, path: Path, palette: list[str] | None = None, title: str | None = None) -> None:
    # 用于 UMAP 等低维坐标图，basis 通常是 "X_umap"。
    coords = adata.obsm[basis]
    values = adata.obs[color_key].astype(str)
    cats = list(pd.Categorical(values).categories)
    color_map = _palette_for_categories(cats, palette)
    fig, ax = plt.subplots(figsize=(7, 7))
    for cat in cats:
        mask = values == cat
        ax.scatter(coords[mask, 0], coords[mask, 1], s=10, c=color_map[cat], label=cat, alpha=0.8, linewidths=0)
    ax.set_title(title or color_key)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.legend(markerscale=2, frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
