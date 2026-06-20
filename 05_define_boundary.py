from __future__ import annotations

# 05：根据 CNVLabel、CNV score、UMAP 位置和空间邻接关系划分 Mal / Bdy / nMal。
# 核心思想：先找高 CNV 的恶性种子，再沿空间邻居逐层扩展，靠近恶性核心的标为 Mal，
# 处在恶性区外缘和正常区之间的标为 Bdy，其余为 nMal。
# 输入：intermediate/04_TumorST_cnv_scored.h5ad，需要 obs 中已有 CNVLabel、cnv_score、seurat_clusters、NormalScore。
# 输出：
# - intermediate/05_TumorST_boundary_defined.h5ad：完整对象，obs["Location"] 为 Mal/Bdy/nMal，供 07/09/10 使用。
# - intermediate/05_TumorST_boundary_subset.h5ad：边界迭代过程中涉及的子集对象，便于复核。
# - output/05_boundary/CRC1_BoundaryDefine*.pdf：边界空间图和 H&E 叠加图。
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.boundary_helpers import cluster_update, compute_interspot_distances, find_neighbors, nbrs
from lib.config_loader import load_config
from lib.plot_utils import save_spatial_image_overlay, save_spatial_scatter

cfg = load_config()
paths = cfg.PATHS
out_dir = paths["output"] / "05_boundary"
out_dir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(paths["intermediate"] / "04_TumorST_cnv_scored.h5ad")
umap = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names, columns=["x", "y"])
pos = adata.obs[["row", "col", "imagerow", "imagecol"]].copy()

# 用 Visium 阵列坐标和像素坐标估计相邻 spot 的距离半径，再构建每个 spot 的空间邻居表。
dists = compute_interspot_distances(pos)
df_j = find_neighbors(pos, dists["radius"])

# 如果用户没有指定 malignant_cnv_labels，就把 CNV score 中位数最高的两个 CNVLabel 作为恶性候选。
mal_labels = cfg.PARAMS["malignant_cnv_labels"]
if mal_labels is None:
    med = adata.obs.groupby("CNVLabel")["cnv_score"].median().sort_values(ascending=False)
    mal_labels = med.index[:2].astype(str).tolist()
else:
    mal_labels = [str(x) for x in mal_labels]

# NormalScore 最高的聚类近似作为正常参考区。
mal_ids = adata.obs_names[adata.obs["CNVLabel"].astype(str).isin(mal_labels)].tolist()
normal_cluster = adata.obs.groupby("seurat_clusters")["NormalScore"].mean().sort_values(ascending=False).index[0]
normal_ids = adata.obs_names[adata.obs["seurat_clusters"].astype(str) == str(normal_cluster)].tolist()

# 找出“多数 spot 都属于恶性 CNVLabel”的空间聚类，作为恶性核心所在 cluster。
tab = pd.crosstab(adata.obs["CNVLabel"].astype(str), adata.obs["seurat_clusters"].astype(str))
cluster_ids = []
for cluster in sorted(adata.obs["seurat_clusters"].astype(str).unique()):
    if tab.reindex(mal_labels).fillna(0)[cluster].sum() > (adata.obs["seurat_clusters"].astype(str) == cluster).sum() * 0.5:
        cluster_ids.append(cluster)

# 在 UMAP 空间中计算恶性核心中心和正常中心；离恶性中心足够近的高 CNV spot 作为初始 Mal seed。
ci_mal = {}
sub_mal = {}
for cluster in cluster_ids:
    ids = [x for x in mal_ids if str(adata.obs.loc[x, "seurat_clusters"]) == cluster]
    if ids:
        sub_mal[cluster] = ids
        ci_mal[cluster] = umap.loc[ids].mean(axis=0)
ci_normal = umap.loc[normal_ids].mean(axis=0)

mal_seed = []
for cluster, ids in sub_mal.items():
    center = ci_mal[cluster]
    for cell in ids:
        rt = np.linalg.norm(umap.loc[cell] - center)
        rn = np.linalg.norm(umap.loc[cell] - ci_normal)
        if rt < (1 / 3) * rn:
            mal_seed.append(cell)

# 孤立的 Mal seed 外邻居用于第一圈边界判断。
lonely = [cell for cell in mal_seed if len([x for x in df_j.get(cell, []) if x in mal_seed]) == 0]
bdy_ids: list[str] = []
first_neighbors = nbrs(df_j, lonely, mal_seed + normal_ids + bdy_ids)
rows = []
for cell, neighbors in first_neighbors.items():
    cluster = str(adata.obs.loc[cell, "seurat_clusters"])
    for nid in neighbors:
        rt = np.linalg.norm(umap.loc[nid] - ci_mal.get(cluster, umap.loc[cell]))
        rn = np.linalg.norm(umap.loc[nid] - ci_normal)
        rows.append((nid, "Mal" if rt < (1 / 3) * rn else "Bdy"))
cluster_l = pd.DataFrame(rows, columns=["CellID", "Location"])
if not cluster_l.empty:
    cluster_l = cluster_l.groupby("CellID")["Location"].apply(lambda x: "Mal" if "Mal" in set(x) else "Bdy").reset_index()

mal_ids = mal_seed + cluster_l.loc[cluster_l["Location"] == "Mal", "CellID"].tolist()
bdy_ids = cluster_l.loc[cluster_l["Location"] == "Bdy", "CellID"].tolist()
mal_new = mal_ids.copy()
cluster_all = pd.concat(
    [
        pd.DataFrame({"CellID": normal_ids, "Location": "Normal"}),
        pd.DataFrame({"CellID": mal_ids, "Location": "Mal"}),
        pd.DataFrame({"CellID": bdy_ids, "Location": "Bdy"}),
    ],
    ignore_index=True,
)

# 迭代 1-6 圈：围绕上一圈新 Mal spot 找未标注邻居，并根据 UMAP 距离更新为 Maln 或 Bdy。
subset_ids = list(dict.fromkeys(normal_ids + mal_ids + bdy_ids))
for step in range(1, 7):
    if len(mal_new) < 3:
        break
    neighbors = nbrs(df_j, mal_new, mal_ids + normal_ids + bdy_ids)
    new_candidates = list(dict.fromkeys([x for values in neighbors.values() for x in values]))
    if len(new_candidates) < 3:
        break
    subset_ids = list(dict.fromkeys(new_candidates + mal_ids + normal_ids + bdy_ids))
    add = cluster_update(mal_new, df_j, umap, normal_ids, bdy_ids, mal_ids, step)
    cluster_all = pd.concat([cluster_all, add], ignore_index=True)
    label_new = cluster_all.drop_duplicates("CellID", keep="last").set_index("CellID")["Location"]
    mal_ids = [x for x in subset_ids if str(label_new.get(x, "")).startswith("Mal")]
    mal_new = [x for x in subset_ids if label_new.get(x, "") == f"Mal{step}"]
    bdy_ids = [x for x in subset_ids if label_new.get(x, "") == "Bdy"]

# 最终折叠为三类：Mal、Bdy、nMal。
mal_barcode = mal_ids
bdy_barcode = bdy_ids
normal_bdy = list(dict.fromkeys([cell for values in nbrs(df_j, mal_barcode, bdy_barcode + mal_barcode).values() for cell in values]))
nmal_barcode = [x for x in adata.obs_names if x not in set(mal_barcode + bdy_barcode + normal_bdy)]
loc = pd.Series(index=adata.obs_names, dtype=object)
loc.loc[mal_barcode] = "Mal"
loc.loc[bdy_barcode + normal_bdy] = "Bdy"
loc.loc[nmal_barcode] = "nMal"
adata.obs["Location"] = pd.Categorical(loc.loc[adata.obs_names], categories=["Mal", "Bdy", "nMal"])

# 输出边界图：一个纯空间散点图，一个叠加 H&E 背景图。
save_spatial_scatter(adata, "Location", out_dir / f"{cfg.SAMPLE_NAME}_BoundaryDefine.pdf", ["#CB181D", "#1f78b4", "#fdb462"])
save_spatial_image_overlay(
    adata,
    "Location",
    out_dir / f"{cfg.SAMPLE_NAME}_BoundaryDefine_HE_overlay.pdf",
    ["#CB181D", "#1f78b4", "#fdb462"],
    title="Boundary on H&E",
)
adata.write_h5ad(paths["intermediate"] / "05_TumorST_boundary_defined.h5ad")
adata[subset_ids].write_h5ad(paths["intermediate"] / "05_TumorST_boundary_subset.h5ad")
