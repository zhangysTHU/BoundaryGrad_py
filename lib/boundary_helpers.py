from __future__ import annotations

# helper：05_define_boundary.py 的空间邻居和边界扩展算法。
# 这里不读取文件，只接收坐标、邻居表和 UMAP 坐标，返回可复用的中间结果。
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist


def compute_interspot_distances(position: pd.DataFrame, scale_factor: float = 1.05) -> dict:
    # 用阵列 row/col 与图像像素 imagerow/imagecol 的线性比例估计相邻 spot 的像素距离。
    xcoef = np.polyfit(position["col"].astype(float), position["imagecol"].astype(float), 1)[0]
    ycoef = np.polyfit(position["row"].astype(float), position["imagerow"].astype(float), 1)[0]
    return {"xdist": xcoef, "ydist": ycoef, "radius": (abs(xcoef) + abs(ycoef)) * scale_factor}


def find_neighbors(position: pd.DataFrame, radius: float, method: str = "cityblock") -> dict[str, list[str]]:
    # 在图像坐标中计算 spot 间距离；cityblock/manhattan 更贴近 Visium 六边形阵列的近邻搜索。
    coords = position.loc[:, ["imagecol", "imagerow"]].to_numpy(float)
    dist = cdist(coords, coords, metric=method)
    names = position.index.astype(str).to_numpy()
    return {name: list(names[(dist[i] <= radius) & (dist[i] > 0)]) for i, name in enumerate(names)}


def nbrs(df_j: dict[str, list[str]], mal_cell_ids: list[str], raw_ids: list[str]) -> dict[str, list[str]]:
    # 对一组恶性边缘 spot，找出尚未属于 raw_ids 的空间邻居，作为下一轮候选。
    raw = set(raw_ids)
    return {cell: [x for x in df_j.get(cell, []) if x not in raw] for cell in mal_cell_ids}


def cluster_update(
    mal_cell_new: list[str],
    df_j: dict[str, list[str]],
    umap: pd.DataFrame,
    normal_ids: list[str],
    bdy_ids: list[str],
    mal_ids: list[str],
    step: int,
) -> pd.DataFrame:
    # 根据当前新恶性层 Mal{step-1} 的邻居，决定未标注候选点是继续成为 Mal{step} 还是 Bdy。
    # 判断标准结合 UMAP 上到恶性中心/边界中心的距离，模拟原 Cottrazm 的逐圈扩展逻辑。
    rows = []
    normal_set, bdy_set, mal_set = set(normal_ids), set(bdy_ids), set(mal_ids)
    for cell in mal_cell_new:
        ncell = df_j.get(cell, [])
        if not ncell:
            continue
        cpos = umap.loc[cell]
        nmal = [x for x in ncell if x in mal_set]
        if nmal:
            ci_mal = pd.concat([umap.loc[nmal], cpos.to_frame().T]).mean(axis=0)
            r_mal = [np.linalg.norm(umap.loc[x] - ci_mal) for x in nmal + [cell]]
        else:
            ci_mal = cpos
            r_mal = [0.0]
        nbdy = [x for x in ncell if x in bdy_set]
        candidates = [x for x in ncell if x not in mal_set and x not in bdy_set and x not in normal_set]
        if not candidates:
            continue
        if len(nbdy) <= 1:
            for cid in candidates:
                p1 = np.linalg.norm(umap.loc[cid] - ci_mal)
                rows.append((cid, "Mal" if p1 <= 0.8 * max(r_mal) else "Bdy"))
        else:
            ci_bdy = umap.loc[nbdy].mean(axis=0)
            r_bdy = [np.linalg.norm(umap.loc[x] - ci_bdy) for x in nbdy]
            for cid in candidates:
                p1 = np.linalg.norm(umap.loc[cid] - ci_mal)
                p2 = np.linalg.norm(umap.loc[cid] - ci_bdy)
                rows.append((cid, f"Mal{step}" if p1 <= 0.8 * max(r_mal) and p2 > max(r_bdy) else "Bdy"))
    if not rows:
        return pd.DataFrame(columns=["CellID", "Location"])
    df = pd.DataFrame(rows, columns=["CellID", "Location"])
    collapsed = df.groupby("CellID")["Location"].apply(lambda x: f"Mal{step}" if any(v.startswith("Mal") for v in x) else "Bdy")
    return collapsed.reset_index()
