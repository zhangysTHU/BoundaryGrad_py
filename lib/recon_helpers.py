from __future__ import annotations

# helper：08_spatial_reconstruction.py 的表达重构。
# 输入 spot 原始 counts、反卷积比例和单细胞 signature；输出 gene x reconstructed_cell 矩阵。
import numpy as np
import pandas as pd


def get_feature_weight(filter_sig: pd.DataFrame) -> pd.DataFrame:
    # 对每个基因，把不同细胞类型 signature 归一化成权重，用于拆分 spot 的该基因 count。
    denom = filter_sig.sum(axis=1).replace(0, np.nan)
    return filter_sig.div(denom, axis=0).fillna(0)


def get_recon_matrix(counts: pd.DataFrame, obs: pd.DataFrame, sig_exp: pd.DataFrame, markers: dict[str, list[str]], decon: pd.DataFrame, locations: list[str]) -> pd.DataFrame:
    # 只重构指定 Location 中的 spot，例如默认的 Bdy。
    spot_ids = obs.index[obs["Location"].isin(locations)].tolist()
    print(f"[08] reconstructing {len(spot_ids)} spots for locations={locations}", flush=True)
    cluster_markers = [g for genes in markers.values() for g in genes if not g.startswith(("IGH", "IGJ", "IGK", "IGL", "RNA", "MT-", "RPS", "RPL"))]
    genes = counts.index.intersection(sig_exp.index).intersection(pd.Index(cluster_markers))
    sub_counts = counts.loc[genes, spot_ids]
    filter_sig = sig_exp.loc[genes]
    weights = get_feature_weight(filter_sig)
    matrices = []
    for idx, spot in enumerate(spot_ids, start=1):
        # 只为该 spot 中比例 > 0 的细胞类型生成 reconstructed columns。
        frac = decon.loc[spot]
        frac = frac[frac > 0]
        if frac.empty:
            continue
        sub = pd.DataFrame(0.0, index=genes, columns=[f"{ct}_{spot}" for ct in frac.index])
        for gene in genes:
            # deno 是该基因在当前 spot 的 signature 加权期望；
            # 若 deno 为 0，则无法按 signature 分配，退回把原 count 复制给各类型。
            sub_w = weights.loc[gene, frac.index]
            deno = float((sub_w * frac).sum())
            if deno == 0:
                sub.loc[gene] = sub_counts.loc[gene, spot]
            else:
                sub.loc[gene] = sub_counts.loc[gene, spot] / deno * sub_w.values
        matrices.append(sub)
        if idx == len(spot_ids) or idx % 50 == 0:
            print(f"[08] reconstructed {idx}/{len(spot_ids)} spots", flush=True)
    return pd.concat(matrices, axis=1) if matrices else pd.DataFrame(index=genes)
