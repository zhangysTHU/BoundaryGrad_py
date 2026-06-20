from __future__ import annotations

# helper：07_spatial_deconvolution.py 的 marker 富集和约束最小二乘反卷积。
import numpy as np
import pandas as pd
from scipy.optimize import minimize


def get_enrich_matrix(filter_sig: pd.DataFrame, markers: dict[str, list[str]]) -> pd.DataFrame:
    # 构造 gene x celltype 的 0/1 marker 矩阵：1 表示该基因是该细胞类型 marker。
    mat = pd.DataFrame(0, index=filter_sig.index, columns=filter_sig.columns, dtype=int)
    for cluster in mat.columns:
        genes = [g for g in markers.get(cluster, []) if g in mat.index]
        mat.loc[genes, cluster] = 1
    return mat


def enrich_analysis(filter_log_expr: pd.DataFrame, enrich_matrix: pd.DataFrame) -> pd.DataFrame:
    # 类似 z-score 的 marker enrichment：
    # 对每个 spot，比较某细胞类型 marker 的平均 fold-change 是否高于全基因背景。
    mean_gene_expr = np.log2(np.mean(np.power(2, filter_log_expr) - 1, axis=1) + 1)
    gene_fold = filter_log_expr.sub(mean_gene_expr, axis=0)
    cell_mean = gene_fold.mean(axis=0)
    cell_sd = gene_fold.std(axis=0).replace(0, np.nan)
    enrichment = pd.DataFrame(index=enrich_matrix.columns, columns=filter_log_expr.columns, dtype=float)
    for cluster in enrich_matrix.columns:
        genes = enrich_matrix.index[enrich_matrix[cluster] == 1]
        sig_mean = gene_fold.loc[genes].mean(axis=0)
        enrichment.loc[cluster] = (sig_mean - cell_mean) * np.sqrt(len(genes)) / cell_sd
    return enrichment.fillna(0)


def solve_nnls_sum_to_one(signature: np.ndarray, bulk: np.ndarray) -> np.ndarray:
    # 求解非负、和为 1 的最小二乘：signature @ proportion 尽量接近 spot 表达。
    n = signature.shape[1]
    x0 = np.ones(n) / n
    bounds = [(0, None)] * n
    cons = {"type": "eq", "fun": lambda x: np.sum(x) - 1}
    result = minimize(lambda x: np.sum((signature @ x - bulk) ** 2), x0, method="SLSQP", bounds=bounds, constraints=cons)
    if not result.success:
        return x0
    return result.x / max(result.x.sum(), 1e-12)


def deconvolve_topic(expr: pd.DataFrame, sig: pd.DataFrame) -> pd.DataFrame:
    # 对同一个 topic 内的每个 spot 分别求比例；返回 celltype x spot 的矩阵。
    genes = expr.index.intersection(sig.index)
    expr = expr.loc[genes]
    sig = sig.loc[genes]
    res = []
    S = sig.to_numpy(float)
    for cell in expr.columns:
        res.append(solve_nnls_sum_to_one(S, expr[cell].to_numpy(float)))
    return pd.DataFrame(np.vstack(res).T, index=sig.columns, columns=expr.columns)
