from __future__ import annotations

# 09：按空间区域 Mal / Bdy / nMal 做差异表达和富集分析。
# 每个区域都与“其它所有区域”比较，输出 DiffGenes_*.xlsx 和 pickle 格式结果。
# 输入：intermediate/05_TumorST_boundary_defined.h5ad，使用 counts 层重新 normalize/log1p。
# 输出：
# - output/09_diff_enrichment/DiffGenes_<Location>.xlsx：每个基因的 Diff、pvalue、Symbol、FDR。
# - intermediate/09_DiffGenes.pkl：所有差异表的 dict，供 10 火山图读取。
# - intermediate/09_Enrichment.pkl：gseapy GO/KEGG 富集结果，供后续复核。
# - output/09_diff_enrichment/enrichment_summary.json：每个 Location 成功产生的富集库摘要。
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.config_loader import load_config

cfg = load_config()
paths = cfg.PATHS
out_dir = paths["output"] / "09_diff_enrichment"
out_dir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(paths["intermediate"] / "05_TumorST_boundary_defined.h5ad")

# 优先从原始 counts 重新 normalize/log1p，避免前面步骤对 X 的临时改写影响差异分析。
if "counts" in adata.layers:
    adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
expr = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
expr = pd.DataFrame(expr, index=adata.obs_names, columns=adata.var_names)
locations = [x for x in ["Mal", "Bdy", "nMal"] if x in set(adata.obs["Location"].astype(str))]

# Welch t-test：对每个 Location，计算该区域 vs 其它区域的均值差和 p 值，再做 FDR 校正。
diff_tables = {}
for loc in locations:
    group = expr.loc[adata.obs["Location"].astype(str) == loc]
    other = expr.loc[adata.obs["Location"].astype(str) != loc]
    diff = group.mean(axis=0) - other.mean(axis=0)
    pvals = ttest_ind(group, other, axis=0, equal_var=False, nan_policy="omit").pvalue
    pvals = np.nan_to_num(pvals, nan=1.0)
    fdr = multipletests(pvals, method="fdr_bh")[1]
    tab = pd.DataFrame({"Diff": diff, "pvalue": pvals, "Symbol": expr.columns, "FDR": fdr}, index=expr.columns)
    diff_tables[loc] = tab
    tab.to_excel(out_dir / f"DiffGenes_{loc}.xlsx")

enrichment = {}
try:
    import gseapy as gp

    for loc, tab in diff_tables.items():
        # 富集只使用上调且显著的非免疫球蛋白/核糖体/线粒体等常见干扰基因。
        genes = tab.loc[
            (tab["Diff"] >= cfg.PARAMS["diff_logfc_cutoff"])
            & (tab["FDR"] <= cfg.PARAMS["diff_fdr_cutoff"])
            & ~tab["Symbol"].str.match(r"^IG[HJKL]|^RNA|^MT-|^RPS|^RPL"),
            "Symbol",
        ].tolist()
        enrichment[loc] = {}
        if genes:
            # Enrichr 查询可能受网络影响；外层 except 会把错误写成 note，不阻断差异表输出。
            enrichment[loc]["GO_Biological_Process_2021"] = gp.enrichr(gene_list=genes, gene_sets="GO_Biological_Process_2021", organism="human", outdir=None).results
            enrichment[loc]["KEGG_2021_Human"] = gp.enrichr(gene_list=genes, gene_sets="KEGG_2021_Human", organism="human", outdir=None).results
except ImportError:
    enrichment["note"] = "gseapy is not installed; enrichment was skipped."
except Exception as exc:
    enrichment["note"] = f"enrichment was skipped after gseapy error: {exc}"

with open(paths["intermediate"] / "09_DiffGenes.pkl", "wb") as fh:
    import pickle

    pickle.dump(diff_tables, fh)
with open(paths["intermediate"] / "09_Enrichment.pkl", "wb") as fh:
    import pickle

    pickle.dump(enrichment, fh)

# enrichment_summary.json 只保存每个区域成功得到哪些富集库，便于快速检查。
(out_dir / "enrichment_summary.json").write_text(json.dumps({k: list(v.keys()) if isinstance(v, dict) else v for k, v in enrichment.items()}, ensure_ascii=False, indent=2), encoding="utf-8")
