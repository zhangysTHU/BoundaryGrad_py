from __future__ import annotations

# 06：准备 07 反卷积需要的单细胞参考。
# 推荐输入：sig_exp.csv（细胞类型 signature 表达矩阵）和 clustermarkers_list.json（各细胞类型 marker）。
# 备用输入：single_cell_reference.h5ad + marker JSON，脚本会按 define_types 指定的 obs 列计算平均表达。
# 输入格式：
# - input/single_cell/sig_exp.csv：行是 gene，列是 cell type，值是该类型平均表达/signature。
# - input/single_cell/clustermarkers_list.json：键是 cell type，值是 marker gene 列表。
# 输出：
# - intermediate/06_sig_exp.csv：07/08 读取的 signature matrix。
# - intermediate/06_clustermarkers_list.json：07 用于 marker enrichment，08 用于选择重构基因。
# - intermediate/06_reference.pkl：同一份参考的 pickle 打包，便于调试或交互读取。
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.config_loader import load_config

cfg = load_config()
paths = cfg.PATHS

# define_types 决定从 h5ad 单细胞参考中按哪一列分组，默认是 Majortypes。
sig_csv = paths["single_cell"] / "sig_exp.csv"
markers_json = paths["single_cell"] / "clustermarkers_list.json"
sc_h5ad = paths["single_cell"] / "single_cell_reference.h5ad"
define_types = "Majortypes"
define_file = paths["single_cell"] / "define_types.txt"
if define_file.exists():
    define_types = define_file.read_text(encoding="utf-8").strip().splitlines()[0]

# 情况 1：已经提供整理好的 signature matrix 和 marker 列表，直接读取。
if sig_csv.exists() and markers_json.exists():
    sig_exp = pd.read_csv(sig_csv, index_col=0)
    markers = json.loads(markers_json.read_text(encoding="utf-8"))
# 情况 2：只有单细胞 h5ad，则用 marker 涉及的基因，按细胞类型求平均表达得到 signature。
elif sc_h5ad.exists() and markers_json.exists():
    sc_adata = sc.read_h5ad(sc_h5ad)
    markers = json.loads(markers_json.read_text(encoding="utf-8"))
    genes = [g for values in markers.values() for g in values if g in sc_adata.var_names]
    expr = sc_adata[:, sorted(set(genes))].X
    expr = expr.toarray() if sparse.issparse(expr) else np.asarray(expr)
    expr = np.power(2, expr) - 1
    expr_df = pd.DataFrame(expr, index=sc_adata.obs_names, columns=sorted(set(genes)))
    sig_exp = expr_df.groupby(sc_adata.obs[define_types]).mean().T
else:
    raise FileNotFoundError(
        "Provide input/single_cell/sig_exp.csv and clustermarkers_list.json, "
        "or single_cell_reference.h5ad plus clustermarkers_list.json."
    )

# 输出三种格式：CSV/JSON 供后续脚本直接读，pkl 方便一次性复用。
sig_exp.to_csv(paths["intermediate"] / "06_sig_exp.csv")
(paths["intermediate"] / "06_clustermarkers_list.json").write_text(json.dumps(markers, ensure_ascii=False, indent=2), encoding="utf-8")
with open(paths["intermediate"] / "06_reference.pkl", "wb") as fh:
    pickle.dump({"sig_exp": sig_exp, "markers": markers}, fh)
