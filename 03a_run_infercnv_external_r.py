from __future__ import annotations

# 03a：Python 流水线中的“外部 R inferCNV”路线。
# 它不直接在 Python 中推断 CNV，而是把 02 写好的矩阵和注释传给 resources/R/run_infercnv_external.R。
# 输入：intermediate/InferCNV/ 下的 matrix/cells/genes/CellAnnotation/reference_cluster，以及 resources/gencode_v38_gene_pos.txt。
# 输出：
# - intermediate/03_cnv_calls.tsv：三列表格 cell_ID、CNVLabel、cnv_score，是 04 的唯一必需输入。
# - output/03_infercnv/external_r/：R inferCNV 的中间对象、图和日志，主要用于排错/复核。
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.config_loader import load_config

cfg = load_config()
paths = cfg.PATHS

# 优先使用 00_config.py 中指定的 Rscript；如果不存在，则从 PATH 搜索 Rscript。
rscript = Path(cfg.R_EXE)
if not rscript.exists():
    found = shutil.which("Rscript")
    if found:
        rscript = Path(found)
    else:
        raise FileNotFoundError("Cannot find Rscript. Edit R_EXE in 00_config.py.")

# 传给 R 脚本的参数依次为：表达矩阵、基因名、细胞名、注释、参考 cluster、
# 基因位置文件、输出目录、CNV calls 输出路径、线程数和 CNV 聚类数。
cmd = [
    str(rscript),
    str(paths["resources"] / "R" / "run_infercnv_external.R"),
    str(paths["intermediate"] / "InferCNV" / "raw_counts_matrix.mtx"),
    str(paths["intermediate"] / "InferCNV" / "genes.tsv"),
    str(paths["intermediate"] / "InferCNV" / "cells.tsv"),
    str(paths["intermediate"] / "InferCNV" / "CellAnnotation.txt"),
    str(paths["intermediate"] / "InferCNV" / "reference_cluster.txt"),
    str(paths["resources"] / "gencode_v38_gene_pos.txt"),
    str(paths["output"] / "03_infercnv" / "external_r"),
    str(paths["intermediate"] / "03_cnv_calls.tsv"),
    str(cfg.PARAMS["infercnv_threads"]),
    str(cfg.PARAMS["cnv_k"]),
]

# check=True 表示 R 端只要报错，Python 也会以异常退出，避免产生半成品被误用。
subprocess.run(cmd, check=True)
