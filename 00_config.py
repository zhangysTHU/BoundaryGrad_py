from __future__ import annotations

import os
from pathlib import Path


# 本文件集中管理 Python 复现流水线的路径、样本名、路线名和算法参数。
# 其它脚本都通过 lib.config_loader.load_config() 动态加载这里的变量，
# 因此改输入目录、样本名、聚类分辨率、CNV/反卷积参数时，优先看这里。
# 注意：本脚本本身不产生生物学输出；它会确保 input/intermediate/output/resources/lib 等目录存在。
SCRIPT_DIR = Path(__file__).resolve().parent

# SAMPLE_NAME 对应 input/ 下的样本目录，也作为默认输出目录名。
# 可通过环境变量 COTTRAZM_SAMPLE_NAME 切换样本。
SAMPLE_NAME = os.environ.get("COTTRAZM_SAMPLE_NAME", "CRC1")

# ROUTE_SUFFIX 用于同一样本的不同路线输出隔离。
# 例如 run_all_03a.ps1 会设置 COTTRAZM_ROUTE_SUFFIX=03a，
# 输出目录和中间目录即变为 output/CRC1_03a/ 与 intermediate/CRC1_03a/。
ROUTE_SUFFIX = os.environ.get("COTTRAZM_ROUTE_SUFFIX", "").strip()
RUN_NAME = f"{SAMPLE_NAME}_{ROUTE_SUFFIX}" if ROUTE_SUFFIX else SAMPLE_NAME

# 所有路径都相对 scripts_format_python 目录组织：
# input/<样本名>/ 存放输入，intermediate/<运行名>/ 存放跨步骤中间结果，
# output/<运行名>/ 存放图表和表格。
# intermediate/ 文件通常是下一步真正读取的文件；output/ 文件主要给人工检查和论文图表使用。
INPUT_ROOT = SCRIPT_DIR / "input"
OUTPUT_ROOT = SCRIPT_DIR / "output"
INTERMEDIATE_ROOT = SCRIPT_DIR / "intermediate"
SAMPLE_INPUT_DIR = INPUT_ROOT / SAMPLE_NAME
RUN_OUTPUT_DIR = OUTPUT_ROOT / RUN_NAME
RUN_INTERMEDIATE_DIR = INTERMEDIATE_ROOT / RUN_NAME

PATHS = {
    "input": SAMPLE_INPUT_DIR,
    "spaceranger": SAMPLE_INPUT_DIR / "spaceranger_outs",
    "single_cell": SAMPLE_INPUT_DIR / "single_cell",
    "intermediate": RUN_INTERMEDIATE_DIR,
    "output": RUN_OUTPUT_DIR,
    "resources": SCRIPT_DIR / "resources",
    "lib": SCRIPT_DIR / "lib",
}

for path in [INPUT_ROOT, OUTPUT_ROOT, INTERMEDIATE_ROOT, *PATHS.values()]:
    path.mkdir(parents=True, exist_ok=True)

# 03a 路线会调用外部 R 脚本运行 inferCNV；如果只用 03b infercnvpy，可不改这里。
R_EXE = r"C:\Program Files\R\R-4.3.1\bin\x64\Rscript.exe"

PARAMS = {
    # 02 形态校正聚类的 Leiden 分辨率，越大通常 cluster 越多。
    "cluster_resolution": 1.5,
    # 03/04 中与 CNV 推断和分组有关的参数。
    "infercnv_assay": "Spatial",
    "infercnv_threads": 30,
    "cnv_k": 8,
    # 若为 None，05 会自动把 CNV score 中位数最高的两个 CNVLabel 视作恶性标签。
    "malignant_cnv_labels": None,
    # 07 反卷积中用来强制补入或辅助判断的细胞类型名称，需与 sig_exp.csv 列名一致。
    "decon_malignant_cluster": "Malignant epithelial cells",
    "decon_tissue_cluster": "Epithelial cells",
    "decon_stromal_cluster": "Fibroblast cells",
    # 08 默认只对边界区 Bdy 做空间重构。
    "recon_locations": ["Bdy"],
    # 09 差异分析和 10 火山图的阈值。
    "diff_logfc_cutoff": 0.25,
    "diff_fdr_cutoff": 0.05,
    "volcano_p_cutoff_log10": 2,
    "volcano_label_n": 10,
    "pie_scale": 0.4,
    "scatterpie_alpha": 0.8,
    "pie_border_color": "grey",
    # 11 LSGI 细胞组分梯度分析和可视化参数。
    "lsgi_n_grids_scale": 10,
    "lsgi_n_cells_per_meta": 50,
    "lsgi_r_squared_thresh": 0.3,
    "lsgi_minimum_fctr": 3,
    "lsgi_arrow_length_scale": 1.4,
    "lsgi_arrow_linewidth": 1.0,
    "lsgi_arrow_head_cm": 0.20,
    "lsgi_arrow_closed": True,
    "lsgi_image_key": "lowres",
}

# 空间聚类和 CNV 标签图使用的离散调色板。
CLUSTER_COLS = [
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
    "#33A02C", "#B2DF8A", "#55B1B1", "#8DD3C7", "#A6761D",
    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00",
]
