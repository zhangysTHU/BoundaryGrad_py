from __future__ import annotations

# helper：文件/目录读写工具，主要解决 10x Visium 不同版本文件名和 inferCNV 输入格式问题。
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


def ensure_dir(path: Path) -> Path:
    # 创建目录并返回 Path，便于一行里初始化输出目录。
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_visium_positions(spatial_dir: Path) -> pd.DataFrame:
    # Space Ranger 老版本常用 tissue_positions_list.csv，新版本常用 tissue_positions.csv；
    # 两者列名/是否有 header 不完全一致，这里统一成后续脚本使用的列名。
    candidates = [
        spatial_dir / "tissue_positions.csv",
        spatial_dir / "tissue_positions_list.csv",
    ]
    for path in candidates:
        if path.exists():
            if path.name == "tissue_positions.csv":
                df = pd.read_csv(path)
                if "barcode" not in df.columns:
                    df = pd.read_csv(path, header=None)
            else:
                df = pd.read_csv(path, header=None)
            break
    else:
        raise FileNotFoundError(f"No tissue_positions file found in {spatial_dir}")

    if list(df.columns[:6]) != ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]:
        df = df.iloc[:, :6]
        df.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    df = df.set_index("barcode")
    df = df.rename(
        columns={
            "array_row": "row",
            "array_col": "col",
            "pxl_row_in_fullres": "imagerow",
            "pxl_col_in_fullres": "imagecol",
        }
    )
    return df


def read_scale_factors(spatial_dir: Path) -> dict:
    # 读取 hires/lowres 图像缩放因子，用于把 full-res 坐标映射到绘图图像。
    path = spatial_dir / "scalefactors_json.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_infercnv_inputs(adata, out_dir: Path, group_key: str = "seurat_clusters") -> None:
    # inferCNV R 端需要 gene x cell 的 MatrixMarket 矩阵、基因名、细胞名和 cell 注释。
    ensure_dir(out_dir)
    matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    matrix = matrix.T
    if not sparse.issparse(matrix):
        matrix = sparse.csr_matrix(matrix)
    mmwrite(out_dir / "raw_counts_matrix.mtx", matrix)
    pd.Series(adata.var_names, name="gene").to_csv(out_dir / "genes.tsv", sep="\t", index=False)
    pd.Series(adata.obs_names, name="cell").to_csv(out_dir / "cells.tsv", sep="\t", index=False)
    pd.DataFrame({"CellID": adata.obs_names, "DefineTypes": adata.obs[group_key].astype(str).values}).to_csv(
        out_dir / "CellAnnotation.txt", sep="\t", index=False, header=False
    )


def save_cnv_calls(path: Path, calls: pd.DataFrame) -> None:
    # 统一 03a/03b 的 CNV 输出格式，04 只认这三列。
    required = {"cell_ID", "CNVLabel", "cnv_score"}
    missing = required.difference(calls.columns)
    if missing:
        raise ValueError(f"CNV calls missing columns: {sorted(missing)}")
    calls.loc[:, ["cell_ID", "CNVLabel", "cnv_score"]].to_csv(path, sep="\t", index=False)
