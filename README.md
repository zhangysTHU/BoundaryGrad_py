# BoundaryGrad Python workflow

本仓库是 Cottrazm/BoundaryGrad 的 Python 顺序脚本工作流。默认样本名为 `CRC1`，所有输入、输出和中间文件都按样本名分目录保存，便于后续加入多个样本。

## 目录放置规则

默认目录如下：

```text
BoundaryGrad_py/
  input/
    CRC1/
      spaceranger_outs/
      single_cell/
  intermediate/
    CRC1/
  output/
    CRC1/
```

`input/CRC1/` 是样本输入目录。`intermediate/CRC1/` 是脚本之间传递的中间文件目录，一般不需要手动修改。`output/CRC1/` 是图、表和人工检查结果目录。

为避免上传大体积测序数据、患者相关数据和可重复生成的运行产物，Git 仓库不跟踪 `input/`、`intermediate/` 与 `output/` 中的内容。克隆后请按下述结构自行放置输入；运行流程会自动创建中间和输出目录。

03a 外部 R inferCNV 路线会使用同一个输入样本，但输出和中间文件放在带后缀的运行目录：

```text
intermediate/CRC1_03a/
output/CRC1_03a/
```

## 空间转录组输入

Space Ranger 结果必须放在：

```text
input/CRC1/spaceranger_outs/
```

该目录至少需要包含以下内容之一：

```text
filtered_feature_bc_matrix.h5
```

或矩阵目录：

```text
filtered_feature_bc_matrix/
  matrix.mtx.gz
  features.tsv.gz
  barcodes.tsv.gz
```

同时必须包含空间图像和坐标目录：

```text
spatial/
  tissue_positions_list.csv
  scalefactors_json.json
  tissue_lowres_image.png
  tissue_hires_image.png
```

如果 Space Ranger 版本输出的是 `tissue_positions.csv`，请确认脚本中的读取函数支持该文件，或将其转换为当前示例使用的 `tissue_positions_list.csv` 格式。

## 单细胞参考输入

单细胞参考文件必须放在：

```text
input/CRC1/single_cell/
```

推荐直接提供整理好的 signature 和 marker 文件：

```text
sig_exp.csv
clustermarkers_list.json
```

`sig_exp.csv` 要求：

- 行名是基因名。
- 列名是细胞类型名。
- 数值是该细胞类型的平均表达或 signature 表达。

`clustermarkers_list.json` 要求：

- 顶层是 JSON object。
- 每个键是细胞类型名。
- 每个值是 marker gene 字符串数组。

也可以提供：

```text
single_cell_reference.h5ad
clustermarkers_list.json
define_types.txt
```

其中 `single_cell_reference.h5ad` 应包含单细胞表达矩阵和细胞类型注释；`define_types.txt` 用于指定或映射细胞类型列名。

## 运行方式

请先准备包含本流程所需科学计算包的 Python 环境。PowerShell 脚本默认调用 PATH 中的 `python`；也可以显式指定解释器：

```powershell
.\run_all.ps1 -Python "C:\path\to\python.exe"
```

03a 路线还需要 R 与 `infercnv` 等依赖。默认调用 PATH 中的 `Rscript`，也可设置环境变量：

```powershell
$env:COTTRAZM_RSCRIPT = "C:\path\to\Rscript.exe"
```

默认运行 03b infercnvpy 路线：

```powershell
.\run_all.ps1
```

运行 03a 外部 R inferCNV 路线：

```powershell
.\run_all_03a.ps1
```

指定其他样本名：

```powershell
.\run_all.ps1 -SampleName CRC2
.\run_all_03a.ps1 -SampleName CRC2
```

此时输入应放在 `input/CRC2/`，输出会写入 `output/CRC2/`；03a 输出会写入 `output/CRC2_03a/`。

## inferCNV 路线说明

Python 版有两条 03 步路线：

- `03a_run_infercnv_external_r.py`：调用 `resources/R/run_infercnv_external.R`，更接近原始 R inferCNV 流程。
- `03b_run_infercnvpy.py`：使用 `infercnvpy` 的纯 Python 路线，结果是近似替代。

两条路线都会生成统一文件：

```text
intermediate/<运行名>/03_cnv_calls.tsv
```

该文件包含：

```text
cell_ID    CNVLabel    cnv_score
```

04 之后的脚本只读取这个统一文件，不需要知道 03 步选择了哪条路线。

## 常见检查点

- `sample_name` 默认来自 `00_config.py`，也可由 `run_all.ps1 -SampleName` 临时指定。
- 输入目录名必须和样本名一致，例如 `CRC1` 对应 `input/CRC1/`。
- 新样本不要覆盖旧样本目录，直接新建 `input/<新样本名>/`。
- 输出目录可以删除后重跑，但 `input/<样本名>/` 不应被流程脚本写入或清空。
