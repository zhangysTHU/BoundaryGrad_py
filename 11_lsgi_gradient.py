from __future__ import annotations

# 11：把 Python 版 Cottrazm 中间结果导出给 R/LSGI，并生成细胞组分空间梯度图。
# 输入：
# - intermediate/05_TumorST_boundary_defined.h5ad
# - intermediate/07_DeconData.tsv
# 输出：
# - intermediate/lsgi_input/*.csv
# - intermediate/11_lsgi_cell_component_result.rds.gz
# - output/11_lsgi_gradient/*
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import anndata as ad
except ImportError:  # pragma: no cover - scanpy env normally provides anndata
    ad = None

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.config_loader import load_config


def _as_bool_arg(value: str | int | bool) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).lower() in {"1", "true", "t", "yes", "y"}


def _find_rscript(configured: str) -> Path:
    rscript = Path(configured)
    if rscript.exists():
        return rscript
    found = shutil.which("Rscript")
    if found:
        return Path(found)
    raise FileNotFoundError("Cannot find Rscript. Edit R_EXE in 00_config.py.")


def _get_lowres_scale_and_image(adata, spaceranger: Path, image_key: str) -> tuple[float, Path | None]:
    image_path = spaceranger / "spatial" / f"tissue_{image_key}_image.png"
    if not image_path.exists() and image_key != "lowres":
        image_path = spaceranger / "spatial" / "tissue_lowres_image.png"
    if not image_path.exists():
        image_path = None

    spatial = adata.uns.get("spatial", {})
    if spatial:
        library_id = next(iter(spatial))
        scalefactors = spatial[library_id].get("scalefactors", {})
        scale = scalefactors.get(f"tissue_{image_key}_scalef")
        if scale is None and image_key != "lowres":
            scale = scalefactors.get("tissue_lowres_scalef")
        if scale is not None:
            return float(scale), image_path

    scale_snapshot = adata.uns.get("spatial_scale_factors", {})
    scale = scale_snapshot.get(f"tissue_{image_key}_scalef")
    if scale is None and image_key != "lowres":
        scale = scale_snapshot.get("tissue_lowres_scalef")
    return float(scale or 1.0), image_path


def _export_lsgi_inputs(adata_path: Path, decon_path: Path, out_dir: Path, spaceranger: Path, image_key: str) -> Path | None:
    if ad is None:
        raise ImportError("anndata is not installed in this Python environment.")

    adata = ad.read_h5ad(adata_path)
    decon = pd.read_csv(decon_path, sep="\t")
    if "cell_ID" not in decon.columns:
        raise ValueError(f"{decon_path} must contain a cell_ID column.")
    if "Location" not in adata.obs.columns:
        raise ValueError(f"{adata_path} must contain obs['Location'] from step 05.")

    scale, image_path = _get_lowres_scale_and_image(adata, spaceranger, image_key)
    if "spatial" in adata.obsm:
        coords = np.asarray(adata.obsm["spatial"], dtype=float) * scale
        coord_df = pd.DataFrame({"cell_ID": adata.obs_names, "X": coords[:, 0], "Y": coords[:, 1]})
    elif {"imagecol", "imagerow"}.issubset(adata.obs.columns):
        coord_df = pd.DataFrame(
            {
                "cell_ID": adata.obs_names,
                "X": adata.obs["imagecol"].astype(float).to_numpy() * scale,
                "Y": adata.obs["imagerow"].astype(float).to_numpy() * scale,
            }
        )
    else:
        raise ValueError(f"{adata_path} must contain obsm['spatial'] or obs imagecol/imagerow coordinates.")

    boundary = pd.DataFrame({"cell_ID": adata.obs_names, "Location": adata.obs["Location"].astype(str).to_numpy()})

    common = pd.Index(coord_df["cell_ID"]).intersection(decon["cell_ID"]).intersection(boundary["cell_ID"])
    if len(common) < 10:
        raise ValueError("Too few matched spots across boundary object and deconvolution output.")

    coord_df = coord_df.set_index("cell_ID").loc[common].reset_index()
    boundary = boundary.set_index("cell_ID").loc[common].reset_index()
    decon = decon.set_index("cell_ID").loc[common].reset_index()

    out_dir.mkdir(parents=True, exist_ok=True)
    coord_df.to_csv(out_dir / "spatial_coords.csv", index=False)
    boundary.to_csv(out_dir / "boundary_labels.csv", index=False)
    decon.to_csv(out_dir / "cell_component_embeddings.csv", index=False)
    return image_path


def main() -> None:
    cfg = load_config()
    paths = cfg.PATHS
    params = cfg.PARAMS

    parser = argparse.ArgumentParser(description="Run LSGI cell-component gradients for the Python Cottrazm workflow.")
    parser.add_argument("--reuse-lsgi", action="store_true", help="Reuse intermediate/11_lsgi_cell_component_result.rds.gz and only redraw outputs.")
    parser.add_argument("--n-grids-scale", type=float, default=params.get("lsgi_n_grids_scale", 10))
    parser.add_argument("--n-cells-per-meta", type=float, default=params.get("lsgi_n_cells_per_meta", 50))
    parser.add_argument("--r-squared-thresh", type=float, default=params.get("lsgi_r_squared_thresh", 0.3))
    parser.add_argument("--minimum-fctr", type=float, default=params.get("lsgi_minimum_fctr", 3))
    parser.add_argument("--arrow-length-scale", type=float, default=params.get("lsgi_arrow_length_scale", 1.4))
    parser.add_argument("--arrow-linewidth", type=float, default=params.get("lsgi_arrow_linewidth", 1.0))
    parser.add_argument("--arrow-head-cm", type=float, default=params.get("lsgi_arrow_head_cm", 0.20))
    parser.add_argument("--arrow-closed", type=_as_bool_arg, default=params.get("lsgi_arrow_closed", True))
    parser.add_argument("--image-key", default=params.get("lsgi_image_key", "lowres"), choices=["lowres", "hires"])
    args = parser.parse_args()

    input_dir = paths["intermediate"] / "lsgi_input"
    output_dir = paths["output"] / "11_lsgi_gradient"
    image_path = _export_lsgi_inputs(
        adata_path=paths["intermediate"] / "05_TumorST_boundary_defined.h5ad",
        decon_path=paths["intermediate"] / "07_DeconData.tsv",
        out_dir=input_dir,
        spaceranger=paths["spaceranger"],
        image_key=args.image_key,
    )

    rscript = _find_rscript(cfg.R_EXE)
    lsgi_root = Path(__file__).resolve().parent.parent / "LSGI-master"
    if not lsgi_root.exists():
        raise FileNotFoundError(f"Cannot find LSGI-master at {lsgi_root}")

    cmd = [
        str(rscript),
        str(paths["resources"] / "R" / "run_lsgi_gradient.R"),
        f"--input-dir={input_dir}",
        f"--output-dir={output_dir}",
        f"--intermediate-dir={paths['intermediate']}",
        f"--lsgi-root={lsgi_root}",
        f"--sample-name={cfg.SAMPLE_NAME}",
        f"--n-grids-scale={args.n_grids_scale}",
        f"--n-cells-per-meta={args.n_cells_per_meta}",
        f"--r-squared-thresh={args.r_squared_thresh}",
        f"--minimum-fctr={args.minimum_fctr}",
        f"--arrow-length-scale={args.arrow_length_scale}",
        f"--arrow-linewidth={args.arrow_linewidth}",
        f"--arrow-head-cm={args.arrow_head_cm}",
        f"--arrow-closed={int(args.arrow_closed)}",
        f"--reuse-lsgi={int(args.reuse_lsgi)}",
    ]
    if image_path is not None:
        cmd.append(f"--image-path={image_path}")

    subprocess.run(cmd, cwd=Path(__file__).resolve().parent, check=True)


if __name__ == "__main__":
    main()
