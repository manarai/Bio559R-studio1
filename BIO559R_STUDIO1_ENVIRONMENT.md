# BIO559R Studio 1 Conda environment

This repository provides **`Bio559R-studio1-env`**, a CPU-oriented Conda environment for the Studio 1 single-cell RNA-seq workflow. **JupyterLab and a Jupyter kernel are included**; no separate Jupyter installation is required.

| Component | Role in Studio 1 |
|---|---|
| **Scanpy** | Python workflow for quality control, preprocessing, visualization, and exploratory single-cell analysis. |
| **AnnData** | Annotated matrix format used by Scanpy (`.h5ad`). |
| **Scrublet** | Optional computational doublet scoring; scores are evidence to inspect, not automatic deletion instructions. |
| **scvi-tools** | Deep generative single-cell models, including batch-aware modeling and integration methods used in later analyses. |
| **JupyterLab** | Browser-based notebook interface included in this environment. |

> **Scope and safety:** This environment is for post-Cell-Ranger count matrices. Keep raw FASTQs and original Cell Ranger outputs unchanged. Do not upload identifiable, protected, unpublished, embargoed, or access-controlled data to a public notebook, Colab, or external AI system.

## 1. Clone the repository

```bash
git clone https://github.com/manarai/Bio559R-studio1.git
cd Bio559R-studio1
```

If you already have the repository, update it first:

```bash
git pull --ff-only
```

## 2. Create and activate the environment

Use one command to create the environment. `mamba` or `micromamba` is generally faster than Conda for solving scientific environments.

```bash
# Preferred, if Mamba is installed
mamba env create -f Bio559R-studio1-env.yml

# Or use Conda
conda env create -f Bio559R-studio1-env.yml
```

Then activate the named environment:

```bash
conda activate Bio559R-studio1-env
```

If `conda activate` is unavailable after a new Conda/Mamba installation, initialize your shell once (for example, `conda init bash`), close and reopen the terminal, then activate the environment again.

## 3. Verify the installation

Run the following before opening a notebook:

```bash
python - <<'PY'
from importlib.metadata import version
import anndata as ad
import scanpy as sc
import scvi
import scrublet

print("Scanpy:", version("scanpy"))
print("AnnData:", version("anndata"))
print("Scrublet:", version("scrublet"))
print("scvi-tools:", version("scvi-tools"))
print("JupyterLab:", version("jupyterlab"))
print("AnnData import:", ad.AnnData)
print("Scanpy import:", sc.__name__)
print("scvi import:", scvi.__name__)
PY
```

## 4. Launch Studio 1 in JupyterLab

```bash
python -m ipykernel install --user \
  --name Bio559R-studio1-env \
  --display-name "Python (BIO559R Studio 1)"

jupyter lab
```

In JupyterLab, select **Kernel → Change Kernel → Python (BIO559R Studio 1)** before running a Studio 1 notebook.

## 5. Compatibility notes

- The environment pins **Python 3.12** because the selected scvi-tools, Scanpy, and AnnData Conda-forge builds support it together.
- The file is a **CPU baseline**. For Nvidia GPU acceleration, install a GPU-compatible PyTorch/CUDA stack according to the official [scvi-tools installation guidance](https://docs.scvi-tools.org/en/stable/installation.html) in a separate tested environment.
- Scrublet should receive raw, unnormalized count matrices. Review scores and experiment-specific expectations before deciding whether any cells should be excluded.
- scvi-tools is not a reason to skip QC or apply batch integration automatically. Inspect sample/batch structure and define a diagnostic plan before selecting an integration model.

## 6. Updating the environment

After changing `Bio559R-studio1-env.yml`, update the existing environment with:

```bash
conda env update --name Bio559R-studio1-env \
  --file Bio559R-studio1-env.yml \
  --prune
```

Re-run the verification commands afterward. For teaching reproducibility, record the resolved package set for a particular class run:

```bash
conda env export --name Bio559R-studio1-env --no-builds > Bio559R-studio1-env-resolved.yml
```
