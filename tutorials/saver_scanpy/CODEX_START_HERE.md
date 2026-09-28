# Codex-assisted student start guide

Use Codex as a **teaching assistant and code reviewer**, not as an automatic analyst. Students remain responsible for experimental design, QC choices, validation, and interpretation.

> **Non-negotiable safeguards:** do not upload identifiable or restricted human data to a public AI service or a personal Colab Drive. Keep raw FASTQs and Cell Ranger outputs unchanged. Do not let an AI agent delete files, change the Conda environment, or invent cell-type labels or scientific conclusions.

## 1. Create and verify the environment yourself

From the repository root, use `mamba` if available (it solves environments faster), otherwise use `conda`.

```bash
# Choose ONE of the next two commands.
mamba env create -f tutorials/saver_scanpy/environment.yml
# conda env create -f tutorials/saver_scanpy/environment.yml

conda activate saver-scanpy

# Confirm the teaching stack before opening the notebook.
python - <<'PY'
import scanpy as sc
import muon as mu
import anndata as ad
print('scanpy:', sc.__version__)
print('muon:', mu.__version__)
print('anndata:', ad.__version__)
PY
Rscript -e 'library(SAVER); print(packageVersion("SAVER"))'

# Register a clearly named kernel and start Jupyter.
python -m ipykernel install --user --name saver-scanpy --display-name 'Python (saver-scanpy)'
jupyter lab tutorials/saver_scanpy/preprocess_saver_scanpy_rna_atac.ipynb
```

If environment solving fails on a managed cluster, ask the HPC team to provide a Conda/Mamba module or a writable project location. Do **not** install packages into the shared `base` environment.

## 2. Stage only derived matrices, never raw FASTQs

The notebook expects one of these local layouts. Make a copy or symbolic link to Cell Ranger / Cell Ranger ARC **filtered** output; do not rename or edit the original output.

```text
data/
├── rna/filtered_feature_bc_matrix.h5             # RNA-only (Cell Ranger)
├── atac/filtered_peak_bc_matrix.h5               # ATAC-only (Cell Ranger ATAC), if applicable
├── atac/atac_fragments.tsv.gz                     # optional but required for fragment-based QC
└── multiome/filtered_feature_bc_matrix.h5         # RNA + ATAC (Cell Ranger ARC), if applicable
```

For multiome output, place `atac_fragments.tsv.gz` beside the HDF5 matrix or update the `MULTIOME_FRAGMENTS` variable in the notebook. The student should choose **one** input route (RNA-only, ATAC-only, or multiome) in the data-contract cell.

## 3. Give Codex this bounded prompt

Open the repository folder in the Codex interface of your choice, paste the prompt below, and retain control of execution. If the interface supports terminal permissions, allow **read-only** access until you understand every suggested command.

```text
You are a BIO559R teaching assistant for a student processing a 10x single-cell dataset.
Use tutorials/saver_scanpy/preprocess_saver_scanpy_rna_atac.ipynb and its README as the source of truth.

Before suggesting code, ask me for: species/reference assembly, assay (RNA-only, ATAC-only, or multiome), Cell Ranger or Cell Ranger ARC version, location of the filtered HDF5 matrix, whether a fragments file exists, sample/batch metadata, and the biological question. Explain each QC metric and propose starting thresholds as values to inspect—not universal rules.

Hard constraints:
- Never modify, delete, or overwrite raw FASTQs, Cell Ranger outputs, or the AnnData counts layer.
- Preserve unnormalized RNA counts in rna.layers['counts']; never pass log-normalized values or ATAC peaks to SAVER.
- Treat SAVER as optional RNA expression recovery after QC. Keep its estimates separate from raw counts and do not use them as a default input for differential expression.
- Do not claim a cluster identity or a biological result without marker evidence and metadata-aware validation.
- Flag all code that is computationally expensive, especially SAVER, fragment-based ATAC QC, or large dense conversions.
- Return a short explanation, a reversible code cell, and a validation check for every step.
- Do not run shell commands, install packages, upload data, or make repository changes; show commands for me to review and run.
```

## 4. Use AI only at the notebook checkpoints

| Checkpoint | Ask Codex to help with | Student must decide / verify |
| --- | --- | --- |
| A. Data contract | Identify input route and expected files | Assay, species, source paths, data-sharing permission |
| B. RNA QC | Explain distributions and propose transparent starting cutoffs | Whether cutoffs remove a plausible cell population; doublet strategy |
| C. RNA analysis | Review layers, neighbors, clustering, and marker plots | Batch handling, resolution, annotation evidence |
| D. Optional SAVER | Estimate runtime/memory and design a small pilot | Whether recovery is warranted; use only post-QC RNA UMI counts |
| E. ATAC QC | Explain fragment, TSS, nucleosome, and LSI diagnostics | Peak/cell filters and genome annotation compatibility |
| F. Reporting | Draft a methods log and figure captions | Final scientific claims and reproducibility record |

## 5. Colab option

Use the notebook's **Open in Colab** badge only for a small, non-sensitive teaching dataset or an approved de-identified subset. The Conda environment and R SAVER workflow are the reproducible route for local/HPC work. Colab can install the Python analysis stack inside the notebook, but has ephemeral storage and is not the default place for full datasets.

## 6. Completion evidence

Before considering the tutorial complete, save and retain:

1. the exact `environment.yml`, resolved package versions, and notebook;
2. a table of sample IDs, chemistry, reference, and input paths;
3. pre- and post-filter RNA/ATAC QC summaries with a written rationale;
4. raw count objects and processed `.h5ad` / `.h5mu` outputs under a new `results/` directory;
5. the exact SAVER command and RDS output only if SAVER was used; and
6. a short statement of what remains exploratory versus validated.
