"""Generate the completed and no-code student notebooks for the BIO559R lesson."""
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parent
KERNEL = {
    "display_name": "Python (bio559r-scrnaseq)",
    "language": "python",
    "name": "bio559r-scrnaseq",
}
LANGUAGE = {"name": "python", "version": "3.11"}


def markdown(source: str):
    return nbf.v4.new_markdown_cell(source.strip() + "\n")


def code(source: str):
    return nbf.v4.new_code_cell(source.strip() + "\n")


def blank_code_cell():
    """An intentionally empty code cell for the no-code student notebook."""
    return nbf.v4.new_code_cell("")


AI_PROMPT_STUDIO = """
## Optional Colab AI prompt studio

> **Use AI as a coach, not an answer key.** Keep the canonical workflow, rerun every suggested change, and judge the result from diagnostics and sensitivity evidence.

**Privacy boundary:** never send protected, identifiable, unpublished, embargoed, or access-controlled data to an AI assistant.

### Prompt 1 — Explain the decision
> I am preprocessing two 10x scRNA-seq count matrices before any embedding or clustering. Explain how per-sample total UMI counts, detected genes, and mitochondrial fraction inform a QC decision. Name assumptions and failure modes. Do not choose thresholds or write my assessed conclusion.

### Prompt 2 — Diagnose before editing
> Here is a minimal non-identifying QC summary or traceback: **[paste summary]**. Give three ranked hypotheses, one discriminating check for each, and the smallest reversible change. Preserve raw inputs, the raw-count layer, and provenance.

### Prompt 3 — Design one bounded variation
> I will change **one** QC threshold for **one** sample. Predict the measurable diagnostic effect, show a compact patch only after explaining it, and state what result would make me retain, narrow, or withdraw the interpretation.

### Prompt 4 — Audit a claim
> Evidence artifact: my pre- and post-filter QC table plus the saved `.h5ad` object. Separate direct observation, inference, limitation, and next validation test. Flag unsupported causal language and unverified citations. Do not invent evidence.

**Verification record:** record the prompt purpose, what you accepted or rejected, the diagnostic you reran, and whether the claim was retained, narrowed, or withdrawn.
"""


complete_cells = [
    markdown(
        """
# BIO559R: Preprocess and combine two scRNA-seq count matrices

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/manarai/Bio559R-studio1/blob/main/tutorials/scrnaseq_preprocessing/scrnaseq_two_samples_preprocessing_complete.ipynb)

**A short, commented Scanpy/scverse workflow for two 10x Genomics `filtered_feature_bc_matrix.h5` files.**

## What this notebook does

1. Loads two Cell Ranger gene-expression HDF5 matrices.
2. checks their genes and sample metadata;
3. combines the matrices while retaining the **sample** label;
4. preserves raw UMI counts in `adata.layers["counts"]`;
5. calculates and inspects quality-control (QC) metrics by sample;
6. applies **documented, sample-specific** QC filters;
7. optionally evaluates doublets without automatically deleting cells;
8. normalizes, log-transforms, and selects batch-aware highly variable genes (HVGs);
9. saves a reproducible preprocessing object and parameter record.

## Deliberate boundary

This lesson **stops after preprocessing and feature selection**. It does **not** run PCA, neighbors, UMAP, clustering, cell-type annotation, differential expression, or batch correction. Concatenating two samples is data combination; it is **not** batch-effect correction or biological integration.

> **Data and AI safety:** use only small, approved, non-sensitive teaching data in Colab. Do not upload identifiable, protected, unpublished, embargoed, or access-controlled data to Colab or an external AI assistant. Keep the original Cell Ranger outputs unchanged.

## Data contract

The two inputs must be 10x Gene Expression HDF5 matrices (usually named `filtered_feature_bc_matrix.h5`) created with compatible genome/reference annotations. This notebook uses `scanpy.read_10x_h5`; it does **not** load `.h5ad` files. If your files are `.h5ad`, stop and choose the correct loader rather than renaming files.
"""
    ),
    markdown(AI_PROMPT_STUDIO),
    markdown(
        """
## 0. Install packages when working in Google Colab

For a local computer, create and activate the Conda environment in this folder **before** opening the notebook. For Colab, run the next cell once per fresh runtime. Package installation may require a runtime restart before the import cell.
"""
    ),
    code(
        """
# Purpose: install the Python packages only when this notebook runs in Google Colab.
# Why: Colab runtimes are temporary and do not use the repository's Conda environment.
# Local/HPC users should activate the bio559r-scrnaseq Conda environment instead.
import sys

IN_COLAB = "google.colab" in sys.modules
if IN_COLAB:
    %pip install --quiet "scanpy==1.11.5" "anndata==0.12.8" "scrublet==0.2.3" python-igraph leidenalg scikit-misc
    print("Installed the Python preprocessing stack for this Colab runtime.")
    print("If the import cell fails, restart the Colab runtime once and run from the top.")
else:
    print("Not in Colab. Activate the 'bio559r-scrnaseq' Conda environment before continuing.")
"""
    ),
    markdown(
        """
## 1. Import tools and record versions

We use **AnnData** as the cell-by-gene data object and **Scanpy** for preprocessing. Recording versions makes a result easier to reproduce later.
"""
    ),
    code(
        """
# Purpose: import the tools used throughout the notebook and record their versions.
# Why: version information is part of the analysis provenance.
from pathlib import Path
import json
import platform
import sys

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 2
sc.set_figure_params(dpi=100, facecolor="white")

print("Python:", sys.version.split()[0])
print("Scanpy:", sc.__version__)
print("AnnData:", ad.__version__)
"""
    ),
    markdown(
        """
## 2. Define the two input samples and their minimal metadata

Give each matrix a clear sample identifier. The identifiers are carried into `adata.obs["sample"]` after the matrices are combined. Replace the placeholder **condition** and **batch** values with your own metadata; these labels are not used to filter cells in this preprocessing lesson.

**Teaching point:** sample labels are essential. Cells are not independent biological replicates, and a two-sample cell-level analysis does not by itself establish a condition-level result.
"""
    ),
    code(
        """
# Purpose: declare input locations, sample IDs, and simple sample-level metadata.
# Why: labels are needed for per-sample QC and for later, metadata-aware analyses.
# Edit this cell before loading data. Do not move or overwrite the original HDF5 files.
SPECIES = "human"  # Choose exactly "human" or "mouse" for the mitochondrial-gene prefix.
INPUT_MODE = "upload" if IN_COLAB else "paths"  # Use "upload" in Colab; use "paths" elsewhere.
SAMPLE_IDS = ["sample_A", "sample_B"]

# Local/HPC path mode only: replace both paths with approved file locations.
LOCAL_H5_PATHS = [
    Path("/path/to/sample_A/filtered_feature_bc_matrix.h5"),
    Path("/path/to/sample_B/filtered_feature_bc_matrix.h5"),
]

# Add only metadata that you can verify. These labels will be copied to every cell.
SAMPLE_INFO = pd.DataFrame(
    {
        "sample": SAMPLE_IDS,
        "condition": ["condition_A", "condition_B"],
        "batch": ["batch_1", "batch_2"],
    }
)

OUTPUT_DIR = Path("results")
OUTPUT_DIR.mkdir(exist_ok=True)

assert SPECIES in {"human", "mouse"}
assert len(SAMPLE_IDS) == 2 and len(set(SAMPLE_IDS)) == 2
assert SAMPLE_INFO["sample"].tolist() == SAMPLE_IDS
SAMPLE_INFO
"""
    ),
    markdown(
        """
## 3. Obtain exactly two HDF5 files

In Colab, the next cell opens a file picker. Upload only the two approved 10x HDF5 matrices. Outside Colab, it checks the two paths from the previous cell.
"""
    ),
    code(
        """
# Purpose: obtain two HDF5 input paths without changing the source files.
# Why: a clear two-file contract prevents accidental use of the wrong matrix or a duplicate file.
if INPUT_MODE == "upload":
    if not IN_COLAB:
        raise RuntimeError("INPUT_MODE='upload' works only in Google Colab. Use INPUT_MODE='paths' locally.")
    from google.colab import files

    uploaded = files.upload()  # Upload exactly two approved 10x Gene Expression .h5 files.
    uploaded_h5 = sorted(Path(name) for name in uploaded if Path(name).suffix.lower() == ".h5")
    if len(uploaded_h5) != 2:
        raise ValueError(f"Expected exactly two .h5 files; received: {[path.name for path in uploaded_h5]}")
    H5_PATHS = uploaded_h5
elif INPUT_MODE == "paths":
    H5_PATHS = LOCAL_H5_PATHS
else:
    raise ValueError("INPUT_MODE must be 'upload' or 'paths'.")

for sample_id, path in zip(SAMPLE_IDS, H5_PATHS):
    if not path.exists():
        raise FileNotFoundError(f"{sample_id}: file not found at {path}. Edit the previous cell; do not guess a path.")
    print(f"{sample_id}: {path}")
"""
    ),
    markdown(
        """
## 4. Load each 10x gene-expression matrix separately

Loading each file first makes it possible to audit dimensions and gene names before combining samples. `var_names_make_unique()` prevents duplicated feature names from breaking later steps.
"""
    ),
    code(
        """
# Purpose: read each 10x Gene Expression HDF5 file into its own AnnData object.
# Why: inspecting inputs separately catches incompatible files before samples are combined.
def read_10x_sample(sample_id: str, path: Path) -> ad.AnnData:
    sample_adata = sc.read_10x_h5(path, gex_only=True)
    sample_adata.var_names_make_unique()
    sample_adata.obs["source_file"] = path.name
    print(f"{sample_id}: {sample_adata.n_obs:,} cells x {sample_adata.n_vars:,} genes")
    return sample_adata

adatas = {
    sample_id: read_10x_sample(sample_id, path)
    for sample_id, path in zip(SAMPLE_IDS, H5_PATHS)
}
"""
    ),
    markdown(
        """
## 5. Audit shared features, then combine the samples

The two inputs should normally come from the same species and compatible Cell Ranger reference. The notebook keeps only the **intersection** of genes (`join="inner"`) so both samples have the same feature space. A low overlap is a reason to stop and verify the sources—not a reason to force a merge.
"""
    ),
    code(
        """
# Purpose: compare gene sets and concatenate the two samples with an explicit sample label.
# Why: combination is only defensible when the matrices share a compatible feature space.
gene_sets = {sample_id: set(sample_adata.var_names) for sample_id, sample_adata in adatas.items()}
shared_genes = set.intersection(*gene_sets.values())

for sample_id, genes in gene_sets.items():
    overlap_percent = 100 * len(shared_genes) / len(genes)
    print(f"{sample_id}: {len(genes):,} genes; shared genes: {len(shared_genes):,} ({overlap_percent:.1f}%)")

if len(shared_genes) == 0:
    raise ValueError("No shared genes. Stop and verify file type, species, and reference annotation.")

adata = ad.concat(adatas, label="sample", join="inner", index_unique="-")
adata.obs = adata.obs.join(SAMPLE_INFO.set_index("sample"), on="sample")
adata.layers["counts"] = adata.X.copy()  # Invariant: this layer contains unnormalized UMI counts.

print(adata)
print("Raw counts preserved:", "counts" in adata.layers)
display(adata.obs[["sample", "condition", "batch", "source_file"]].head())
"""
    ),
    markdown(
        """
## 6. Calculate QC metrics before filtering

We examine three standard diagnostics **by sample**:

- `total_counts`: total UMIs/library size;
- `n_genes_by_counts`: number of detected genes;
- `pct_counts_mt`: percent of counts from mitochondrial genes.

High mitochondrial fraction can reflect damaged cells, but it can also depend on tissue and biology. High count/gene values can identify potential doublets, but they are not proof. Plot first; choose thresholds second.
"""
    ),
    code(
        """
# Purpose: mark mitochondrial genes and calculate cell-level QC metrics on raw counts.
# Why: QC values must be inspected before any filtering or normalization changes the working matrix.
MITO_PREFIX = "MT-" if SPECIES == "human" else "mt-"
adata.var["mt"] = adata.var_names.str.startswith(MITO_PREFIX)

if adata.var["mt"].sum() == 0:
    raise ValueError(
        f"No mitochondrial genes matched the prefix {MITO_PREFIX!r}. "
        "Check species, gene naming, and reference annotation before continuing."
    )

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True,
)

qc_columns = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
qc_summary = adata.obs.groupby("sample", observed=True)[qc_columns].describe().round(2)
display(qc_summary)
"""
    ),
    code(
        """
# Purpose: visualize QC distributions overall and by sample.
# Why: distributions—not universal thresholds—should guide a transparent filtering decision.
sc.pl.violin(
    adata,
    keys=["total_counts", "n_genes_by_counts", "pct_counts_mt"],
    groupby="sample",
    jitter=0.25,
    rotation=30,
)

sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="sample")
"""
    ),
    markdown(
        """
## 7. Decide and document sample-specific QC rules

The values below are deliberately conservative **starting points**, not universal recommendations. The completed code runs with only a minimum detected-gene rule unless you choose additional limits after reviewing the plots. Set `None` when you are not applying a filter.

> Record the biological/technical reason for every active cutoff in a methods log. Never choose a cutoff solely because it is common in an online tutorial.
"""
    ),
    code(
        """
# Purpose: define transparent per-sample QC rules after inspecting the preceding diagnostics.
# Why: different samples can have different QC distributions, and every active threshold should be justified.
QC_RULES = {
    "sample_A": {"min_genes": 200, "max_genes": None, "max_pct_mt": None},
    "sample_B": {"min_genes": 200, "max_genes": None, "max_pct_mt": None},
}
MIN_CELLS_PER_GENE = 3

if set(QC_RULES) != set(SAMPLE_IDS):
    raise ValueError("QC_RULES must contain exactly one entry for each sample ID.")

pd.DataFrame(QC_RULES).T
"""
    ),
    code(
        """
# Purpose: apply the documented cell filters and a low-information gene filter.
# Why: filtering occurs only after QC review and preserves the raw-count layer for retained cells.
keep = np.ones(adata.n_obs, dtype=bool)
filter_report = []

for sample_id, rule in QC_RULES.items():
    sample_mask = adata.obs["sample"].to_numpy() == sample_id
    sample_obs = adata.obs.loc[sample_mask]
    sample_keep = np.ones(sample_mask.sum(), dtype=bool)

    sample_keep &= sample_obs["n_genes_by_counts"].to_numpy() >= rule["min_genes"]
    if rule["max_genes"] is not None:
        sample_keep &= sample_obs["n_genes_by_counts"].to_numpy() <= rule["max_genes"]
    if rule["max_pct_mt"] is not None:
        sample_keep &= sample_obs["pct_counts_mt"].to_numpy() <= rule["max_pct_mt"]

    keep[sample_mask] = sample_keep
    filter_report.append(
        {"sample": sample_id, "before_cells": int(sample_mask.sum()), "after_cell_filters": int(sample_keep.sum())}
    )

adata = adata[keep].copy()
sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)

filter_report = pd.DataFrame(filter_report)
filter_report["after_gene_filter_total_cells"] = adata.n_obs
filter_report["genes_after_filter"] = adata.n_vars
display(filter_report)
print("Raw counts preserved after subsetting:", "counts" in adata.layers)
print(adata)
"""
    ),
    markdown(
        """
## 8. Optional: evaluate doublets, then make a separate decision

Scanpy's Scrublet wrapper should receive **unnormalized counts** and, for two samples, uses the sample label as `batch_key`. It adds a score and a predicted call; it does not prove that a cell is a doublet. The cell below is disabled by default. If you enable it, inspect the score distributions and record why you retain or remove any cells.
"""
    ),
    code(
        """
# Purpose: optionally calculate Scrublet doublet scores on unnormalized counts.
# Why: doublet calls are evidence to inspect, not an automatic deletion rule.
RUN_SCRUBLET = False
EXPECTED_DOUBLET_RATE = 0.05  # Replace only with a rate justified by the experiment/loading design.

if RUN_SCRUBLET:
    sc.pp.scrublet(
        adata,
        batch_key="sample",
        expected_doublet_rate=EXPECTED_DOUBLET_RATE,
        random_state=0,
    )
    display(adata.obs.groupby("sample", observed=True)[["doublet_score", "predicted_doublet"]].agg(["mean", "sum"]))
    sc.pl.violin(adata, keys="doublet_score", groupby="sample", jitter=0.25)
else:
    print("Scrublet skipped. Set RUN_SCRUBLET=True only after reviewing the design and expected rate.")
"""
    ),
    code(
        """
# Purpose: make an explicit, reversible decision about predicted doublets.
# Why: a predicted label should not silently change the analysis population.
REMOVE_PREDICTED_DOUBLETS = False

if REMOVE_PREDICTED_DOUBLETS:
    if "predicted_doublet" not in adata.obs:
        raise RuntimeError("Run the optional Scrublet cell before removing predicted doublets.")
    before_doublet_filter = adata.n_obs
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    print(f"Cells retained after doublet decision: {adata.n_obs:,} / {before_doublet_filter:,}")
else:
    print("No cells removed because of a doublet prediction.")
"""
    ),
    markdown(
        """
## 9. Normalize and log-transform the working matrix

`adata.X` changes in this step, while `adata.layers["counts"]` remains raw UMI counts. Normalization makes cell library sizes more comparable for exploratory expression analyses; it does not make two biological samples equivalent.
"""
    ),
    code(
        """
# Purpose: normalize the working matrix and apply a log(1 + x) transform.
# Why: downstream feature selection expects a comparable, transformed expression scale.
TARGET_SUM = 10_000

sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
sc.pp.log1p(adata)

print("Working matrix: normalized and log-transformed.")
print("Invariant check — raw counts still available:", "counts" in adata.layers)
"""
    ),
    markdown(
        """
## 10. Select batch-aware highly variable genes (HVGs)

HVGs are informative features for a later dimension-reduction/integration lesson. The `batch_key="sample"` argument helps avoid selecting genes that appear variable only because of one sample. We **mark** HVGs but do not run PCA or subset the object in this lesson.
"""
    ),
    code(
        """
# Purpose: mark highly variable genes while accounting for the two sample labels.
# Why: these features are a reproducible handoff to a later PCA/integration lesson.
N_TOP_HVGS = min(2_000, adata.n_vars)

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=N_TOP_HVGS,
    batch_key="sample",
    flavor="seurat",
)

sc.pl.highly_variable_genes(adata)
print(f"HVGs selected: {int(adata.var['highly_variable'].sum()):,}")
if "highly_variable_nbatches" in adata.var:
    display(adata.var.loc[adata.var["highly_variable"], ["highly_variable_nbatches"]].head())
"""
    ),
    markdown(
        """
## 11. Save the preprocessing handoff and provenance

This writes a **new** `.h5ad` file and JSON parameter record in `results/`. It does not modify the original two HDF5 files. The saved AnnData object retains raw counts, QC metrics, sample labels, normalized/log-transformed `X`, and the HVG annotation.
"""
    ),
    code(
        """
# Purpose: save the derived AnnData object and a human-readable provenance record.
# Why: a later analysis can start from the same documented preprocessing state.
provenance = {
    "input_files": {sample_id: str(path) for sample_id, path in zip(SAMPLE_IDS, H5_PATHS)},
    "species": SPECIES,
    "scanpy_version": sc.__version__,
    "anndata_version": ad.__version__,
    "python_version": sys.version.split()[0],
    "platform": platform.platform(),
    "qc_rules": QC_RULES,
    "min_cells_per_gene": MIN_CELLS_PER_GENE,
    "scrublet_run": RUN_SCRUBLET,
    "predicted_doublets_removed": REMOVE_PREDICTED_DOUBLETS,
    "target_sum": TARGET_SUM,
    "n_top_hvgs_requested": N_TOP_HVGS,
    "scope": "Preprocessing only; no PCA, neighbor graph, UMAP, clustering, annotation, DE, or batch correction.",
}

adata.uns["bio559r_preprocessing"] = {
    "counts_layer": "counts",
    "sample_ids": SAMPLE_IDS,
    "scope": provenance["scope"],
}

output_h5ad = OUTPUT_DIR / "two_sample_preprocessed.h5ad"
output_json = OUTPUT_DIR / "two_sample_preprocessing_parameters.json"
adata.write(output_h5ad, compression="gzip")
output_json.write_text(json.dumps(provenance, indent=2))

print("Saved:", output_h5ad)
print("Saved:", output_json)
"""
    ),
    markdown(
        """
## Completion checklist and handoff

- [ ] I confirmed both inputs are compatible 10x Gene Expression HDF5 matrices.
- [ ] I recorded species, reference build, Cell Ranger version, chemistry, sample IDs, and approved data location.
- [ ] I preserved raw UMI counts in `adata.layers["counts"]`.
- [ ] I inspected QC distributions separately for both samples before filtering.
- [ ] I documented every active QC threshold and any doublet decision.
- [ ] I saved the `.h5ad` object and JSON provenance record.
- [ ] I did **not** interpret an embedding, cluster, or cell type—none was created in this lesson.

### Next lesson (intentionally not performed here)

Use the saved object to inspect principal components and decide whether a batch-correction/integration method is warranted. Any correction should be selected and assessed with sample metadata, diagnostics, and a clear biological question—not applied simply because there are two files.
"""
    ),
]

student_steps = [
    (
        "0",
        "Set the learning and privacy boundary",
        "State that the task is to preprocess and combine two approved 10x Gene Expression HDF5 matrices, ending after highly variable genes. List excluded methods: PCA, UMAP, clustering, annotation, differential expression, and batch correction. Record whether the data may be uploaded to Colab or summarized for an AI assistant.",
        "I am building a two-sample scRNA-seq preprocessing notebook for a class. Explain the difference between concatenating two samples and batch correction. List the raw-data, privacy, and claim-boundary safeguards I should write before any code. Do not generate analysis code yet.",
    ),
    (
        "1",
        "Create the environment and import cell",
        "Create one small cell that detects Colab, installs the approved Scanpy/AnnData stack only in Colab, imports AnnData, Scanpy, pandas, NumPy, matplotlib, and path utilities, then prints Python and package versions.",
        "Generate one short, well-commented Jupyter code cell for this step. It must detect Google Colab, install Scanpy, AnnData, Scrublet, python-igraph, leidenalg, and scikit-misc only in Colab, then print versions after imports. Explain each section before the code. Do not combine later analysis steps into this cell.",
    ),
    (
        "2",
        "Declare the data contract and sample metadata",
        "Create a configuration cell with two meaningful sample IDs, species (human or mouse), two local paths or a Colab-upload mode, an output directory, and a two-row metadata table containing sample, condition, and batch. Add assertions for exactly two unique IDs.",
        "Generate one small, well-commented configuration cell. The input is exactly two approved 10x filtered Gene Expression HDF5 files, not H5AD files. Include a species variable and a two-row sample metadata table. Explain why sample-level metadata is needed and why cells are not biological replicates.",
    ),
    (
        "3",
        "Obtain and verify exactly two input files",
        "Create one cell that opens a Colab uploader when selected or checks local paths otherwise. It must reject anything other than two `.h5` files and print the sample-ID-to-path mapping without editing source files.",
        "Generate one reversible code cell that obtains exactly two HDF5 input paths in either Colab upload mode or local-path mode. It must validate file existence and never rename, overwrite, or move the original files. Before code, explain the data-contract check.",
    ),
    (
        "4",
        "Load each matrix separately",
        "Create one cell with a small helper function that uses the Scanpy 10x HDF5 reader for gene expression, makes feature names unique, records the source filename, and prints cells by genes for each sample.",
        "Generate one small cell that reads each sample separately with Scanpy's 10x HDF5 reader. Explain why it is safer to inspect each AnnData object before combining samples. Include a check that would help identify an incompatible input, but do not concatenate yet.",
    ),
    (
        "5",
        "Audit genes, combine samples, and preserve raw counts",
        "Create one cell that reports gene-set overlap, stops if there are no shared genes, concatenates using the gene intersection, creates an explicit sample label, joins sample metadata, and copies unnormalized UMI counts to a named counts layer.",
        "Generate one short, well-commented cell that compares the two gene sets, combines compatible matrices using only shared genes, carries a sample label into cell metadata, joins the two-row metadata table, and preserves raw UMI counts in a counts layer. Explain why combining is not batch correction.",
    ),
    (
        "6",
        "Calculate QC metrics on raw counts",
        "Create one cell that marks mitochondrial genes using the selected species prefix, calculates total counts, detected genes, and mitochondrial percentage, errors clearly if no mitochondrial genes are found, and shows a per-sample summary.",
        "Generate one code cell to calculate Scanpy QC metrics before normalization. It must use the species setting to choose the mitochondrial prefix and report QC distributions by sample. Explain what total counts, detected genes, and mitochondrial percentage can and cannot indicate.",
    ),
    (
        "7",
        "Visualize QC before filtering",
        "Create one plotting cell that shows total counts, detected genes, and mitochondrial percentage by sample, plus a counts-versus-genes scatter plot colored once by mitochondrial percentage and once by sample.",
        "Generate one compact plotting cell for pre-filter QC. Include per-sample violin plots and two diagnostic scatter plots. Explain what patterns could indicate low-quality cells or possible doublets, and why no universal threshold follows from a plot alone.",
    ),
    (
        "8",
        "Choose documented per-sample filters",
        "Write a separate configuration cell containing a transparent per-sample threshold dictionary: minimum detected genes, optional maximum detected genes, optional maximum mitochondrial percentage, plus a minimum cells-per-gene rule. Use `None` for inactive filters.",
        "Do not choose my thresholds. Instead, give me a template for a per-sample QC-rule dictionary and explain how to record the biological or technical rationale for every active filter. Include a bounded variation: change just one threshold for one sample and predict which diagnostic count should change.",
    ),
    (
        "9",
        "Apply filters and inspect retention",
        "Create one cell that applies only the documented rules, subsets safely, filters low-information genes, displays cells retained per sample, and verifies the raw counts layer still exists.",
        "Generate one reversible filtering cell that applies the QC-rule dictionary without hard-coding extra thresholds, filters genes expressed in too few cells, reports before/after counts by sample, and confirms that the raw counts layer is preserved. Explain how to undo or revise this decision by rerunning earlier cells.",
    ),
    (
        "10",
        "Optionally assess doublets",
        "Create a disabled-by-default cell that runs Scanpy Scrublet on unnormalized counts with the sample label as the batch key, reports scores/calls by sample, and does not remove cells automatically. Then create a separate explicit decision cell for optional removal.",
        "Generate two small cells: one disabled-by-default Scrublet scoring cell for raw counts with a sample batch key, and one separate, explicit decision cell for optional removal. Explain why a predicted doublet is not proof and why the decision must be recorded. Do not enable removal by default.",
    ),
    (
        "11",
        "Normalize and log-transform the working matrix",
        "Create one cell that normalizes total counts to a stated target sum and applies log1p. It must explain that the working matrix changes while the preserved counts layer does not.",
        "Generate one short cell to normalize total counts and apply log1p to the working matrix. Explain the purpose, state the target sum, and verify that the raw-count layer remains unchanged. Do not add scaling, PCA, or an embedding.",
    ),
    (
        "12",
        "Select batch-aware highly variable genes",
        "Create one cell that marks a bounded number of HVGs using the sample label as a batch key, plots the result, and does not run PCA or subset the object.",
        "Generate one short cell to mark highly variable genes after log normalization. It must use the sample label as a batch key, show the HVG plot, and stop before PCA. Explain why HVG selection is a handoff to a later analysis rather than biological integration.",
    ),
    (
        "13",
        "Save and audit the handoff",
        "Create one cell that saves a new compressed H5AD file plus a JSON provenance record containing inputs, package versions, QC rules, doublet choice, normalization target, and scope boundary. End with a written claim audit.",
        "Generate one well-commented save/provenance cell for a preprocessed AnnData object and a JSON parameter file. Then help me write a four-part audit: direct observation, inference, limitation, and next validation test. Do not claim that preprocessing corrected batch effects or identified cell types.",
    ),
]

student_cells = [
    markdown(
        """
# BIO559R: Build a two-sample scRNA-seq preprocessing notebook — student prompt studio

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/manarai/Bio559R-studio1/blob/main/tutorials/scrnaseq_preprocessing/scrnaseq_two_samples_preprocessing_student_guided.ipynb)

**This notebook intentionally contains no solution code.** Build one small cell at a time with your AI assistant, read every line, run it, and verify the diagnostic before moving on.

## Goal

Starting from two compatible 10x `filtered_feature_bc_matrix.h5` files, build a preprocessing-only Scanpy workflow that preserves raw counts, checks QC by sample, applies transparent filters, optionally assesses doublets, normalizes/log-transforms data, selects batch-aware HVGs, and saves a reusable `.h5ad` handoff.

## Boundary

Do **not** add PCA, neighbors, UMAP, clustering, annotation, differential expression, or batch correction. Combining samples is not batch correction.

> **Privacy rule:** Do not upload identifiable, protected, unpublished, embargoed, or access-controlled data to Colab or an external AI assistant. Give an AI only approved, non-identifying summaries and code/error text.

## How to use this notebook

1. Read one step below.
2. Paste its prompt into an AI assistant.
3. Ask for a **single small cell**, not an entire notebook.
4. Read and explain the proposed code before running it.
5. Paste the code into the empty cell directly beneath the step.
6. Record what you accepted or rejected and the diagnostic you reran.

### Required verification record

For each decision, write: **prompt purpose → accepted/rejected suggestion → diagnostic rerun → claim disposition (retained/narrowed/withdrawn)**.
"""
    ),
    markdown(AI_PROMPT_STUDIO),
]

for number, title, task, prompt in student_steps:
    student_cells.append(
        markdown(
            f"""
## Step {number}: {title}

### Your task
{task}

### AI prompt to paste

> {prompt}

### Before moving on

- Can you explain each line in the generated cell?
- What file, object, or metadata does it change?
- What diagnostic or assertion proves the step worked?
- Does the step preserve original inputs and raw UMI counts?
"""
        )
    )
    student_cells.append(blank_code_cell())

student_cells.append(
    markdown(
        """
## Final self-audit

Write short answers before submitting or continuing to the next lesson:

1. **Direct observation:** What changed in the number of cells and genes for each sample after your chosen QC rules?
2. **Inference:** What does the QC evidence support about the retained analysis population?
3. **Limitation:** Which technical or biological uncertainty remains? Why does preprocessing not correct batch effects?
4. **Next test:** What PCA/batch diagnostic would you inspect next, before selecting any integration approach?

### Completion checklist

- [ ] I used two compatible 10x Gene Expression HDF5 matrices.
- [ ] I recorded species, reference, chemistry, source paths, and sample metadata.
- [ ] I stored unnormalized counts in a named layer before preprocessing.
- [ ] I inspected QC by sample before filtering.
- [ ] I documented every active QC and doublet decision.
- [ ] I saved a derived `.h5ad` object and parameter record.
- [ ] I did not create or interpret clusters, embeddings, or cell types.
"""
    )
)


def write_notebook(path: Path, cells):
    notebook = nbf.v4.new_notebook(cells=cells, metadata={"kernelspec": KERNEL, "language_info": LANGUAGE})
    nbf.write(notebook, path)


write_notebook(ROOT / "scrnaseq_two_samples_preprocessing_complete.ipynb", complete_cells)
write_notebook(ROOT / "scrnaseq_two_samples_preprocessing_student_guided.ipynb", student_cells)
print("Wrote both notebooks to", ROOT)
