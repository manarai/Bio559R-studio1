# SAVER + Scanpy/Muon single-cell preprocessing module

This **stand-alone teaching module** adds a reproducible path from Cell Ranger count matrices to a quality-controlled exploratory analysis in **Scanpy** (RNA) and **Muon** (ATAC / multiome). It is designed for BIO559R students who have completed primary processing and need to begin downstream single-cell work.

> **Scope clarification:** in this module, **SAVER** means *Single-cell Analysis Via Expression Recovery*, the R package for optional expression recovery on **post-QC UMI RNA counts**. SAVER is **not** an ATAC method and is not a substitute for quality control, normalization, batch-aware modeling, or differential-expression design. The tutorial never overwrites raw counts; SAVER-derived values are kept separate and used for sensitivity checks or visualization only.

## Start here

| Item | Purpose |
| --- | --- |
| [Tutorial README](tutorials/saver_scanpy/README.md) | Data contracts, prerequisites, workflow map, installation, scientific guardrails, and progress checklist |
| [Conda environment](tutorials/saver_scanpy/environment.yml) | One environment containing Scanpy, Muon, Jupyter, and `r-saver` |
| [Studio 1 Conda environment](Bio559R-studio1-env.yml) | CPU baseline with Scanpy, AnnData, Scrublet, scvi-tools, and JupyterLab; see the [installation guide](BIO559R_STUDIO1_ENVIRONMENT.md). |
| [Interactive notebook](tutorials/saver_scanpy/preprocess_saver_scanpy_rna_atac.ipynb) | Student tutorial for RNA, ATAC, or 10x multiome output; includes AI checkpoints |
| [Codex start guide](tutorials/saver_scanpy/CODEX_START_HERE.md) | Safe prompts and commands for an AI-assisted, student-owned analysis |

## What this addition does—and does not do

- **Does:** provide an opt-in analysis module after Cell Ranger / Cell Ranger ARC output; preserve count matrices; support 10x HDF5 upload paths; document separate RNA, ATAC, and multiome routes.
- **Does not:** alter the original course README, Cell Ranger workflow, SLURM scripts, raw FASTQs, or upstream count generation.

For the primary 10x Gene Expression FASTQ-to-matrix workflow, continue to use the repository's existing [main README](README.md). This module begins **after** a filtered count matrix has been produced.
