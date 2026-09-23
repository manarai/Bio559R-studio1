# BIO559R quick start: run `10x_A054` with Cell Ranger 9.0.1 on SLURM

This is a **student-ready test run** for the supplied 10x mouse scRNA-seq FASTQs. It uses the paths verified from the supplied terminal screenshots and is deliberately configured as a **single-node SLURM job**. It does not modify the raw FASTQs.

> **Do not run Cell Ranger on `login01`.** Submit the batch file with `sbatch`; SLURM will run it on a compute node.

## Verified inputs

| Item | Verified path or value |
|---|---|
| Paired FASTQ folder | `/home/terooatt/groups/fslg_dnasc/nobackup/archive/BIO559R/Data/test/scRNA/10x_A054` |
| R1 input | `10x_A054_S1_L003_R1_001.fastq.gz` |
| R2 input | `10x_A054_S1_L003_R2_001.fastq.gz` |
| Index FASTQ present | `10x_A054_S1_L003_I1_001.fastq.gz` |
| Cell Ranger used | `/home/terooatt/groups/fslg_dnasc/cellranger/cellranger-9.0.1/cellranger` |
| Mouse reference used | `/home/terooatt/groups/fslg_dnasc/reference_genome/mouse/refdata-gex-GRCm39-2024-A` |
| Cell Ranger sample selector | `10x_A054` |

The screenshot also shows `cellranger-3.0.2` and the legacy `refdata-gex-mm10-2020-A` reference. This test intentionally uses **Cell Ranger 9.0.1** with the newer **GRCm39-2024-A** mouse reference. Do not combine output matrices generated with GRCm39 and the legacy mm10 reference without a documented harmonization procedure.[1]

The `I1` file is an index-read FASTQ generated during demultiplexing. It is checked for file integrity but is not supplied as a separate input argument to `cellranger count`; Cell Ranger discovers the matching R1/R2 files from `--fastqs` and `--sample`.[2]

## 1. Copy the two batch scripts to the cluster

The two scripts supplied with this tutorial are:

- `01_check_fastqs_10x_A054.sbatch`, which confirms that all three gzip-compressed FASTQs are readable and records checksums.
- `02_cellranger_count_10x_A054.sbatch`, which runs `cellranger count` on the R1/R2 pair.

Copy the scripts to a directory on the supercomputer, for example:

```bash
mkdir -p ~/BIO559R_10x_A054_scripts
# Transfer the two .sbatch files to this directory using the course-approved method.
cd ~/BIO559R_10x_A054_scripts
```

Alternatively, create the files directly on the cluster by copying their contents from this tutorial package. Do not copy the raw FASTQs; the scripts use the existing shared paths.

## 2. Set up a writable run directory and submit the integrity check

The scripts are configured to write to:

```text
/home/terooatt/groups/fslg_dnasc/nobackup/archive/BIO559R/Data/test/scRNA/cellranger_runs
```

Before submitting, check that the directory is writable. If it is not writable, edit **only** the `RUN_PARENT` line in both scripts to a course-approved writable project or scratch directory. Never set `RUN_PARENT` inside `10x_A054`, because that directory contains the original source FASTQs.

```bash
export RUN_PARENT=/home/terooatt/groups/fslg_dnasc/nobackup/archive/BIO559R/Data/test/scRNA/cellranger_runs
mkdir -p "$RUN_PARENT"/{logs,provenance}
test -w "$RUN_PARENT" && echo "Run directory is writable"
```

The scripts contain generic SLURM requests. If the supercomputer requires an account, partition, or QoS, uncomment and complete the relevant `#SBATCH` lines in both scripts before submission. Do **not** invent a partition or account name.

Submit the file check. Explicit `--output` and `--error` arguments ensure that the log directory already exists before SLURM opens the log files.

```bash
cd ~/BIO559R_10x_A054_scripts
CHECK_JOB=$(sbatch \
  --output="$RUN_PARENT/logs/check_10x_A054-%j.out" \
  --error="$RUN_PARENT/logs/check_10x_A054-%j.err" \
  01_check_fastqs_10x_A054.sbatch | awk '{print $4}')
echo "FASTQ check job: $CHECK_JOB"

squeue -j "$CHECK_JOB"
```

After the job finishes, confirm that it passed:

```bash
sacct -j "$CHECK_JOB" --format=JobID,JobName,State,ExitCode,Elapsed
cat "$RUN_PARENT/provenance/10x_A054.fastq_inventory.txt"
```

The expected success message is `FASTQ gzip checks passed`. A failure means the user should stop and report the exact scheduler log and missing/unreadable filename; do not run Cell Ranger on incomplete FASTQs.

## 3. Submit the Cell Ranger test run

The Cell Ranger script requests one node, one task, 8 CPU cores, 64 GB RAM, and one day. These values meet the published general minimum for Cell Ranger but must comply with the local cluster policy and may need to increase for a different dataset. The script caps Cell Ranger at the job allocation: `--localcores=$SLURM_CPUS_PER_TASK` and `--localmem=56` GB. This prevents it from using more CPUs or memory than SLURM granted.[3] [4]

The test run intentionally uses `--create-bam=false` to reduce disk use and elapsed time. The expected outputs still include the Cell Ranger HTML web summary, metrics CSV, filtered feature-by-barcode matrix, raw matrix, and Cell Ranger secondary analysis. A BAM is not needed for routine Seurat or Scanpy matrix-based analysis.[5]

```bash
cd ~/BIO559R_10x_A054_scripts
COUNT_JOB=$(sbatch \
  --output="$RUN_PARENT/logs/cr_10x_A054-%j.out" \
  --error="$RUN_PARENT/logs/cr_10x_A054-%j.err" \
  02_cellranger_count_10x_A054.sbatch | awk '{print $4}')
echo "Cell Ranger job: $COUNT_JOB"

squeue -j "$COUNT_JOB"
```

The `--sample=10x_A054` argument is correct because it matches the beginning of both files:

```text
10x_A054_S1_L003_R1_001.fastq.gz
10x_A054_S1_L003_R2_001.fastq.gz
```

The script supplies the FASTQ **directory**, not individual R1/R2 filenames. This lets Cell Ranger collect all matching lanes for that sample if more lanes are added later.[2]

## 4. Monitor the run and inspect results

While the job is pending or running:

```bash
squeue -j "$COUNT_JOB" -o '%.18i %.20j %.2t %.10M %.6D %R'
tail -f "$RUN_PARENT/logs/cr_10x_A054-${COUNT_JOB}.out"
```

After it is complete:

```bash
sacct -j "$COUNT_JOB" --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS

export RUN_ID=10x_A054_CR9_GRCm39_2024A
ls -lh "$RUN_PARENT/$RUN_ID/outs"
```

A successful run ends with **`Pipestance completed successfully!`** in the Cell Ranger output. The first file to review is:

```text
$RUN_PARENT/10x_A054_CR9_GRCm39_2024A/outs/web_summary.html
```

Copy this HTML file to a local computer or open it using a course-approved browser method. Then review the following outputs:

| Output | Use |
|---|---|
| `outs/web_summary.html` | Interactive Cell Ranger QC summary and alerts; review first. |
| `outs/metrics_summary.csv` | Tabular metrics for a lab notebook or class submission. |
| `outs/filtered_feature_bc_matrix/` or `.h5` | Standard input for downstream Seurat or Scanpy analysis. |
| `outs/raw_feature_bc_matrix/` or `.h5` | Includes background/uncalled barcodes; do not substitute it for the filtered matrix without a specific method. |
| `outs/analysis/` | Cell Ranger PCA, clusters, UMAP/t-SNE, and exploratory differential-expression files. |

## 5. What to do if it fails

Do not delete the run directory and do not submit a second job with the same `RUN_ID` while the first job is active. Instead, collect the following information:

```bash
sacct -j "$COUNT_JOB" --format=JobID,State,ExitCode,Elapsed,ReqMem,MaxRSS
cat "$RUN_PARENT/logs/cr_10x_A054-${COUNT_JOB}.err"
find "$RUN_PARENT/10x_A054_CR9_GRCm39_2024A/log" -type f -maxdepth 2 -print
```

A memory or time-limit error should be reported with `sacct` output to the instructor or HPC team. A FASTQ, reference, or permission error should be reported with the relevant error log. After the cause is fixed and no pipestance is running, resubmit the same command to allow Cell Ranger to resume.[6]

## References

[1]: https://www.10xgenomics.com/support/software/cell-ranger/10.0/release-notes/cr-reference-release-notes "10x Genomics: Cell Ranger reference release notes and compatibility"
[2]: https://www.10xgenomics.com/support/software/cell-ranger/8.0/analysis/cr-specifying-fastqs "10x Genomics: Specifying input FASTQ files"
[3]: https://www.10xgenomics.com/support/software/cell-ranger/10.0/getting-started/cr-system-requirements "10x Genomics: Cell Ranger system requirements"
[4]: https://www.10xgenomics.com/support/software/cell-ranger/10.0/advanced/cr-job-submission-mode "10x Genomics: Cell Ranger job submission mode"
[5]: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-gex-overview "10x Genomics: Gene Expression output files"
[6]: https://www.10xgenomics.com/support/software/cell-ranger/10.0/advanced/cr-troubleshooting "10x Genomics: Cell Ranger troubleshooting"
