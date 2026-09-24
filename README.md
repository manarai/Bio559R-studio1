# Mouse scRNA-seq with Cell Ranger on a SLURM Supercomputer


**Scope:** A reproducible starter workflow for **10x Genomics mouse Gene Expression** FASTQs on a university SLURM cluster. The worked reference is `refdata-gex-GRCm39-2024-A` (mouse GRCm39). Always confirm that the selected reference is compatible with the installed Cell Ranger release before running production data.[1] [2]

> **BIO559R students:** For the supplied `10x_A054` test data, Cell Ranger 9.0.1 installation, and mouse GRCm39 reference, begin with [the course-specific quick start](BIO559R_10x_A054_SLURM_QUICKSTART.md) and its two SLURM scripts.

## 1. What this tutorial does

This workflow installs or loads the needed programs, checks the original demultiplexed FASTQ files, optionally performs **approved, Read-2-only** trimming, runs `cellranger count` in a SLURM allocation, and explains the main output files. It assumes that each library is a conventional 10x 3' or 5' Gene Expression library represented by paired FASTQs. It does **not** cover Cell Ranger ARC, V(D)J, Flex, or a multiplexed Feature Barcode experiment; those assays require different, version-specific commands.[3]

> **Important default:** Do **not** externally trim standard 10x Gene Expression FASTQs merely because FastQC reports a warning. In particular, do not trim, crop, or filter **Read 1 (R1)**. R1 contains cell-barcode and UMI sequence, and inconsistent R1 lengths can cause Cell Ranger failures. Cell Ranger already trims template-switch oligo and poly-A sequence from Gene Expression Read 2 before alignment.[4] [5]

Use external trimming only when the sequencing core and the protocol for the exact chemistry document a genuine R2 adapter or terminal-quality problem. Preserve the raw FASTQs unchanged, use **one** trimming program rather than both, rerun QC, and document the reason and command.

## 2. Before starting: cluster requirements and project layout

Ask the HPC team for the required SLURM account, partition, quality of service (QoS), maximum job count, software-module names, shared scratch location, backup/retention rules, and whether compute nodes may submit child jobs. Do not run Cell Ranger, FastQC, or MultiQC on a login node.

Cell Ranger's published general baseline is at least 8 AVX-capable CPU cores, 64 GB RAM, and 1.5 TB of free disk. These are **not** a universal allocation. The required resources depend on read count, cells, assay, whether a BAM is created, storage performance, and the Cell Ranger version. Start with one representative library and use SLURM accounting to set later requests.[6] [7]

Create a project directory on the site-approved high-performance filesystem. Keep raw data read-only and separate from derived files.

```bash
export PROJECT=/path/on/scratch_or_project/mouse_scrnaseq_2026
mkdir -p "$PROJECT"/{raw_fastqs,reference,software,qc,trimmed_fastqs,runs,logs,provenance}

# Copy or link, but do not alter, demultiplexed input FASTQs.
# Expected Cell Ranger-style names include:
# SampleA_S1_L001_R1_001.fastq.gz
# SampleA_S1_L001_R2_001.fastq.gz
```

A Cell Ranger FASTQ directory may contain multiple lanes for the same library. Use `--fastqs` to point to the directory and `--sample` to select the Illumina sample prefix; do not concatenate or casually rename lane FASTQs.[8]

Record the software version, reference name, FASTQ paths, checksums when practical, command line, and SLURM job ID in `$PROJECT/provenance/`.

## 3. Install or load the required software

### 3.1 Preferred method: use site-managed modules

First see whether the HPC center already provides supported versions. Module names vary by cluster, so the following is a discovery step rather than a copy-and-paste command.

```bash
module spider cellranger
module spider fastqc
module spider multiqc
module spider cutadapt
module spider trimmomatic
module spider samtools
```

Load the exact versions selected for the project and save them in the run log.

```bash
module load CellRanger/<VERSION>
module load FastQC/<VERSION>
module load MultiQC/<VERSION>
module load Cutadapt/<VERSION>
module load Trimmomatic/<VERSION>
module load SAMtools/<VERSION>

{ module list; cellranger --version; fastqc --version; multiqc --version; \
  cutadapt --version; samtools --version; } 2>&1 | tee "$PROJECT/provenance/software_versions.txt"
```

### 3.2 User-space alternative for FastQC, MultiQC, Cutadapt, Trimmomatic, and SAMtools

If the cluster provides Miniconda, Miniforge, Conda, or micromamba, create one isolated environment for the lightweight QC tools. Ask the HPC team before installing into a shared location. Cell Ranger itself should be installed from the official 10x archive, as described in the next section.

```bash
# Example only: replace this with the site-approved Conda or micromamba module.
module load Miniforge3/<VERSION>

conda create -y -p "$PROJECT/software/scrna-qc" \
  -c conda-forge -c bioconda \
  fastqc multiqc cutadapt trimmomatic samtools
conda activate "$PROJECT/software/scrna-qc"

fastqc --version
multiqc --version
cutadapt --version
samtools --version
# Trimmomatic is launched with Java; locate the JAR after installation if needed.
```

| Tool | Role in this workflow | Required for `cellranger count`? |
|---|---|---|
| Cell Ranger | 10x alignment, barcode/UMI processing, cell calling, feature-by-barcode matrices, and web summary | Yes |
| FastQC | Per-FASTQ diagnostic report | No, but strongly recommended |
| MultiQC | Combined report across FastQC and other logs | No, but strongly recommended |
| Cutadapt | Optional paired-end adapter or quality trimming | No; use only when justified |
| Trimmomatic | Alternative optional trimming program | No; use only when justified |
| SAMtools | Optional post-run BAM integrity and alignment summary checks | No; only if a BAM was created |

## 4. Install Cell Ranger and download the mouse reference

### 4.1 Download and install Cell Ranger

Choose a Cell Ranger release from the [official 10x Genomics downloads page][1]. License acceptance or authenticated download may be required. Download the Linux archive in a browser and transfer it to the cluster with the institution-approved method, or download from a compute-accessible transfer node if permitted.

```bash
# Example assumes the archive has been transferred to the shared application directory.
cd /shared/apps                              # Replace with the approved site location.
tar -xzf cellranger-<CELLRANGER_VERSION>.tar.gz

export PATH="/shared/apps/cellranger-<CELLRANGER_VERSION>:$PATH"
cellranger --version

# Run once in a writable test directory before analyzing real data.
cd "$PROJECT/runs"
cellranger testrun --id=cellranger_testrun
```

A shared installation must be readable at the **same absolute path** from every node if cluster mode will be used. Pin one Cell Ranger version per project. Do not use the GRCm39-2024-A reference with Cell Ranger v5.0.1 or older; use the reference compatibility notes or, preferably, upgrade Cell Ranger.[2]

### 4.2 Download and verify the pre-built mouse reference

For a conventional mouse Gene Expression experiment, download the 10x pre-built GRCm39 reference. The archive below is `refdata-gex-GRCm39-2024-A`; at the time of writing, 10x publishes MD5 `37c51137ccaeabd4d151f80dc86ce0b3` for this file. Confirm the checksum against the download page before use.[1] [9]

```bash
cd "$PROJECT/reference"
curl -fLO https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz

echo '37c51137ccaeabd4d151f80dc86ce0b3  refdata-gex-GRCm39-2024-A.tar.gz' | md5sum -c -
tar -xzf refdata-gex-GRCm39-2024-A.tar.gz

export MOUSE_REF="$PROJECT/reference/refdata-gex-GRCm39-2024-A"
test -f "$MOUSE_REF/fasta/genome.fa" && echo "Mouse reference is present"
```

The 2024-A package uses the **GRCm39** assembly. Do not merge its matrices or gene identifiers with a legacy **mm10/GRCm38** reference without a documented harmonization procedure.[2]

## 5. Inspect raw FASTQs before any modification

Use FastQC and MultiQC on original files first. The FastQC reports are diagnostics, not automatic pass/fail criteria for barcoded single-cell libraries. Review R1 and R2 separately, compare lanes, and look for a consistent issue across a library or run. R1's fixed barcode/UMI structure can produce composition or duplication warnings that are expected for this assay.[4] [10] [11]

Save the following as `01_fastq_qc.sbatch`. The resource values are deliberately modest **examples for QC only**. Replace account, partition, time, CPUs, and memory with values approved by the HPC team.

```bash
#!/usr/bin/env bash
#SBATCH --job-name=fastq_qc
#SBATCH --account=<ACCOUNT>
#SBATCH --partition=<PARTITION>
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/path/to/project/logs/%x-%j.out
#SBATCH --error=/path/to/project/logs/%x-%j.err

set -euo pipefail

PROJECT=/path/to/project
FASTQ_DIR="$PROJECT/raw_fastqs"
QC_DIR="$PROJECT/qc/raw"
mkdir -p "$QC_DIR/fastqc" "$QC_DIR/multiqc"

module load FastQC/<VERSION> MultiQC/<VERSION>  # Or activate the Conda environment.

mapfile -t FASTQS < <(find "$FASTQ_DIR" -maxdepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) | sort)
(( ${#FASTQS[@]} > 0 )) || { echo "No gzipped FASTQs found" >&2; exit 2; }

gzip -t "${FASTQS[@]}"                         # Confirms gzip streams are readable.
fastqc --threads "$SLURM_CPUS_PER_TASK" --outdir "$QC_DIR/fastqc" "${FASTQS[@]}"
multiqc --force --outdir "$QC_DIR/multiqc" --filename multiqc_raw.html "$QC_DIR/fastqc"

printf '%s\n' "${FASTQS[@]}" > "$PROJECT/provenance/raw_fastq_paths.txt"
sha256sum "${FASTQS[@]}" > "$PROJECT/provenance/raw_fastq.sha256"
```

Submit it with `sbatch 01_fastq_qc.sbatch`. Open `qc/raw/multiqc/multiqc_raw.html` in a browser after transferring a copy from the cluster, or use an institution-approved remote browser. Review the complete report before deciding whether trimming is warranted.

## 6. Optional FASTQ trimming: decision rule and safe methods

### 6.1 Decision rule for a 10x library

External trimming is usually unnecessary for 10x Gene Expression. If FastQC identifies a suspected adapter or low-quality tail, first ask the sequencing core whether the issue is real, whether it is in R2, and which adapter sequence and chemistry apply. A universal Illumina adapter sequence, Phred threshold, or minimum length is not scientifically defensible for all 10x libraries.[4] [5]

If trimming is approved, follow all of these rules:

1. Copy outputs to `$PROJECT/trimmed_fastqs`; do not overwrite raw data.
2. Keep **R1 unmodified**. Preserve its bases and length.
3. Preserve synchronized R1/R2 pairs and original Cell Ranger-compatible sample, lane, and read naming.
4. Run FastQC and MultiQC again on the paired trimmed outputs.
5. Give Cell Ranger only the paired outputs. Never give it Trimmomatic's unpaired reads.
6. Use either Cutadapt **or** Trimmomatic for a given derived dataset, not both.

### 6.2 Cutadapt: preferred controlled R2-only trimming

Cutadapt can target R2 independently in paired-end mode. In the template below, `-A` removes a confirmed **3' adapter from R2**, `-Q 20` applies a 3' quality cutoff to R2, `-q 0` leaves R1 quality trimming off, and `-m :20` retains only pairs with an R2 at least 20 bases long. The adapter sequence and cutoff are placeholders: obtain them from the sequencing core and exact library documentation.[12]

```bash
# Run inside a scheduled SLURM job. Do not replace placeholders with generic values.
module load Cutadapt/<VERSION>

R1_IN=/path/raw_fastqs/SampleA_S1_L001_R1_001.fastq.gz
R2_IN=/path/raw_fastqs/SampleA_S1_L001_R2_001.fastq.gz
OUT_DIR=/path/trimmed_fastqs
mkdir -p "$OUT_DIR"

R1_OUT="$OUT_DIR/SampleA_S1_L001_R1_001.fastq.gz"
R2_OUT="$OUT_DIR/SampleA_S1_L001_R2_001.fastq.gz"

cutadapt --cores="$SLURM_CPUS_PER_TASK" \
  -A '<CONFIRMED_R2_3PRIME_ADAPTER_SEQUENCE>' \
  -q 0 -Q <R2_3PRIME_PHRED_CUTOFF> -m :<R2_MINIMUM_LENGTH> \
  -o "$R1_OUT" -p "$R2_OUT" "$R1_IN" "$R2_IN" \
  > "$OUT_DIR/SampleA_S1_L001.cutadapt.log" 2>&1

gzip -t "$R1_OUT" "$R2_OUT"
```

For a confirmed R2 quality issue without an adapter, omit `-A`. For a confirmed adapter issue without quality trimming, omit `-Q`. Repeat the operation for every matching lane/chunk pair; do not concatenate files before Cell Ranger.

### 6.3 Trimmomatic: an alternative, with a 10x-specific limitation

Trimmomatic's `LEADING`, `TRAILING`, and `SLIDINGWINDOW` operations are valuable for ordinary paired-end data because they remove low-quality ends or windows. However, these quality operations affect both mates in a standard paired-end command. Therefore, **do not apply a generic Trimmomatic quality-trimming command to a standard 10x Gene Expression R1/R2 pair**: it can alter R1 barcode/UMI sequence.[13]

For an approved **R2 adapter-only** intervention, create an adapter FASTA containing the core-confirmed R2 adapter sequence with a FASTA header ending in `/2`. The `/2` suffix limits adapter matching to the second read. This command deliberately omits `LEADING`, `TRAILING`, `SLIDINGWINDOW`, `CROP`, `HEADCROP`, and `MINLEN` to protect R1.

```bash
# r2_adapters.fa example. Replace the sequence with the confirmed sequence.
# >confirmed_adapter/2
# <CONFIRMED_R2_3PRIME_ADAPTER_SEQUENCE>

module load Java/<VERSION> Trimmomatic/<VERSION>

java -jar <PATH_TO_TRIMMOMATIC_JAR> PE -threads "$SLURM_CPUS_PER_TASK" -phred33 -validatePairs \
  -summary "$OUT_DIR/SampleA_S1_L001.trimmomatic.summary.txt" \
  "$R1_IN" "$R2_IN" \
  "$R1_OUT" "$OUT_DIR/SampleA_S1_L001_R1_unpaired.fastq.gz" \
  "$R2_OUT" "$OUT_DIR/SampleA_S1_L001_R2_unpaired.fastq.gz" \
  "ILLUMINACLIP:/path/r2_adapters.fa:2:30:10" \
  > "$OUT_DIR/SampleA_S1_L001.trimmomatic.log" 2>&1
```

> **For non-10x or explicitly approved non-barcode data only:** a conventional Trimmomatic quality-trimming pattern is `LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36`. Do not apply that pattern to standard 10x GEX R1/R2 FASTQs without written assay-specific guidance.

After either trimming method, rerun the QC job using `$PROJECT/trimmed_fastqs` as `FASTQ_DIR` and compare the new MultiQC report with the raw report.

## 7. Run `cellranger count` in a SLURM job

### 7.1 Recommended starting point: one Cell Ranger run within one allocated node

For most users, start with Cell Ranger's **job-submission (local) mode** in one allocated compute-node job. Set `--ntasks=1`, set `--cpus-per-task` to the number of cores planned for Cell Ranger, set `--localcores` to that same number, and set `--localmem` to a value below the requested node memory. Tenx documentation illustrates using about 90% of the SLURM memory allocation for `--localmem`.[7]

Save this as `02_cellranger_count.sbatch`. Replace every angle-bracket placeholder. The illustrative 16 CPUs and 128 GB are **not a required recommendation**; select values from a pilot and local HPC policy.

```bash
#!/usr/bin/env bash
#SBATCH --job-name=mouse_gex
#SBATCH --account=<ACCOUNT>
#SBATCH --partition=<PARTITION>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=/path/to/project/logs/%x-%j.out
#SBATCH --error=/path/to/project/logs/%x-%j.err
#SBATCH --no-requeue

set -euo pipefail

PROJECT=/path/to/project
CELLRANGER=/shared/apps/cellranger-<CELLRANGER_VERSION>/cellranger
MOUSE_REF="$PROJECT/reference/refdata-gex-GRCm39-2024-A"
FASTQ_DIR="$PROJECT/raw_fastqs"              # Or $PROJECT/trimmed_fastqs if justified and QCed.
SAMPLE_PREFIX=SampleA                           # Prefix before _S1_L001_R1_001.fastq.gz.
RUN_ID=SampleA_GEX_GRCm39_2024A
CREATE_BAM=false                                # Use true only when a BAM is needed.
LOCAL_MEM_GB=115                                # Must be lower than requested --mem=128G.

mkdir -p "$PROJECT/logs" "$PROJECT/runs"
cd "$PROJECT/runs"

"$CELLRANGER" count \
  --id="$RUN_ID" \
  --transcriptome="$MOUSE_REF" \
  --fastqs="$FASTQ_DIR" \
  --sample="$SAMPLE_PREFIX" \
  --create-bam="$CREATE_BAM" \
  --jobmode=local \
  --localcores="$SLURM_CPUS_PER_TASK" \
  --localmem="$LOCAL_MEM_GB"

"$CELLRANGER" --version | tee "$PROJECT/provenance/${RUN_ID}.cellranger_version.txt"
printf 'SLURM_JOB_ID=%s\n' "$SLURM_JOB_ID" | tee "$PROJECT/provenance/${RUN_ID}.slurm_job.txt"
```

Submit with:

```bash
sbatch 02_cellranger_count.sbatch
squeue -u "$USER"
sacct -j <JOB_ID> --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS
```

The run is successful only when Cell Ranger reports `Pipestance completed successfully!`. The `--id` must be unique; an existing ID directory is treated as a resumable pipestance. If a run fails, preserve its logs, correct the cause, and rerun the original command to resume only after confirming no instance remains active.[3] [14]

For several independent GEM wells, submit one job per library. A SLURM array can be useful, but throttle it to the HPC-approved concurrent-job limit, for example `#SBATCH --array=1-<N>%<ALLOWED_CONCURRENT_JOBS>`. Do not use one `cellranger count` run to combine biologically independent libraries.

### 7.2 Advanced mode: Cell Ranger submits many SLURM stage jobs

Use Cell Ranger cluster mode only with an administrator-tested Martian SLURM template. It can submit hundreds or thousands of jobs. The Cell Ranger installation and project filesystem must be shared and visible at the same paths on all nodes; the site must allow the parent process to run and submit `sbatch` child jobs.[15]

Start with the `slurm.template.example` supplied with the exact installed Cell Ranger version. The administrator should adapt and validate directives while retaining Martian placeholders such as `__MRO_JOB_NAME__`, `__MRO_THREADS__`, `__MRO_MEM_GB__`, `__MRO_STDOUT__`, `__MRO_STDERR__`, and `__MRO_CMD__`.

```bash
# Example only after HPC validation of /shared/apps/cellranger-slurm.template
cellranger count \
  --id=SampleA_GEX_GRCm39_2024A \
  --transcriptome="$MOUSE_REF" \
  --fastqs="$FASTQ_DIR" \
  --sample=SampleA \
  --create-bam=false \
  --jobmode=/shared/apps/cellranger-slurm.template
```

Do not increase `--maxjobs` or decrease `--jobinterval` without HPC approval. If nested job submission is disallowed, use the single-node local-mode script in Section 7.1 instead.[15]

## 8. Inspect the completed Cell Ranger output

A successful run creates `$PROJECT/runs/<RUN_ID>/outs/`. Start with `web_summary.html`; it provides sample-level metrics and alerts. Treat the automated clustering and dimensionality reduction as exploratory QC, not as a final biological conclusion.[16]

| Output | What it contains | Typical next use |
|---|---|---|
| `web_summary.html` | Interactive QC summary and alerts | First review of sequencing, mapping, cell, and library metrics |
| `metrics_summary.csv` | Machine-readable summary metrics | Project QC spreadsheet or reproducibility record |
| `filtered_feature_bc_matrix/` and `.h5` | UMI matrix for Cell Ranger-called cell barcodes | Standard input to Seurat, Scanpy, or other downstream analysis |
| `raw_feature_bc_matrix/` and `.h5` | Matrix for all detected barcodes, including ambient/background barcodes | Specialized methods; not a direct substitute for filtered data |
| `analysis/` | Cell Ranger PCA, clusters, UMAP/t-SNE, and differential-expression outputs | Initial visualization and QC only |
| `possorted_genome_bam.bam` and `.bai` | Aligned reads, only when `--create-bam=true` | Alignment review or BAM-based downstream tools |
| `molecule_info.h5` | Molecule-level information | Retain if later using `cellranger aggr` |
| `cloupe.cloupe` | Loupe Browser file | Optional interactive exploration in Loupe Browser |

The main metrics to review are the estimated number of cells, mean reads per cell, median genes per cell, median UMIs per cell, fraction of reads in cells, sequencing saturation, valid barcodes, and confident transcriptome mapping. Interpret these metrics against the experiment's cell-recovery target, chemistry, tissue quality, and comparable libraries—not against a single universal cutoff.[16]

If a BAM was produced, perform optional integrity and summary checks on a compute node:

```bash
module load SAMtools/<VERSION>
OUTS="$PROJECT/runs/$RUN_ID/outs"
samtools quickcheck -v "$OUTS/possorted_genome_bam.bam"
samtools flagstat -@ "$SLURM_CPUS_PER_TASK" "$OUTS/possorted_genome_bam.bam" \
  > "$PROJECT/qc/${RUN_ID}.flagstat.txt"
```

`samtools quickcheck` checks the BAM header and EOF block; it is not a complete corruption test. SAMtools is not needed to create the Cell Ranger matrix and should not be used to regenerate the original 10x FASTQs for initial analysis.[17] [18]

## 9. Troubleshooting checklist

| Symptom | First action |
|---|---|
| `FASTQ` files are not found | Confirm the directory, `--sample` prefix, lane naming, and that every R1 has its matching R2. |
| Mixed or short R1 length error | Return to original demultiplexed FASTQs. Check for trimming during demultiplexing or manual preprocessing; do not repair R1 by ad hoc cropping. [5] |
| SLURM memory or time-limit termination | Inspect `sacct` for `State`, `ExitCode`, `ReqMem`, and `MaxRSS`; increase only the constraint that failed and rerun the same invocation if the pipestance is inactive. [7] [14] |
| Slow run or lack of disk space | Verify scratch capacity, I/O policy, optional BAM choice, and filesystem quotas before resubmission. |
| Cluster mode jobs stay pending or fail to submit | Check local account/QoS/pending reason codes and ask the HPC team to validate the Martian template. Do not raise job throttles independently. [15] [19] |
| FastQC warning looks severe | Compare R1/R2 and all lanes, verify the actual chemistry and sequencing configuration with the core, and inspect Cell Ranger web-summary metrics before trimming. |

For a failed Cell Ranger run, preserve `<RUN_ID>/log`, scheduler stdout/stderr, and `<RUN_ID>.mri.tgz` before changing anything. These are the materials required for local troubleshooting or a 10x support request.[14]

## 10. Recommended student learning links

The following sources explain the workflow in greater depth. Begin with the Cell Ranger download, `count` pipeline, system-requirements, and output pages; then consult the tool manuals before performing optional trimming.

1. [Cell Ranger downloads and pre-built references][1]
2. [Cell Ranger `count` pipeline][3]
3. [Cell Ranger job-submission mode on SLURM][7]
4. [Cell Ranger cluster mode and templates][15]
5. [Cell Ranger Gene Expression outputs][16]
6. [FastQC documentation][10]
7. [MultiQC documentation][11]
8. [Cutadapt quality trimming and paired-end guide][12]
9. [Trimmomatic documentation][13]
10. [SAMtools manuals][17]

## References

[1]: https://www.10xgenomics.com/support/software/cell-ranger/downloads "10x Genomics: Cell Ranger downloads and pre-built references"
[2]: https://www.10xgenomics.com/support/software/cell-ranger/10.0/release-notes/cr-reference-release-notes "10x Genomics: Cell Ranger reference release notes and compatibility"
[3]: https://www.10xgenomics.com/support/software/cell-ranger/10.0/analysis/running-pipelines/cr-gex-count "10x Genomics: The cellranger count pipeline"
[4]: https://www.10xgenomics.com/support/software/cell-ranger/10.0/algorithms-overview/cr-gex-algorithm "10x Genomics: Gene Expression algorithm and Read 2 trimming"
[5]: https://www.10xgenomics.com/support/software/cell-ranger/10.0/resources/cr-error-codes "10x Genomics: Cell Ranger error codes including mixed R1 lengths"
[6]: https://www.10xgenomics.com/support/software/cell-ranger/10.0/getting-started/cr-system-requirements "10x Genomics: Cell Ranger system requirements"
[7]: https://www.10xgenomics.com/support/software/cell-ranger/10.0/advanced/cr-job-submission-mode "10x Genomics: Cell Ranger job submission mode"
[8]: https://www.10xgenomics.com/support/software/cell-ranger/8.0/analysis/cr-specifying-fastqs "10x Genomics: Specifying input FASTQ files"
[9]: https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz "10x Genomics: refdata-gex-GRCm39-2024-A reference archive"
[10]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ "FastQC: A quality-control tool for high-throughput sequence data"
[11]: https://docs.seqera.io/multiqc/getting_started/running_multiqc "MultiQC documentation: Running MultiQC"
[12]: https://cutadapt.readthedocs.io/en/stable/guide.html "Cutadapt user guide: quality trimming and paired-end processing"
[13]: https://github.com/usadellab/Trimmomatic "Trimmomatic: flexible read trimming for Illumina data"
[14]: https://www.10xgenomics.com/support/software/cell-ranger/10.0/advanced/cr-troubleshooting "10x Genomics: Cell Ranger troubleshooting"
[15]: https://www.10xgenomics.com/support/software/cell-ranger/10.0/advanced/cr-cluster-mode "10x Genomics: Cell Ranger cluster mode"
[16]: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-gex-overview "10x Genomics: Gene Expression output files"
[17]: https://www.htslib.org/doc/samtools-quickcheck.html "SAMtools quickcheck manual"
[18]: https://www.htslib.org/doc/samtools-fasta.html "SAMtools fastq/fasta manual"
[19]: https://slurm.schedmd.com/job_reason_codes.html "Slurm Workload Manager: job reason codes"
