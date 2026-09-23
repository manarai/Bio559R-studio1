#!/usr/bin/env bash
#SBATCH --job-name=cr_10x_A054
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --no-requeue
## If required by the cluster, uncomment and replace the following two lines.
## #SBATCH --account=<YOUR_ACCOUNT>
## #SBATCH --partition=<YOUR_PARTITION>

# Run Cell Ranger 9.0.1 on the BIO559R 10x_A054 mouse GEX test FASTQs.
# Submit with the command in BIO559R_10x_A054_SLURM_QUICKSTART.md so the
# log directory exists before Slurm opens standard output and standard error.

set -euo pipefail

# ---- Verified project inputs ------------------------------------------------
CELLRANGER_DIR="/home/terooatt/groups/fslg_dnasc/cellranger/cellranger-9.0.1"
CELLRANGER="${CELLRANGER_DIR}/cellranger"
MOUSE_REFERENCE="/home/terooatt/groups/fslg_dnasc/reference_genome/mouse/refdata-gex-GRCm39-2024-A"
FASTQ_DIR="/home/terooatt/groups/fslg_dnasc/nobackup/archive/BIO559R/Data/test/scRNA/10x_A054"
SAMPLE_ID="10x_A054"

# ---- Derived output paths ----------------------------------------------------
# This is deliberately outside the raw FASTQ directory.
RUN_PARENT="/home/terooatt/groups/fslg_dnasc/nobackup/archive/BIO559R/Data/test/scRNA/cellranger_runs"
RUN_ID="10x_A054_CR9_GRCm39_2024A"
LOCAL_MEM_GB=56   # Below the 64-GiB Slurm allocation; leaves memory for the OS.

mkdir -p "${RUN_PARENT}/provenance"
cd "$RUN_PARENT"

# ---- Read-only preflight checks ---------------------------------------------
[[ -x "$CELLRANGER" ]] || { echo "Cell Ranger executable not found: $CELLRANGER" >&2; exit 2; }
[[ -d "$MOUSE_REFERENCE" ]] || { echo "Mouse reference not found: $MOUSE_REFERENCE" >&2; exit 2; }
[[ -r "$FASTQ_DIR/${SAMPLE_ID}_S1_L003_R1_001.fastq.gz" ]] || { echo "R1 FASTQ is missing or unreadable" >&2; exit 2; }
[[ -r "$FASTQ_DIR/${SAMPLE_ID}_S1_L003_R2_001.fastq.gz" ]] || { echo "R2 FASTQ is missing or unreadable" >&2; exit 2; }

printf 'Started: %s\n' "$(date --iso-8601=seconds)" | tee "${RUN_PARENT}/provenance/${RUN_ID}.run_info.txt"
printf 'SLURM_JOB_ID=%s\n' "${SLURM_JOB_ID:-not_set}" | tee -a "${RUN_PARENT}/provenance/${RUN_ID}.run_info.txt"
"$CELLRANGER" --version | tee -a "${RUN_PARENT}/provenance/${RUN_ID}.run_info.txt"
printf 'Reference=%s\nFASTQ_DIR=%s\nSample=%s\n' "$MOUSE_REFERENCE" "$FASTQ_DIR" "$SAMPLE_ID" \
  | tee -a "${RUN_PARENT}/provenance/${RUN_ID}.run_info.txt"

# Cell Ranger uses exactly the CPU and memory limits granted by this Slurm job.
# --sample=10x_A054 matches 10x_A054_S1_L003_R1_001.fastq.gz and R2.
# The I1 FASTQ is an index read and is not an analysis input for cellranger count.
# --create-bam=false avoids a large aligned BAM for this teaching/test run.
"$CELLRANGER" count \
  --id="$RUN_ID" \
  --transcriptome="$MOUSE_REFERENCE" \
  --fastqs="$FASTQ_DIR" \
  --sample="$SAMPLE_ID" \
  --create-bam=false \
  --jobmode=local \
  --localcores="$SLURM_CPUS_PER_TASK" \
  --localmem="$LOCAL_MEM_GB"

printf 'Completed: %s\n' "$(date --iso-8601=seconds)" | tee -a "${RUN_PARENT}/provenance/${RUN_ID}.run_info.txt"
printf 'Review: %s/outs/web_summary.html\n' "$RUN_PARENT/$RUN_ID" \
  | tee -a "${RUN_PARENT}/provenance/${RUN_ID}.run_info.txt"
