#!/usr/bin/env bash
#SBATCH --job-name=check_10x_A054
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
## If required by the cluster, uncomment and replace the following two lines.
## #SBATCH --account=<YOUR_ACCOUNT>
## #SBATCH --partition=<YOUR_PARTITION>

# Lightweight integrity check for the raw FASTQs. This script reads only the
# three supplied .fastq.gz files and writes a small manifest in the run area.

set -euo pipefail

FASTQ_DIR="/home/terooatt/groups/fslg_dnasc/nobackup/archive/BIO559R/Data/test/scRNA/10x_A054"
SAMPLE_ID="10x_A054"
RUN_PARENT="/home/terooatt/groups/fslg_dnasc/nobackup/archive/BIO559R/Data/test/scRNA/cellranger_runs"
MANIFEST_DIR="${RUN_PARENT}/provenance"
mkdir -p "$MANIFEST_DIR"

R1="$FASTQ_DIR/${SAMPLE_ID}_S1_L003_R1_001.fastq.gz"
R2="$FASTQ_DIR/${SAMPLE_ID}_S1_L003_R2_001.fastq.gz"
I1="$FASTQ_DIR/${SAMPLE_ID}_S1_L003_I1_001.fastq.gz"

for file in "$R1" "$R2" "$I1"; do
  [[ -r "$file" ]] || { echo "Missing or unreadable: $file" >&2; exit 2; }
done

gzip -t "$R1" "$R2" "$I1"
ls -lh "$R1" "$R2" "$I1" | tee "$MANIFEST_DIR/${SAMPLE_ID}.fastq_inventory.txt"
sha256sum "$R1" "$R2" "$I1" > "$MANIFEST_DIR/${SAMPLE_ID}.fastq.sha256"
printf 'FASTQ gzip checks passed at %s\n' "$(date --iso-8601=seconds)" \
  | tee -a "$MANIFEST_DIR/${SAMPLE_ID}.fastq_inventory.txt"
