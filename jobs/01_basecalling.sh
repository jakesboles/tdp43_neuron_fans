#!/bin/bash
#SBATCH -A p32095
#SBATCH -p gengpu
#SBATCH --gres=gpu:h100:1
#SBATCH -t 48:00:00
#SBATCH --mem=72G
#SBATCH --array=1-12
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --job-name=01_basecalling
#SBATCH -o /projects/b1042/Gate_Lab/boles/tdp43_neuron_fans/logs/%x_%A_%a.log

# Basecall one sequencing run (one folder under pod5/) per array task.
# Samples are read from config/samples.tsv (one row per pod5/ folder), so the
# array size must equal the number of rows in that file (currently 12).
# Pools are combined later, in 02_epi2me_single_cell.sh.
#
# NOTE: the logs/ directory must exist BEFORE sbatch, or SLURM can't write the
# log and the job dies silently:
#   mkdir -p /projects/b1042/Gate_Lab/boles/tdp43_neuron_fans/logs

set -euo pipefail

PROJECT_DIR="/projects/b1042/Gate_Lab/boles/tdp43_neuron_fans"
SAMPLE_SHEET="${PROJECT_DIR}/config/samples.tsv"
POD5_ROOT="${PROJECT_DIR}/pod5"
OUTPUT_DIR="${PROJECT_DIR}/basecalled"
DORADO="/projects/b1169/jackie/software/dorado-0.9.0-linux-x64/bin/dorado"
DORADO_MODEL="/projects/b1169/benney/ont_ref/dorado_model/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"

cd "${PROJECT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# +1 skips the header row of the sample sheet
SAMPLE=$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR == n { print $1 }' "${SAMPLE_SHEET}")
if [[ -z "${SAMPLE}" ]]; then
  echo "ERROR: no row ${SLURM_ARRAY_TASK_ID} in ${SAMPLE_SHEET}" >&2
  exit 1
fi

# Each sample folder holds one MinKNOW run folder, which holds pod5/.
# Find it rather than assuming exactly one entry from `ls`.
mapfile -t POD5_DIRS < <(find "${POD5_ROOT}/${SAMPLE}" -mindepth 2 -maxdepth 2 -type d -name pod5)
if [[ ${#POD5_DIRS[@]} -ne 1 ]]; then
  echo "ERROR: expected exactly 1 pod5/ directory under ${POD5_ROOT}/${SAMPLE}, found ${#POD5_DIRS[@]}" >&2
  printf '  %s\n' "${POD5_DIRS[@]}" >&2
  exit 1
fi
POD5_DIR="${POD5_DIRS[0]}"

N_POD5=$(find "${POD5_DIR}" -maxdepth 1 -name '*.pod5' | wc -l)
if [[ ${N_POD5} -eq 0 ]]; then
  echo "ERROR: no .pod5 files in ${POD5_DIR} (transfer incomplete?)" >&2
  exit 1
fi

OUT_FASTQ="${OUTPUT_DIR}/${SAMPLE}.fastq.gz"
TMP_FASTQ="${OUT_FASTQ}.partial"

if command -v pigz >/dev/null 2>&1; then
  COMPRESS="pigz -p 8"
else
  COMPRESS="gzip"
fi

echo "Sample           = ${SAMPLE}"
echo "Dorado           = ${DORADO}"
${DORADO} --version 2>&1 || true
echo "Dorado model     = ${DORADO_MODEL}"
echo "POD5 directory   = ${POD5_DIR} (${N_POD5} files)"
echo "Output           = ${OUT_FASTQ}"

# --no-trim: wf-single-cell locates the 10x adapter / barcode / UMI / polyT
# itself and trims them after extraction. Leave the reads untouched here so
# dorado's adapter/primer trimming can't clip sequence that barcode
# identification relies on.
${DORADO} basecaller "${DORADO_MODEL}" "${POD5_DIR}" \
  --min-qscore 10 \
  --no-trim \
  --emit-fastq \
  | ${COMPRESS} > "${TMP_FASTQ}"

# Only expose the final file once basecalling + compression finished cleanly,
# so a failed/timed-out job never leaves a truncated FASTQ that 02 would use.
mv "${TMP_FASTQ}" "${OUT_FASTQ}"
echo "Done: ${OUT_FASTQ}"
