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
#SBATCH -o /projects/b1169/boles/tdp43_neuron_fans/logs/%x_%A_%a.log

# Basecall one sequencing run (one folder under pod5/) per array task.
# Runs are read from config/samples.tsv: array task N = data row N, so the
# --array range must match the number of rows (currently 12). To redo or add
# specific runs, e.g. the 158-A/158-B runs (rows 5-12):
#   sbatch --array=5-12 jobs/01_basecalling.sh
#
# Submit from the repo root (config/ is read from there).
#
# Output: plain (uncompressed) FASTQ at basecalled/<run>.fastq on b1042.
# wf-single-cell reads .fastq and .fastq.gz equally, so there is no need to
# spend GPU-job time compressing.
#
# NOTE: the logs/ directory must exist BEFORE sbatch, or SLURM can't write the
# log and the job dies silently:
#   mkdir -p /projects/b1169/boles/tdp43_neuron_fans/logs

set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
SAMPLE_SHEET="${REPO_DIR}/config/samples.tsv"
POD5_ROOT="/projects/b1042/Gate_Lab/boles/tdp43_neuron_fans/pod5"
BASECALL_DIR="/projects/b1042/Gate_Lab/boles/tdp43_neuron_fans/basecalled"
DORADO="/projects/b1169/jackie/software/dorado-0.9.0-linux-x64/bin/dorado"
DORADO_MODEL="/projects/b1169/benney/ont_ref/dorado_model/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"

if [[ ! -f "${SAMPLE_SHEET}" ]]; then
  echo "ERROR: ${SAMPLE_SHEET} not found; submit this job from the repo root" >&2
  exit 1
fi

# +1 skips the header row of the sample sheet
RUN=$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR == n { print $1 }' "${SAMPLE_SHEET}")
if [[ -z "${RUN}" ]]; then
  echo "ERROR: no row ${SLURM_ARRAY_TASK_ID} in ${SAMPLE_SHEET}" >&2
  exit 1
fi

OUT_FASTQ="${BASECALL_DIR}/${RUN}.fastq"
TMP_FASTQ="${OUT_FASTQ}.partial"
if [[ -s "${OUT_FASTQ}" ]]; then
  echo "${RUN} already basecalled: ${OUT_FASTQ} (delete it to redo)"
  exit 0
fi

# Each run folder holds one MinKNOW run folder, which holds pod5/.
mapfile -t POD5_DIRS < <(find "${POD5_ROOT}/${RUN}" -mindepth 2 -maxdepth 2 -type d -name pod5 2>/dev/null)
if [[ ${#POD5_DIRS[@]} -ne 1 ]]; then
  echo "ERROR: expected exactly 1 pod5/ folder under ${POD5_ROOT}/${RUN}, found ${#POD5_DIRS[@]}" >&2
  exit 1
fi
POD5_DIR="${POD5_DIRS[0]}"

N_POD5=$(find "${POD5_DIR}" -maxdepth 1 -name '*.pod5' | wc -l)
if [[ ${N_POD5} -eq 0 ]]; then
  echo "ERROR: no .pod5 files in ${POD5_DIR} (transfer incomplete?)" >&2
  exit 1
fi

mkdir -p "${BASECALL_DIR}"

echo "Run              = ${RUN}"
echo "Dorado           = ${DORADO}"
${DORADO} --version 2>&1 || true
echo "Dorado model     = ${DORADO_MODEL}"
echo "POD5 folder      = ${POD5_DIR} (${N_POD5} files)"
echo "Output           = ${OUT_FASTQ}"
date

# --no-trim: wf-single-cell locates the 10x adapter / barcode / UMI / polyT
# itself and trims them after extraction. Leave the reads untouched here so
# dorado's adapter/primer trimming can't clip sequence that barcode
# identification relies on.
${DORADO} basecaller "${DORADO_MODEL}" "${POD5_DIR}" \
  --min-qscore 10 \
  --no-trim \
  --emit-fastq \
  > "${TMP_FASTQ}"

# Only expose the final file once dorado exits cleanly, so a failed or
# timed-out task never leaves a truncated FASTQ that 02 would pick up.
mv "${TMP_FASTQ}" "${OUT_FASTQ}"
date
echo "Done: ${OUT_FASTQ}"
