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
# Output: unaligned BAM (dorado's default, already compressed) at
# basecalled/<run>.bam on b1042. wf-single-cell takes it directly with --bam,
# and it keeps dorado's run/model metadata that FASTQ drops.
#
# Timeouts: while dorado runs, output goes to <run>.bam.partial. If a task is
# killed (e.g. at the 48 h limit), just resubmit the same array index: the
# partial file is moved to <run>.bam.resume and dorado --resume-from copies
# those reads over and basecalls only the rest.
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

OUT_BAM="${BASECALL_DIR}/${RUN}.bam"
PARTIAL_BAM="${OUT_BAM}.partial"   # written while dorado runs
RESUME_BAM="${OUT_BAM}.resume"     # reads kept from an interrupted run
if [[ -s "${OUT_BAM}" ]]; then
  echo "${RUN} already basecalled: ${OUT_BAM} (delete it to redo)"
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

# ---- resume an interrupted run, if there is one -------------------------
# A resumed run copies the old reads into its new .partial first, so after
# repeated timeouts keep whichever of .partial / .resume holds more data.
size() { if [[ -f "$1" ]]; then stat -c %s "$1"; else echo 0; fi; }
if [[ $(size "${PARTIAL_BAM}") -gt $(size "${RESUME_BAM}") ]]; then
  mv -f "${PARTIAL_BAM}" "${RESUME_BAM}"
fi
rm -f "${PARTIAL_BAM}"
RESUME_ARGS=()
if [[ -s "${RESUME_BAM}" ]]; then
  echo "Resuming from ${RESUME_BAM} ($(size "${RESUME_BAM}") bytes)"
  RESUME_ARGS=(--resume-from "${RESUME_BAM}")
fi

echo "Run              = ${RUN}"
echo "Dorado           = ${DORADO}"
${DORADO} --version 2>&1 || true
echo "Dorado model     = ${DORADO_MODEL}"
echo "POD5 folder      = ${POD5_DIR} (${N_POD5} files)"
echo "Output           = ${OUT_BAM}"
date

# --no-trim: wf-single-cell locates the 10x adapter / barcode / UMI / polyT
# itself and trims them after extraction. Leave the reads untouched here so
# dorado's adapter/primer trimming can't clip sequence that barcode
# identification relies on.
${DORADO} basecaller "${DORADO_MODEL}" "${POD5_DIR}" \
  --min-qscore 10 \
  --no-trim \
  ${RESUME_ARGS[@]+"${RESUME_ARGS[@]}"} \
  > "${PARTIAL_BAM}"

# Check the BAM is complete (has its EOF block) if samtools is available
if command -v samtools >/dev/null 2>&1; then
  samtools quickcheck "${PARTIAL_BAM}"
fi

# Only expose the final file once dorado exits cleanly, so a failed or
# timed-out task never leaves a truncated BAM that 02 would pick up.
mv "${PARTIAL_BAM}" "${OUT_BAM}"
rm -f "${RESUME_BAM}"
date
echo "Done: ${OUT_BAM}"
