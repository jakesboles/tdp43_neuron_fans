#!/bin/bash
#SBATCH -A p32095
#SBATCH -p gengpu
#SBATCH --gres=gpu:h100:1
#SBATCH -t 48:00:00
#SBATCH --mem=72G
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --job-name=01_basecalling
#SBATCH -o /projects/b1042/Gate_Lab/boles/tdp43_neuron_fans/logs/%x_%A_%a.log

# Basecall ONE chunk of POD5 files per array task.
# Chunks come from the manifest written by jobs/00_make_basecall_chunks.sh,
# which also prints the --array range to submit with, e.g.
#   sbatch --array=1-42 jobs/01_basecalling.sh
#
# Output: plain (uncompressed) FASTQ at
#   basecalled/<run>/<run>_partNNN.fastq
# wf-single-cell reads .fastq and .fastq.gz equally; compressing here would
# only add CPU time inside the 48 h GPU window.
#
# NOTE: the logs/ directory must exist BEFORE sbatch, or SLURM can't write the
# log and the job dies silently:
#   mkdir -p /projects/b1042/Gate_Lab/boles/tdp43_neuron_fans/logs

set -euo pipefail

BASECALL_DIR="/projects/b1042/Gate_Lab/boles/tdp43_neuron_fans/basecalled"
MANIFEST="${BASECALL_DIR}/chunks.tsv"
DORADO="/projects/b1169/jackie/software/dorado-0.9.0-linux-x64/bin/dorado"
DORADO_MODEL="/projects/b1169/benney/ont_ref/dorado_model/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: submit as an array, e.g. sbatch --array=1-N jobs/01_basecalling.sh" >&2
  echo "       (jobs/00_make_basecall_chunks.sh prints the range)" >&2
  exit 1
fi

# manifest columns: chunk_id run n_pod5 pod5_list fastq
read -r CHUNK RUN N_POD5 POD5_LIST OUT_FASTQ < <(
  awk -F'\t' -v id="${SLURM_ARRAY_TASK_ID}" 'NR > 1 && $1 == id { print $1, $2, $3, $4, $5 }' "${MANIFEST}"
) || true
if [[ -z "${CHUNK:-}" ]]; then
  echo "ERROR: chunk ${SLURM_ARRAY_TASK_ID} not in ${MANIFEST}" >&2
  exit 1
fi

if [[ -s "${OUT_FASTQ}" ]]; then
  echo "Chunk ${CHUNK} already done: ${OUT_FASTQ}"
  exit 0
fi

# dorado takes a directory, so give it one holding only this chunk's files
CHUNK_POD5_DIR="${POD5_LIST%.txt}_pod5"
rm -rf "${CHUNK_POD5_DIR}"
mkdir -p "${CHUNK_POD5_DIR}" "$(dirname "${OUT_FASTQ}")"
while read -r f; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: ${f} listed in ${POD5_LIST} no longer exists" >&2
    exit 1
  fi
  ln -s "${f}" "${CHUNK_POD5_DIR}/"
done < "${POD5_LIST}"

TMP_FASTQ="${OUT_FASTQ}.partial"

echo "Chunk            = ${CHUNK} (${RUN}, ${N_POD5} pod5 files)"
echo "Dorado           = ${DORADO}"
${DORADO} --version 2>&1 || true
echo "Dorado model     = ${DORADO_MODEL}"
echo "Output           = ${OUT_FASTQ}"
date

# --no-trim: wf-single-cell locates the 10x adapter / barcode / UMI / polyT
# itself and trims them after extraction. Leave the reads untouched here so
# dorado's adapter/primer trimming can't clip sequence that barcode
# identification relies on.
${DORADO} basecaller "${DORADO_MODEL}" "${CHUNK_POD5_DIR}" \
  --min-qscore 10 \
  --no-trim \
  --emit-fastq \
  > "${TMP_FASTQ}"

# Only expose the final file once dorado exits cleanly, so a failed or
# timed-out task never leaves a truncated FASTQ that 02 would pick up.
mv "${TMP_FASTQ}" "${OUT_FASTQ}"
rm -rf "${CHUNK_POD5_DIR}"
date
echo "Done: ${OUT_FASTQ}"
