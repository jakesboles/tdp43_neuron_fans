#!/bin/bash
# Split every sequencing run's POD5 files into chunks for 01_basecalling.sh.
#
# Basecalling a whole PromethION run on one GPU can exceed the 48 h gengpu
# limit, so each 01 array task basecalls one chunk of POD5 files instead.
# Run this on a login node (it only lists files), from the repo root:
#
#   bash jobs/00_make_basecall_chunks.sh [FILES_PER_CHUNK]   # default 10
#
# It writes a manifest (chunks.tsv) plus one file list per chunk under
# ${BASECALL_DIR}/chunks/, and prints the sbatch command for 01 with the
# right --array range.
#
# Re-run it whenever POD5s are added (e.g. once the 158-A/158-B transfers
# finish) -- but NOT while 01 or 02 jobs are running, since they read the
# manifest. Chunks whose FASTQ already exists are listed as done and skipped
# in the suggested --array range.

set -euo pipefail

FILES_PER_CHUNK="${1:-10}"
REPO_DIR="$(pwd)"
SAMPLE_SHEET="${REPO_DIR}/config/samples.tsv"
POD5_ROOT="/projects/b1042/Gate_Lab/boles/tdp43_neuron_fans/pod5"
BASECALL_DIR="/projects/b1042/Gate_Lab/boles/tdp43_neuron_fans/basecalled"

if [[ ! -f "${SAMPLE_SHEET}" ]]; then
  echo "ERROR: ${SAMPLE_SHEET} not found; run this from the repo root" >&2
  exit 1
fi
if ! [[ "${FILES_PER_CHUNK}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: FILES_PER_CHUNK must be a positive integer, got '${FILES_PER_CHUNK}'" >&2
  exit 1
fi

CHUNK_DIR="${BASECALL_DIR}/chunks"
MANIFEST="${BASECALL_DIR}/chunks.tsv"
mkdir -p "${CHUNK_DIR}"
rm -f "${CHUNK_DIR}"/chunk_*.txt

printf 'chunk_id\trun\tn_pod5\tpod5_list\tfastq\n' > "${MANIFEST}"

chunk_id=0
todo=()
while IFS=$'\t' read -r run pool; do
  [[ "${run}" == "sample" ]] && continue   # header

  mapfile -t pod5_dirs < <(find "${POD5_ROOT}/${run}" -mindepth 2 -maxdepth 2 -type d -name pod5 2>/dev/null)
  if [[ ${#pod5_dirs[@]} -ne 1 ]]; then
    echo "SKIP ${run}: expected 1 pod5/ folder, found ${#pod5_dirs[@]}" >&2
    continue
  fi
  mapfile -t files < <(find "${pod5_dirs[0]}" -maxdepth 1 -name '*.pod5' | sort -V)
  if [[ ${#files[@]} -eq 0 ]]; then
    echo "SKIP ${run}: no .pod5 files in ${pod5_dirs[0]}" >&2
    continue
  fi

  n_run_chunks=0
  for ((i = 0; i < ${#files[@]}; i += FILES_PER_CHUNK)); do
    chunk_id=$((chunk_id + 1))
    n_run_chunks=$((n_run_chunks + 1))
    list="${CHUNK_DIR}/chunk_$(printf '%04d' "${chunk_id}").txt"
    printf '%s\n' "${files[@]:i:FILES_PER_CHUNK}" > "${list}"
    n=$(wc -l < "${list}")
    fastq="${BASECALL_DIR}/${run}/${run}_part$(printf '%03d' "${n_run_chunks}").fastq"
    printf '%s\t%s\t%s\t%s\t%s\n' "${chunk_id}" "${run}" "${n}" "${list}" "${fastq}" >> "${MANIFEST}"
    [[ -s "${fastq}" ]] || todo+=("${chunk_id}")
  done
  echo "${run} (pool ${pool}): ${#files[@]} pod5 files -> ${n_run_chunks} chunk(s)"
done < "${SAMPLE_SHEET}"

echo
echo "Manifest: ${MANIFEST} (${chunk_id} chunks, ${#todo[@]} not yet basecalled)"
if [[ ${#todo[@]} -gt 0 ]]; then
  # compress e.g. 1 2 3 7 8 -> 1-3,7-8 for --array
  range=$(printf '%s\n' "${todo[@]}" | awk 'NR==1{s=p=$1; next} $1==p+1{p=$1; next} {printf "%s,", (s==p?s:s"-"p); s=p=$1} END{print (s==p?s:s"-"p)}')
  echo "Submit with:"
  echo "  sbatch --array=${range} jobs/01_basecalling.sh"
fi
