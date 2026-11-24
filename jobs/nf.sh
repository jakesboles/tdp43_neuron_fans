#!/bin/bash
#SBATCH -A b1169
#SBATCH -p b1169
#SBATCH -t 24:00:00
#SBATCH --mem=1G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=nf_10x_lr
#SBATCH -o /projects/b1169/boles/tdp43_neuron_fans/logs/%x-%j.out
#SBATCH --mail-user=jake.boles@northwestern.edu
#SBATCH --mail-type=END,FAIL

# Usage (from the repo root, one job per 10x pool):
#   sbatch jobs/nf.sh 149
#   sbatch jobs/nf.sh 158-A
#   sbatch jobs/nf.sh 158-B
# The pool's sequencing runs come from config/samples.tsv and the forced cell
# number from config/pools.tsv. This is the 10x_LR_preprocessing pipeline,
# an alternative to jobs/00-02 (it does its own basecalling); it is NOT needed
# to run 00-02.
#
# NOTE: /projects/b1169/boles/tdp43_neuron_fans/logs must exist before sbatch.

# note: no need to change the mem above, all resource allocation is handled in nextflow.config
# the 1G is used to just submit the jobs

# ==================================================================================================
# INTRO/USAGE
# ==================================================================================================
# - nextflow manages all resources per process via the nextflow.config
# - before submitting, open the config (in this folder) and skim it:
#        - check that you have access to/are happy w the accounts/paritions/time/memory allocs
#        - good practice to be super generous with resources, LR data gets unpredictable w its usage
# - processes in this pipeline: basecalling -> nanoplot -> scanner -> assigner -> curator
# - to run: sbatch jobs/nf.sh <pool> from the repo root (only need 1G of mem, but set it to the max b1169 time possible)
# - if a step fails: fix whatever caused it, then just re-submit the same script with
#   the same pool (do not change inputs/outputs), then -resume picks up where it left off
# - heads up: your samples are huuuuge, may take a couple tries to get the right resource amounts per 
#   process, but that's chill -- adjust in the config, re-submit, let -resume handle the rest.
# ==================================================================================================
#


module purge all
module load nextflow/24.04.4
# ==================================================================================================

# PER-POOL SETTINGS (read from config/, nothing to edit here)
# ==================================================================================================
REPO_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
SAMPLE_SHEET="${REPO_DIR}/config/samples.tsv"
POOL_SHEET="${REPO_DIR}/config/pools.tsv"
POD5_ROOT="/projects/b1042/Gate_Lab/boles/tdp43_neuron_fans/pod5"
ANALYSIS_DIR="/projects/b1169/boles/tdp43_neuron_fans"
SCRATCH_DIR="/projects/b1042/Gate_Lab/boles/tdp43_neuron_fans"

POOL="${1:-}"
if [[ -z "${POOL}" ]]; then
  echo "ERROR: give the pool as an argument, e.g. sbatch jobs/nf.sh 149" >&2
  exit 1
fi
for f in "${SAMPLE_SHEET}" "${POOL_SHEET}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: ${f} not found; submit this job from the repo root" >&2
    exit 1
  fi
done

FORCED_CELL_NUM=$(awk -F'\t' -v p="${POOL}" 'NR > 1 && $1 == p { print $2 }' "${POOL_SHEET}")
if ! [[ "${FORCED_CELL_NUM}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: expected_cells for pool '${POOL}' is '${FORCED_CELL_NUM}'; set it in ${POOL_SHEET}" >&2
  exit 1
fi
export FORCED_CELL_NUM

mapfile -t RUNS < <(awk -F'\t' -v p="${POOL}" 'NR > 1 && $2 == p { print $1 }' "${SAMPLE_SHEET}")
if [[ ${#RUNS[@]} -eq 0 ]]; then
  echo "ERROR: no sequencing runs for pool '${POOL}' in ${SAMPLE_SHEET}" >&2
  exit 1
fi

# one pod5/ folder per sequencing run; fail if any run is missing or empty
POD5_DIRS=()
for RUN in "${RUNS[@]}"; do
  mapfile -t found < <(find "${POD5_ROOT}/${RUN}" -mindepth 2 -maxdepth 2 -type d -name pod5 2>/dev/null)
  if [[ ${#found[@]} -ne 1 ]] || [[ -z "$(find "${found[0]}" -maxdepth 1 -name '*.pod5' -print -quit)" ]]; then
    echo "ERROR: run ${RUN} (pool ${POOL}) has no single, non-empty pod5/ folder under ${POD5_ROOT}/${RUN}" >&2
    exit 1
  fi
  POD5_DIRS+=("${found[0]}")
done

export SAMPLE_NO="JSB${POOL}"
export IS_GEM_X=true
export OUTPUT_DIR="${ANALYSIS_DIR}/nf_outs/${POOL}"
WORK_DIR="${SCRATCH_DIR}/nf_work/${POOL}"   # large intermediate files stay on b1042

# scanner params -- check this!
export SCANNING_REGION=100
export A5_SEQ="AAGCAGTGGTATCAACGCAGAGTACAT"
export A3_SEQ="CTACACGACGCTCTTCCGATCT"
export PT_SEQ="TTTTTTTTTTTT"
export UMI_LENGTH=12
export CB_LENGTH=16 # barcode length
# ==================================================================================================
# DONE EDITING
# ==================================================================================================


# REFERENCES
export REF_GENOME="/projects/p32095/refs/genomes/hg38.fa"
export IDX_GENOME="/projects/b1169/benney/ont_ref/ref_hg38.mmi"
POD5_DIRS_CSV=$(IFS=,; echo "${POD5_DIRS[*]}")

mkdir -p "$OUTPUT_DIR" "$WORK_DIR"

echo "Pool = ${POOL}; runs = ${RUNS[*]}"
printf '  %s\n' "${POD5_DIRS[@]}"
echo "Forced cell number = ${FORCED_CELL_NUM}"
echo "Output = ${OUTPUT_DIR}; work = ${WORK_DIR}"
cd "$OUTPUT_DIR"

CONFIG_FILE_PATH="/projects/p32095/nf/10x_LR_preprocessing/nextflow.config"
nextflow -log "${OUTPUT_DIR}/nextflow.log" run /projects/p32095/nf/10x_LR_preprocessing/main.nf \
  -c "$CONFIG_FILE_PATH" \
  -work-dir "$WORK_DIR" \
  --sample_no "$SAMPLE_NO" \
  --is_GEM_X "$IS_GEM_X" \
  --pod5_dirs "$POD5_DIRS_CSV" \
  --scanning_region "$SCANNING_REGION" \
  --a5_seq "$A5_SEQ" \
  --a3_seq "$A3_SEQ" \
  --pT_seq "$PT_SEQ" \
  --lUMI "$UMI_LENGTH" \
  --lCB "$CB_LENGTH" \
  --forced_cell_num "$FORCED_CELL_NUM" \
  --ref_genome "$REF_GENOME" \
  --idx_genome "$IDX_GENOME" \
  --base_output "$OUTPUT_DIR" \
  -profile quest_slurm \
  -resume