#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH --mem=225G
#SBATCH --array=1-3
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --job-name=02_epi2me
#SBATCH -o /projects/b1169/boles/tdp43_neuron_fans/logs/%x_%A_%a.log

# Run epi2me-labs/wf-single-cell once per 10x pool (one array task per row of
# config/pools.tsv). The basecalled BAMs of every sequencing run in the pool
# (config/samples.tsv) are fed in together as ONE sample, because they share a
# single cell barcode space.
#
# Submit from the repo root (config/ is read from there), after all 01 tasks
# have finished, e.g.
#   sbatch --dependency=afterok:<01 job id> jobs/02_epi2me_single_cell.sh
# or for one pool only:  sbatch --array=2 jobs/02_epi2me_single_cell.sh
#
# Storage:
#   b1169 (kept):    epi2me_outs/<pool>   results
#                    epi2me_runs/<pool>   nextflow launch dir, history, log
#   b1042 (scratch): epi2me_work/<pool>   nextflow work dir (large, disposable)
#                    epi2me_store/<pool>  cached reference index
# To resume a failed pool just resubmit it; -resume works as long as the
# b1042 work dir is still there.
#
# NOTE: /projects/b1169/boles/tdp43_neuron_fans/logs must exist before sbatch.

module load nextflow/25.10.3
module load singularity

# after module loads: Lmod's own scripts can trip over `set -u`
set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
SAMPLE_SHEET="${REPO_DIR}/config/samples.tsv"
POOL_SHEET="${REPO_DIR}/config/pools.tsv"
SCRATCH_DIR="/projects/b1042/Gate_Lab/boles/tdp43_neuron_fans"
ANALYSIS_DIR="/projects/b1169/boles/tdp43_neuron_fans"
BASECALL_DIR="${SCRATCH_DIR}/basecalled"
REF_GENOME_DIR="/projects/p31535/boles/cellranger_references/refdata-gex-GRCh38-2024-A"
WF_REVISION="v3.3.4"
KIT="3prime:v4" # GEM-X 3' v4

export NXF_SINGULARITY_CACHEDIR=/projects/p32095/software/wf-single-cell/singularity_cache
export SINGULARITY_TMPDIR=/projects/p32095/software/.scratch/$USER/singularity_tmp
export SINGULARITY_CACHEDIR=/projects/p32095/software/.scratch/$USER/singularity_cache
mkdir -p "$SINGULARITY_TMPDIR" "$SINGULARITY_CACHEDIR"

for f in "${SAMPLE_SHEET}" "${POOL_SHEET}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: ${f} not found; submit this job from the repo root" >&2
    exit 1
  fi
done

# ---- which pool is this task? ----------------------------------------------
# +1 skips the header row
read -r POOL CELL_NUM < <(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR == n { print $1, $2 }' "${POOL_SHEET}") || true
if [[ -z "${POOL:-}" ]]; then
  echo "ERROR: no row ${SLURM_ARRAY_TASK_ID} in ${POOL_SHEET}" >&2
  exit 1
fi
if ! [[ "${CELL_NUM}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: expected_cells for pool ${POOL} is '${CELL_NUM}'; set it in ${POOL_SHEET}" >&2
  exit 1
fi
SAMPLE_NAME="JSB${POOL}"

RUN_DIR="${ANALYSIS_DIR}/epi2me_runs/${POOL}"     # nextflow launch dir (.nextflow/ history lives here)
POOL_BAM_DIR="${RUN_DIR}/bam"                     # symlinks to this pool's unaligned BAMs only
OUTPUT_DIR="${ANALYSIS_DIR}/epi2me_outs/${POOL}"
WORK_DIR="${SCRATCH_DIR}/epi2me_work/${POOL}"
STORE_DIR="${SCRATCH_DIR}/epi2me_store/${POOL}"   # per pool: concurrent tasks must not write the same storeDir

mkdir -p "${RUN_DIR}" "${WORK_DIR}" "${STORE_DIR}" "${OUTPUT_DIR}"

# ---- collect this pool's BAMs ----------------------------------------------
# wf-single-cell treats a flat directory of (u)BAMs as a single sample.
# (A directory of sub-directories would instead be read as ONT-barcoded
# samples, which is not what we want here.)
rm -rf "${POOL_BAM_DIR}"
mkdir -p "${POOL_BAM_DIR}"
mapfile -t RUNS < <(awk -F'\t' -v p="${POOL}" 'NR > 1 && $2 == p { print $1 }' "${SAMPLE_SHEET}")
if [[ ${#RUNS[@]} -eq 0 ]]; then
  echo "ERROR: no sequencing runs listed for pool ${POOL} in ${SAMPLE_SHEET}" >&2
  exit 1
fi
for RUN in "${RUNS[@]}"; do
  BAM="${BASECALL_DIR}/${RUN}.bam"
  if [[ ! -s "${BAM}" ]]; then
    # Fail rather than silently running the pool without one of its runs.
    echo "ERROR: missing basecalled BAM for ${RUN} (pool ${POOL}): ${BAM}" >&2
    exit 1
  fi
  ln -s "${BAM}" "${POOL_BAM_DIR}/"
done

# ---- cap nextflow's local executor at what SLURM gave us ----------------
# With -profile singularity every process runs inside this one allocation.
# By default the local executor assumes it owns the whole node's memory.
MEM_GB=$(( ${SLURM_MEM_PER_NODE:-230400} / 1024 * 9 / 10 ))   # leave ~10% for the nextflow JVM
cat > "${RUN_DIR}/local_executor.config" <<EOF
executor {
  cpus = ${SLURM_CPUS_ON_NODE:-32}
  memory = '${MEM_GB} GB'
}
EOF

echo "Pool             = ${POOL}"
echo "Sample name      = ${SAMPLE_NAME}"
echo "Expected cells   = ${CELL_NUM}"
echo "Runs             = ${RUNS[*]}"
echo "Executor limits  = ${SLURM_CPUS_ON_NODE:-32} cpus, ${MEM_GB} GB"
echo "Work dir         = ${WORK_DIR}"
echo "Output           = ${OUTPUT_DIR}"

cd "${RUN_DIR}"

nextflow -log "${RUN_DIR}/nextflow.log" run epi2me-labs/wf-single-cell \
  -revision "${WF_REVISION}" \
  -profile singularity \
  -c "${RUN_DIR}/local_executor.config" \
  -work-dir "${WORK_DIR}" \
  --bam "${POOL_BAM_DIR}" \
  --sample "${SAMPLE_NAME}" \
  --kit "${KIT}" \
  --expected_cells "${CELL_NUM}" \
  --ref_genome_dir "${REF_GENOME_DIR}" \
  --store_dir "${STORE_DIR}" \
  --threads "${SLURM_CPUS_ON_NODE:-32}" \
  --out_dir "${OUTPUT_DIR}" \
  -resume
