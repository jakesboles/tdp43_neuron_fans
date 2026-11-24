#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH --mem=225G
#SBATCH --array=1-3
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --job-name=02_epi2me
#SBATCH -o /projects/b1042/Gate_Lab/boles/tdp43_neuron_fans/logs/%x_%A_%a.log

# Run epi2me-labs/wf-single-cell once per 10x pool (one array task per row of
# config/pools.tsv). All basecalled FASTQs from the sequencing runs that belong
# to a pool (config/samples.tsv) are fed in together as ONE sample, because
# they share a single cell barcode space.
#
# Run after every 01_basecalling task for the pool has finished, e.g.
#   sbatch --dependency=afterok:<01 job id> jobs/02_epi2me_single_cell.sh
#
# To resume a failed pool, just resubmit: each pool launches nextflow from its
# own directory (epi2me_runs/<pool>), so -resume finds the right run history
# no matter where sbatch was called from, and the 3 pools never share a lock.
#
# NOTE: logs/ must exist before sbatch (see 01_basecalling.sh).

module load nextflow/25.10.3
module load singularity

# after module loads: Lmod's own scripts can trip over `set -u`
set -euo pipefail

PROJECT_DIR="/projects/b1042/Gate_Lab/boles/tdp43_neuron_fans"
SAMPLE_SHEET="${PROJECT_DIR}/config/samples.tsv"
POOL_SHEET="${PROJECT_DIR}/config/pools.tsv"
FASTQ_DIR="${PROJECT_DIR}/basecalled"
REF_GENOME_DIR="/projects/p31535/boles/cellranger_references/refdata-gex-GRCh38-2024-A"
WF_REVISION="v3.3.4"
KIT="3prime:v4" # GEM-X 3' v4

export NXF_SINGULARITY_CACHEDIR=/projects/p32095/software/wf-single-cell/singularity_cache
export SINGULARITY_TMPDIR=/projects/p32095/software/.scratch/$USER/singularity_tmp
export SINGULARITY_CACHEDIR=/projects/p32095/software/.scratch/$USER/singularity_cache
mkdir -p "$SINGULARITY_TMPDIR" "$SINGULARITY_CACHEDIR"

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

RUN_DIR="${PROJECT_DIR}/epi2me_runs/${POOL}"      # nextflow launch dir (.nextflow/ history lives here)
POOL_FASTQ_DIR="${RUN_DIR}/fastq"                 # symlinks to this pool's FASTQs only
WORK_DIR="${RUN_DIR}/work"
STORE_DIR="${RUN_DIR}/store"                      # per pool: concurrent tasks must not write the same storeDir
OUTPUT_DIR="${PROJECT_DIR}/epi2me_outs/${POOL}"

mkdir -p "${RUN_DIR}" "${WORK_DIR}" "${STORE_DIR}" "${OUTPUT_DIR}"

# ---- collect this pool's FASTQs --------------------------------------------
# wf-single-cell treats a flat directory of FASTQs as a single sample.
# (A directory of sub-directories would instead be read as ONT-barcoded
# samples, which is not what we want here.)
rm -rf "${POOL_FASTQ_DIR}"
mkdir -p "${POOL_FASTQ_DIR}"
mapfile -t RUNS < <(awk -F'\t' -v p="${POOL}" 'NR > 1 && $2 == p { print $1 }' "${SAMPLE_SHEET}")
if [[ ${#RUNS[@]} -eq 0 ]]; then
  echo "ERROR: no sequencing runs listed for pool ${POOL} in ${SAMPLE_SHEET}" >&2
  exit 1
fi
for RUN in "${RUNS[@]}"; do
  FQ="${FASTQ_DIR}/${RUN}.fastq.gz"
  if [[ ! -s "${FQ}" ]]; then
    # Fail rather than silently running the pool without one of its runs.
    echo "ERROR: missing basecalled FASTQ for ${RUN} (pool ${POOL}): ${FQ}" >&2
    exit 1
  fi
  ln -s "${FQ}" "${POOL_FASTQ_DIR}/${RUN}.fastq.gz"
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
echo "Output           = ${OUTPUT_DIR}"

cd "${RUN_DIR}"

nextflow -log "${RUN_DIR}/nextflow.log" run epi2me-labs/wf-single-cell \
  -revision "${WF_REVISION}" \
  -profile singularity \
  -c "${RUN_DIR}/local_executor.config" \
  -work-dir "${WORK_DIR}" \
  --fastq "${POOL_FASTQ_DIR}" \
  --sample "${SAMPLE_NAME}" \
  --kit "${KIT}" \
  --expected_cells "${CELL_NUM}" \
  --ref_genome_dir "${REF_GENOME_DIR}" \
  --store_dir "${STORE_DIR}" \
  --threads "${SLURM_CPUS_ON_NODE:-32}" \
  --out_dir "${OUTPUT_DIR}" \
  -resume
