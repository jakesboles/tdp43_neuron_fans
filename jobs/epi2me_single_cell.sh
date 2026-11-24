#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH --mem=225G
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --job-name=JSB149_full_epi2me
#SBATCH -o /projects/b1042/Gate_Lab/boles/sc_lr/logs/%x-%j.out

cd /projects/b1042/Gate_Lab/boles/sc_lr/fastqs/JSB149
# annoying - need to submit the job from the same place every time so 
# nextflow can find its own run history to resume a job
# in the future make sure a starting point is set deliberately

module load cuda/12.0.1-gcc-12.3.0
module load nextflow/25.10.3
module load singularity

# ignore
export SINGULARITYENV_CUDA_VISIBLE_DEVICES=0

export NXF_HOME=/home/egn9175/.nextflow/assets/epi2me-labs/wf-single-cell
# export NXF_SINGULARITY_CACHEDIR=/gpfs/projects/p31535/boles/epi2me/wf-single-cell/singularity_cache
# export SINGULARITY_TMPDIR=/gpfs/projects/p31535/boles/epi2me/singularity_tmp
# export SINGULARITY_CACHEDIR=/gpfs/projects/p31535/boles/epi2me/singularity_cache

# export NXF_HOME=/projects/p32095/software/wf-single-cell/nxf_assets
export NXF_SINGULARITY_CACHEDIR=/projects/p32095/software/wf-single-cell/singularity_cache
export SINGULARITY_TMPDIR=/projects/p32095/software/.scratch/$USER/singularity_tmp
export SINGULARITY_CACHEDIR=/projects/p32095/software/.scratch/$USER/singularity_cache

export 'PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:512'
mkdir -p "$SINGULARITY_TMPDIR" "$SINGULARITY_CACHEDIR"

# update base dir, and have your fastq (.gz accepted) in $BASE/input folder
BASE="/projects/b1042/Gate_Lab/boles/sc_lr"
OUTPUT_DIR="$BASE/epi2me_outs/JSB149_merged"
WORK_DIR="$BASE/work"
STORE_DIR="$BASE/store"
FASTQ_DIR="/projects/b1042/Gate_Lab/boles/sc_lr/fastqs/JSB149"

# params
KIT="3prime:v4" # GEM-X
CELL_NUM=18512
SAMPLE_NAME="JSB149"

mkdir -p "$OUTPUT_DIR" "$BASE/logs" "$WORK_DIR" "$STORE_DIR"


# can use your own ref genome here, but I just used this prebuilt 10x reference genome bundle for convenience
#          to use your own, replace --epi2me_resource_bundle with --ref_genome_dir. expected format:
#          └── ref_genome_dir
#                ├── fasta
#                │   └── genome.fa (or genome.fa.gz)
#                └── genes
#                    └── genes.gtf (or genes.gtf.gz)


nextflow -log "$BASE/logs/nextflow.log" run epi2me-labs/wf-single-cell \
  -revision v3.3.4 \
  -profile singularity \
  -work-dir ${WORK_DIR} \
  --fastq ${FASTQ_DIR} \
  --kit ${KIT} \
  --expected_cells ${CELL_NUM} \
  --ref_genome_dir /projects/p31535/boles/cellranger_references/refdata-gex-GRCh38-2024-A \
  --store_dir ${STORE_DIR} \
  --sample ${SAMPLE_NAME} \
  --threads 32 \
  --out_dir ${OUTPUT_DIR} \
  -resume