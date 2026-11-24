#!/bin/bash
#SBATCH -A b1169
#SBATCH -p b1169
#SBATCH -t 24:00:00
#SBATCH --mem=1G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=nf-AN73-2
#SBATCH -o /projects/b1042/Gate_Lab/boles/sc_lr/logs/%j-nf-73-2.out
#SBATCH --mail-user=jake.boles@northwestern.edu
#SBATCH --mail-type=END,FAIL


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
# - to run: just submit this script (only need 1G of mem, but set it to the max b1169 time possible)
# - if a step fails: fix whatever caused it, then just re-submit the same script with
#   the same output dir (do not change inputs/outputs), then -resume picks up where it left off
# - heads up: your samples are huuuuge, may take a couple tries to get the right resource amounts per 
#   process, but that's chill -- adjust in the config, re-submit, let -resume handle the rest.
# ==================================================================================================
#


module purge all
module load nextflow/24.04.4
# ==================================================================================================

# EDIT PER SAMPLE -- just change this and then sbatch /projects/p32095/nf/10x_LR_preprocessing/nf.sh
# ==================================================================================================
export SAMPLE_NO="73-2"
export IS_GEM_X=true
POD5_DIRS=(
  "/projects/b1042/Gate_Lab/boles/sc_lr/pod5s/73-2/20260703_1950_1C_PBI52126_73f0fca7/pod5"
)
export OUTPUT_DIR="/projects/b1042/Gate_Lab/boles/sc_lr/nf_outs/${SAMPLE_NO}"
export FORCED_CELL_NUM=18000 # change this to the appropriate cell number for this sample

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

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

CONFIG_FILE_PATH="/projects/p32095/nf/10x_LR_preprocessing/nextflow.config"
nextflow -log "${OUTPUT_DIR}/nextflow.log" run /projects/p32095/nf/10x_LR_preprocessing/main.nf \
  -c "$CONFIG_FILE_PATH" \
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