#!/bin/bash

#SBATCH --job-name=02_doubletfinder
#SBATCH --array=1-4
#SBATCH --output=/projects/p31535/boles/JSB149/logs/%x_%A_%a.log
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --account=p31535
#SBATCH --partition=normal

cd /projects/p31535/boles/JSB149

module load R/4.4.0

PARAMS_FILE="02_doubletfinder_params.txt"

# this is for the one that failed for no reason
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f1 -d,)
CELL_COUNT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f2 -d,)

echo "Running DoubletFinder in JSB${SAMPLE} (${CELL_COUNT} cells)"

Rscript r_scripts/02_doubletfinder.R "${SAMPLE}" "${CELL_COUNT}"