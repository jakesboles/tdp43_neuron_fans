#!/bin/bash
#SBATCH --job-name=02_doubletfinder
#SBATCH --array=1-4
#SBATCH --output=/projects/b1169/boles/tdp43_neuron_fans/logs/%x_%A_%a.log
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --account=p31535
#SBATCH --partition=normal

# One task per sample, in sorted orig.ident order from
# data/sc_data/01_qc/metadata.rds. The --array range must match the number
# of samples across all pools (4 per pool: 1-4 for pool 149 alone, 1-12 once
# all three pools are in). The R script stops with a clear error if a task
# id is out of range.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/tdp43_neuron_fans/scripts/sc_scripts/02_doubletfinder.R
