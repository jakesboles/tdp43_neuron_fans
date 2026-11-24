#!/bin/bash
#SBATCH -A p31535
#SBATCH -p gengpu
#SBATCH --gres=gpu:h100:1
#SBATCH -t 48:00:00
#SBATCH --mem=72G
#SBATCH --array=1-2
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --job-name=basecalling
#SBATCH -o /projects/b1042/Gate_Lab/boles/sc_lr/logs/%x_%A_%a.log

cd /projects/b1042/Gate_Lab/boles/sc_lr

SAMPLE=$(ls pod5s | sed -n "${SLURM_ARRAY_TASK_ID}p")
POD5_DIR=$(ls "pod5s/${SAMPLE}")
OUTPUT_DIR="fastqs"
DORADO_VERSION="/projects/b1169/jackie/software/dorado-0.9.0-linux-x64/bin/dorado"
DORADO_MODEL="/projects/b1169/benney/ont_ref/dorado_model/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"

echo "Dorado version = ${DORADO_VERSION}"
echo "Dorado model = ${DORADO_MODEL}"
echo "POD5 directory = ${POD5_DIR}"
echo "Output directory = ${OUTPUT_DIR}"

${DORADO_VERSION} basecaller "${DORADO_MODEL}" "pod5s/${SAMPLE}/${POD5_DIR}/pod5" --min-qscore 10 --emit-fastq > "${OUTPUT_DIR}/${SAMPLE}_calls.fastq"