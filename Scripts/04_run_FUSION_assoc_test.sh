#!/bin/bash
# Run FUSION association tests chromosome-by-chromosome using a SLURM array.
#SBATCH --job-name=twas_fusion
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --array=1-22
#SBATCH --output=./SLURM_OUT/%x_%A_%a.out
#SBATCH --error=./SLURM_OUT/%x_%A_%a.err

# --- Environment setup ---
# Adjust module loads or conda activation based on your cluster environment
conda activate pheast

CHR=$SLURM_ARRAY_TASK_ID
GWAS="path/to/GWAS/file"
WEIGHT="path/to/weight.pos"
WEIGHT_DIR="/path/to/weight"
LD="prefix_to_LD_files_by_chr"
OUT="path/to/FUSION/output"
MODALITY="modality_name"

mkdir -p "${OUT}"

echo "Running FUSION association test for ${MODALITY} modality, Chromosome ${CHR}"

Rscript /path/to/FUSION.assoc_test.R \
  --sumstats ${GWAS} \
  --weights ${WEIGHT} \
  --weights_dir ${WEIGHT_DIR} \
  --ref_ld_chr ${LD} \
  --chr ${CHR} \
  --out ${OUT}/${MODALITY}_chr${CHR}