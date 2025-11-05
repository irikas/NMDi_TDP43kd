#!/bin/bash

#SBATCH --job-name=eclip_tdp43
#SBATCH --time=72:00:00
#SBATCH --partition=parallel
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=100GB

# ---------------------------
# Define variables
# ---------------------------
MAIN_DIR="/data/jling2/irika/clipdata/251028_Tollervey_TDPCLIP_E-MTAB-530"
GENOME_DIR="/data/jling2/irika/bin/genomes/GRCh38_v115"
OUT_DIR="${MAIN_DIR}/output/$(date +"%Y%m%d_%H%M%S")"
DESIGN_FILE="${MAIN_DIR}/input/251105_eclip_design.csv"
SING_DIR="$HOME/singularity"
export SINGULARITY_TMPDIR=$SING_DIR/tmp
export SINGULARITY_CACHEDIR=$SING_DIR/singularity_cache
mkdir -p $SING_DIR $SINGULARITY_TMPDIR $SINGULARITY_CACHEDIR $OUT_DIR

# ---------------------------
# Load modules
# ---------------------------
ml load singularity/3.8.7

# ---------------------------
# Run nf-core/clipseq
# ---------------------------
NXF_VER=22.10.6 /data/jling2/irika/bin/nextflow run nf-core/clipseq --input "${DESIGN_FILE}" --outdir "${OUT_DIR}" \
  --fasta "${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
  --fai "${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai" \
  --gtf "${GENOME_DIR}/Homo_sapiens.GRCh38.115.gtf " \
  --star_index "${GENOME_DIR}/indices" \
  --peakcaller "iCount" --motif true \
  --email isinha1@jh.edu -profile singularity
