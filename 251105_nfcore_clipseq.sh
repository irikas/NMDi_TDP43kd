#!/bin/bash

#SBATCH --job-name=eclip_tdp43
#SBATCH --time=72:00:00
#SBATCH --partition=parallel
#SBATCH --ntasks=30
#SBATCH --mem=150G

# ---------------------------
# Define variables
# ---------------------------
MAIN_DIR="/data/jling2/irika/clipdata/251028_Tollervey_TDPCLIP_E-MTAB-530"
GENOME_DIR="/data/jling2/irika/bin/genomes/GRCh38_v115"
OUT_DIR="${MAIN_DIR}/output/$(date +"%Y%m%d_%H%M%S")"
SING_DIR="$HOME/singularity"
export SINGULARITY_TMPDIR=$SING_DIR/tmp
export SINGULARITY_CACHEDIR=$SING_DIR/singularity_cache
mkdir -p $SING_DIR $SINGULARITY_TMPDIR $SINGULARITY_CACHEDIR $OUT_DIR ${OUT_DIR}/ERR03985{0..5}
# ---------------------------
# Load modules
# ---------------------------
ml load Java/11.0.2 singularity/3.8.7

# ---------------------------
# Run nf-core/clipseq
# GRCh37 bc this version of iCount is not compatible with GRCh38
# ---------------------------
cd $MAIN_DIR/input

for file in 251106_*; do \
  fName=${file##251106_eclip_design_}
  fName=${fName%%.csv}
  NXF_VER=22.10.6 /data/jling2/irika/bin/nextflow run nf-core/clipseq --input "${MAIN_DIR}/input/${file}" --outdir "${OUT_DIR}/ERR0398${fName}" \
    --genome "GRCh37" \
    --peakcaller "iCount" --motif true --max_cpus 28 --max_time 72.h \
    --email isinha1@jh.edu -profile singularity
done


