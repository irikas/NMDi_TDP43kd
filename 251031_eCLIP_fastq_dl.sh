#!/bin/bash

#SBATCH --job-name=eclip_tdp43  
#SBATCH --time=72:00:00
#SBATCH --partition=parallel
#SBATCH --ntasks=6
#SBATCH --mem-per-cpu=100GB

# ---------------------------
# Define variables
# ---------------------------
BASE_DIR="/data/jling2/irika/clipdata/251028_Tollervey_TDPCLIP_E-MTAB-530"
BIN_DIR="/data/jling2/irika/bin"
GENOME_DIR="${BIN_DIR}/genomes/STAR_hg38_indices"
DATE=$(date +"%Y%m%d_%H%M%S")
OUT_DIR="${BASE_DIR}/output/${DATE}"
SING_DIR="$HOME/singularity"
SING_IMAGES="$SING_DIR/singularity_images"
export SINGULARITY_TMPDIR=$SING_DIR/tmp
export SINGULARITY_CACHEDIR=$SING_DIR/singularity_cache
mkdir -p $SING_DIR $SING_IMAGES $SINGULARITY_TMPDIR $SINGULARITY_CACHEDIR $OUT_DIR

# ---------------------------
# Load modules
# ---------------------------
ml load sra-tools/3.0.3 singularity/3.8.7 fastqc/0.11.9 

python3 /data/jling2/bin/email/sendemail.py -s "clip_tdp43" -c "clip_tdp43 START" -r isinha1@jh.edu
# ---------------------------
# Singularity images
# ---------------------------
cd $SING_IMAGES
singularity pull --name trim-galore.sif docker://quay.io/biocontainers/trim-galore:0.6.9--hdfd78af_0

# ---------------------------
# Download and zip fastq 
# ---------------------------
cd $BASE_DIR/data
# fasterq-dump ERR039851 ERR039852 ERR039855 ERR039854 ERR039853 ERR039850
# gzip *

# ---------------------------
# fastqc for all
# ---------------------------
mkdir $OUT_DIR/fastqc
fastqc *.fastq.gz -o $OUT_DIR/fastqc/

python3 /data/jling2/bin/email/sendemail.py -s "clip_tdp43" -c "clip_tdp43 fastqc complete" -r isinha1@jh.edu

# ---------------------------
# Read trimming using Trim Galore!
# Use universal Illumina Adapter as IDed by FastQC
# ---------------------------

cd $BASE_DIR/data

for file in *.fastq.gz; do
  singularity exec $SING_IMAGES/trim-galore.sif trim_galore -o "${OUT_DIR}/Trimmed/" $file
done

python3 /data/jling2/bin/email/sendemail.py -s "clip_tdp43" -c "clip_tdp43 Trim Galore! complete" -r isinha1@jh.edu

# ---------------------------
# STAR for mapping to genome
# ---------------------------

# ---------------------------
# Load modules
# ---------------------------
ml load STAR/2.7.10a samtools/1.10 

# Generate index
# STAR --runMode genomeGenerate --runThreadN 6 --genomeDir $GENOME_DIR/indices --genomeFastaFiles $GENOME_DIR/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile $GENOME_DIR/Homo_sapiens.GRCh38.114.gtf --sjdbOverhang 100
# python3 /data/jling2/bin/email/sendemail.py -s "clip_tdp43" -c "index finished running" -r isinha1@jh.edu

# Align reads
mkdir $OUT_DIR/STARoutput

for file in ${OUT_DIR}/Trimmed/*trimmed.fq.gz;  \
  do STAR --runMode alignReads --runThreadN 6  \
  --genomeDir $GENOME_DIR/STAR_hg38_indices  \
  --readFilesIn $file --readFilesCommand gunzip  \
  --outFileNamePrefix $OUT_DIR/STARoutput/${file%.fastq.gz}  \
  --outSAMtype BAM Unsorted
done

# Sort using samtools
cd $OUT_DIR/STARoutput

for file in *.Aligned.out.bam; do \
  samtools sort -@ 6 -o ${file%.Aligned.out.bam}.Aligned.sortedByCoord.out.bam $file
done

for file in *.Aligned.sortedByCoord.out.bam; do \
  samtools index -@ 6 -o $file
done
    
## Basic tools complete

# ---------------------------
# Ending the run
# ---------------------------
python3 /data/jling2/bin/email/sendemail.py -s "clip_tdp43" -c "clip_tdp43 fastq basic tools finished running" -r isinha1@jh.edu
python3 /data/jling2/irika/bin/Scripts/KillSession.py
