#!/bin/bash

#SBATCH --job-name=eclip_tdp43  
#SBATCH --time=72:00:00
#SBATCH --partition=parallel
#SBATCH --ntasks=6
#SBATCH --mem-per-cpu=100GB

# ---------------------------
# Load modules
# ---------------------------
ml load sra-tools/3.0.3 singularity/3.8.7 fastqc/0.11.9 STAR/2.7.10a samtools/1.10

# ---------------------------
# Define variables
# ---------------------------
BASE_DIR="/data/jling2/irika/clipdata/251028_Tollervey_TDPCLIP_E-MTAB-530"
BIN_DIR="/data/jling2/irika/bin"
GENOME_DIR="${BIN_DIR}/genomes/STAR_hg38_indices"
OUT_DIR="${BASE_DIR}/output"
SING_DIR="${BASE_DIR}/scripts/singularity_images"
INPUT_CSV="${BASE_DIR}/input/251028_eclip_info.csv"
DATE=date +"%Y%m%d_%H.%M.%S"
OUT_DIR="${OUT_DIR}/$DATE"

mkdir $OUT_DIR

# ---------------------------
# Singularity images
# ---------------------------
cd $SING_DIR
singularity pull --name cutadapt_5.2.sif docker://quay.io/biocontainers/cutadapt:5.2--py310h1fe012e_0

# ---------------------------
# Download and zip fastq 
# ---------------------------
cd $BASE_DIR/data
fasterq-dump ERR039851 ERR039852 ERR039855 ERR039854 ERR039853 ERR039850
gzip *

# ---------------------------
# fastqc for all
# ---------------------------
mkdir $OUT_DIR/fastqc
fastqc *.fastq.gz -o $OUT_DIR/fastqc/

python3 /data/jling2/bin/email/sendemail.py -s "clip_tdp43" -c "clip_tdp43 fastqc complete" -r isinha1@jh.edu

# ---------------------------
# Read trimming using CutAdapt
# Use universal Illumina Adapter as IDed by FastQC
# ---------------------------
cd $BASE_DIR/data
mkdir $OUT_DIR/CutAdapt

for file in *.fastq.gz; do
  singularity exec $SING_DIR/cutadapt_5.2.sif -j 0 -m 12 -a AGATCGGAAGAG -A AGATCGGAAGAG \\
  --info-file info_${file}.txt \\
  -o $OUT_DIR/CutAdapt/trimmed_${file} $file \\
  --json=$OUT_DIR/CutAdapt/${file}.cutadapt.json
done

python3 /data/jling2/bin/email/sendemail.py -s "clip_tdp43" -c "clip_tdp43 CutAdapt complete" -r isinha1@jh.edu

# ---------------------------
# STAR for mapping to genome
# ---------------------------

# Generate index
STAR --runMode genomeGenerate --runThreadN 6 --genomeDir $GENOME_DIR/indices --genomeFastaFiles $GENOME_DIR/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile $GENOME_DIR/Homo_sapiens.GRCh38.114.gtf --sjdbOverhang 100

python3 /data/jling2/bin/email/sendemail.py -s "clip_tdp43" -c "index finished running" -r isinha1@jh.edu

# Align reads
mkdir $OUT_DIR/STARoutput

for file in $OUT_DIR/CutAdapt/*; do \\
  STAR --runMode alignReads --runThreadN 6 --genomeDir $GENOME_DIR/STAR_hg38_indices \\
    --readFilesIn $file --readFilesCommand gunzip --outFileNamePrefix $OUT_DIR/STARoutput/${file%.fastq.gz} \\
    --outSAMtype BAM Unsorted
done

# Sort using samtools
cd $OUT_DIR/STARoutput

for file in *.Aligned.out.bam; do \\
  samtools sort -@ 6 -o ${file%.Aligned.out.bam}.Aligned.sortedByCoord.out.bam $file
done

for file in *.Aligned.sortedByCoord.out.bam; do \\
  samtools index -@ 6 -o $file
done
    
## Basic tools complete

# ---------------------------
# Ending the run
# ---------------------------
python3 /data/jling2/bin/email/sendemail.py -s "clip_tdp43" -c "clip_tdp43 fastq basic tools finished running" -r isinha1@jh.edu
python3 /data/jling2/irika/bin/Scripts/KillSession.py
