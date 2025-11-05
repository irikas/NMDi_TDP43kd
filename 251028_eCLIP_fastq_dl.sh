#!/bin/bash

#SBATCH --job-name=eclip_tdp43  
#SBATCH --time=72:00:00
#SBATCH --partition=parallel
#SBATCH --ntasks=6
#SBATCH --mem-per-cpu=100GB

# ---------------------------
# Load modules
# ---------------------------
ml load sra-tools/3.0.3 Java/11.0.2 singularity/3.8.7 fastqc/0.11.9 cutadapt/3.2 

# ---------------------------
# Define variables
# ---------------------------
BASE_DIR="/data/jling2/irika/clipdata/251028_Tollervey_TDPCLIP_E-MTAB-530"
BIN_DIR="/data/jling2/irika/bin"
GENOME_DIR="${BIN_DIR}/genomes"
OUT_DIR="${BASE_DIR}/output"
INPUT_CSV="${BASE_DIR}/input/251028_eclip_info.csv"
DATE=date +"%Y%m%d_%H.%M.%S"

# ---------------------------
# Download and zip fastq 
# ---------------------------
cd $BASE_DIR/data
if 
fasterq-dump ERR039851 ERR039852 ERR039855 ERR039854 ERR039853 ERR039850
gzip *

# ---------------------------
# fastqc for all
# ---------------------------
mkdir $OUT_DIR/$DATE/fastqc
fastqc *.fastq.gz -o $OUT_DIR/$DATE/fastqc/

# ---------------------------
# Read trimming using CutAdapt
# Use universal Illumina Adapter as IDed by FastQC
# ---------------------------
mkdir $OUT_DIR/$DATE/CutAdapt

for file in *.fastq.gz; do
  cutadapt -j 0 -m 12 -a AGATCGGAAGAG -A AGATCGGAAGAG \
  --info-file info_${file}.txt \
  -o trimmed_${file} $file \
  > ${file}_cutadapt.log
done


# cd /data/jling2/irika/clipdata/251028_Tollervey_TDPCLIP_E-MTAB-530/data/
# fasterq-dump ERR039850
# fasterq-dump ERR039851

# gzip ERR039850.fastq
# gzip ERR039851.fastq

NXF_VER=22.10.6 /data/jling2/irika/bin/nextflow run nf-core/clipseq \
	--input /data/jling2/irika/clipdata/251028_Tollervey_TDPCLIP_E-MTAB-530/input/251028_eclip_design.csv \
	--outdir /data/jling2/irika/clipdata/251028_Tollervey_TDPCLIP_E-MTAB-530/output/ \
	--fasta /data/jling2/irika/bin/genomes/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
	--gtf /data/jling2/irika/bin/genomes/Homo_sapiens.GRCh38.114.gtf.gz \
	--email "isinha1@jh.edu" -profile singularity

module purge

exit 0

#ending the run
#python3 /data/jling2/bin/email/sendemail.py -s "eclip_tdp43" -c "eclip_tdp43 fastq finished running" -r isinha1@jh.edu
#python3 /data/jling2/irika/bin/Scripts/KillSession.py
