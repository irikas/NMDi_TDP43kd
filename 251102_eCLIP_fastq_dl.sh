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
GENOME_DIR="/data/jling2/irika/bin/genomes/GRCh38_v115"
OUT_DIR="${BASE_DIR}/output/$(date +"%Y%m%d_%H%M%S")"
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
singularity pull --name umi_tools.sif docker://quay.io/biocontainers/umi_tools:1.1.6--py310h1fe012e_0
singularity pull --name icount.sif docker://quay.io/biocontainers/icount:2.0.0--py_1

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
# umi_tools: Extract the UMIs
# ---------------------------
mkdir $OUT_DIR/UMIextracted
cd $BASE_DIR/data

singularity exec $SING_IMAGES/umi_tools.sif umi_tools extract --stdin=ERR039850.fastq.gz --bc-pattern=NNNCCCCNN --log=$OUT_DIR/UMIextracted/ERR039850_processed.log --stdout $OUT_DIR/UMIextracted/ERR039850_UMI.fastq.gz 
singularity exec $SING_IMAGES/umi_tools.sif umi_tools extract --stdin=ERR039851.fastq.gz --bc-pattern=NNCCNNCCC --log=$OUT_DIR/UMIextracted/ERR039851_processed.log --stdout $OUT_DIR/UMIextracted/ERR039851_UMI.fastq.gz 
singularity exec $SING_IMAGES/umi_tools.sif umi_tools extract --stdin=ERR039852.fastq.gz --bc-pattern=NNCCNNCCC --log=$OUT_DIR/UMIextracted/ERR039852_processed.log --stdout $OUT_DIR/UMIextracted/ERR039852_UMI.fastq.gz 
singularity exec $SING_IMAGES/umi_tools.sif umi_tools extract --stdin=ERR039853.fastq.gz --bc-pattern=NNNCCCCNN --log=$OUT_DIR/UMIextracted/ERR039853_processed.log --stdout $OUT_DIR/UMIextracted/ERR039853_UMI.fastq.gz 
singularity exec $SING_IMAGES/umi_tools.sif umi_tools extract --stdin=ERR039854.fastq.gz --bc-pattern=NNNCCCCNN --log=$OUT_DIR/UMIextracted/ERR039854_processed.log --stdout $OUT_DIR/UMIextracted/ERR039854_UMI.fastq.gz 
singularity exec $SING_IMAGES/umi_tools.sif umi_tools extract --stdin=ERR039855.fastq.gz --bc-pattern=CCNNN --log=$OUT_DIR/UMIextracted/ERR039855_processed.log --stdout $OUT_DIR/UMIextracted/ERR039855_UMI.fastq.gz 

# ---------------------------
# Read trimming using Trim Galore!
# Use universal Illumina Adapter as IDed by FastQC
# Make new folder named "Trimmed"
# ---------------------------

cd $OUT_DIR/UMIextracted

for file in *.fastq.gz; do
  singularity exec $SING_IMAGES/trim-galore.sif trim_galore -j 6 --fastqc --dont_gzip -o "${OUT_DIR}/Trimmed/" $file
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
mkdir indices
STAR --runMode genomeGenerate --runThreadN 6 --genomeDir $GENOME_DIR/indices --genomeFastaFiles $GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile $GENOME_DIR/Homo_sapiens.GRCh38.115.gtf --sjdbOverhang 100
python3 /data/jling2/bin/email/sendemail.py -s "clip_tdp43" -c "index finished running" -r isinha1@jh.edu

# Align reads
mkdir $OUT_DIR/STARoutput
cd $OUT_DIR/Trimmed

for file in *trimmed.fq; do \
  STAR --runMode alignReads --runThreadN 6  \
  --genomeDir $GENOME_DIR/indices  \
  --readFilesIn $file  \
  --outFileNamePrefix $OUT_DIR/STARoutput/${file%.fq}_  \
  --outSAMtype BAM Unsorted 
done

# Sort using samtools
cd $OUT_DIR/STARoutput

for file in *_Aligned.out.bam; do \
  samtools sort -@ 6 -o ${file%_Aligned.out.bam}_Aligned.sortedByCoord.out.bam $file
done

for file in *_Aligned.sortedByCoord.out.bam; do \
  samtools index -@ 6 $file
done

# ---------------------------
# UMI deduplication
# ---------------------------
mkdir $OUT_DIR/UMIdedup
cd $OUT_DIR/STARoutput
ml load singularity/3.8.7 

for file in *_Aligned.sortedByCoord.out.bam; do \
  singularity exec $SING_IMAGES/umi_tools.sif umi_tools dedup -I $file --output-stats=deduplicated -S $OUT_DIR/UMIdedup/${file%_Aligned.sortedByCoord.out.bam}_AlignSortDedup.bam
done

## Basic tools complete
python3 /data/jling2/bin/email/sendemail.py -s "clip_tdp43" -c "clip_tdp43 fastq basic tools finished running" -r isinha1@jh.edu

# ---------------------------
# HTseq: Prepare annotation
# HTseq: createSlidingWindows
# HTseq: mapToID
# ---------------------------
mkdir $GENOME_DIR/HTseq-CLIP
cd $GENOME_DIR/HTseq-CLIP

singularity exec $SING_IMAGES/htseq-clip.sif htseq-clip annotation -g $GENOME_DIR/Homo_sapiens.GRCh38.115.gff3 -o Homo_sapiens.GRCh38.115_flat.bed --unsorted
singularity exec $SING_IMAGES/htseq-clip.sif htseq-clip annotation -g $GENOME_DIR/gencode.v49.filtered.sorted.gtf  -o Homo_sapiens.GRCh38.115_flat.bed 

singularity exec $SING_IMAGES/htseq-clip.sif htseq-clip createSlidingWindows -i Homo_sapiens.GRCh38.115_flat.bed -o Homo_sapiens.GRCh38.115_windows.bed
singularity exec $SING_IMAGES/htseq-clip.sif htseq-clip mapToId -i Homo_sapiens.GRCh38.115_flat.bed -o Homo_sapiens.GRCh38.115_mapped.txt

# ---------------------------
# HTseq: Extract crosslink sites
# HTseq: Count crosslink sites
# ---------------------------
mkdir $OUT_DIR/HTseq-CLIP
cd $OUT_DIR/UMIdedup

for file in *_UMI_trimmed_AlignSortDedup.bam; do \
  singularity exec $SING_IMAGES/htseq-clip.sif htseq-clip extract -i $file -c 6 -o $OUT_DIR/HTseq-CLIP/${file%_UMI_trimmed_AlignSortDedup.bam}_XLink_Extraction.bed
  singularity exec $SING_IMAGES/htseq-clip.sif htseq-clip count -i $OUT_DIR/HTseq-CLIP/${file%_UMI_trimmed_AlignSortDedup.bam}_XLink_Extraction.bed -a $GENOME_DIR/HTseq-CLIP/Homo_sapiens.GRCh38.115_windows.bed -o $OUT_DIR/HTseq-CLIP/${file%_UMI_trimmed_AlignSortDedup.bam}_XLink_Counts.txt
done

# ---------------------------
# Ending the run
# ---------------------------
python3 /data/jling2/bin/email/sendemail.py -s "clip_tdp43" -c "clip_tdp43 HTSeq-clip finished running" -r isinha1@jh.edu
exit 0
