#!/bin/bash

#SBATCH --job-name=Rhine_clip_tdp43
#SBATCH --time=72:00:00
#SBATCH --partition=parallel
#SBATCH --ntasks=6
#SBATCH --mem-per-cpu=100GB

# ---------------------------
# Define variables
# ---------------------------
MAIN_DIR="/data/jling2/irika/clipdata/251231_Rhine_ipsc"
EXP_NAME="Rhine iPSCs"
GENOME_DIR="/data/jling2/irika/bin/genomes/GRCh38_v115"
OUT_DIR="${MAIN_DIR}/output/$(date +"%Y%m%d_%H%M%S")"
SING_DIR="$HOME/singularity"
SING_IMAGES="$SING_DIR/singularity_images"
export SINGULARITY_TMPDIR=$SING_DIR/tmp
export SINGULARITY_CACHEDIR=$SING_DIR/singularity_cache
mkdir -p $SING_DIR $SING_IMAGES $SINGULARITY_TMPDIR $SINGULARITY_CACHEDIR $MAIN_DIR/data \
  $OUT_DIR $OUT_DIR/fastqc $OUT_DIR/UMIextracted $OUT_DIR/STARoutput $OUT_DIR/UMIdedup $OUT_DIR/iCount

# ---------------------------
# Load modules
# ---------------------------
ml load sra-tools/3.0.3 singularity/3.8.7 fastqc/0.11.9
python3 /data/jling2/bin/email/sendemail.py -s "${EXP_NAME}: clip_tdp43" -c "clip_tdp43 START" -r isinha1@jh.edu

# ---------------------------
# Singularity images
# ---------------------------
cd $SING_IMAGES

if [ -e "trim-galore.sif" ]; then
  echo "Trim Galore! Singularity image already exists."
else
  singularity pull --name trim-galore.sif docker://quay.io/biocontainers/trim-galore:0.6.9--hdfd78af_0
fi

if [ -e "umi_tools.sif" ]; then
  echo "umi-tools Singularity image already exists."
else
  singularity pull --name umi_tools.sif docker://quay.io/biocontainers/umi_tools:1.1.6--py310h1fe012e_0
fi

if [ -e "icount.sif" ]; then
  echo "iCount Singularity image already exists."
else
  singularity pull --name icount.sif docker://quay.io/biocontainers/icount:2.0.0--py_1
fi

# ---------------------------
# Download and zip fastq
# ---------------------------
fastq files have been downloaded

# cd $MAIN_DIR/data
# fasterq-dump SRR30639402 SRR30639403 SRR30639404 SRR30639405 SRR30639406 SRR30639407
# gzip *

# ---------------------------
# fastqc for all
# ---------------------------
cd $MAIN_DIR/data
fastqc *.fastq.gz -o $OUT_DIR/fastqc/

python3 /data/jling2/bin/email/sendemail.py -s "${EXP_NAME}: clip_tdp43" -c "clip_tdp43 fastqc complete" -r isinha1@jh.edu

# ---------------------------
# umi_tools: Extract the UMIs
# ---------------------------
cd $MAIN_DIR/data

for file in *.fastq.gz; do
  singularity exec $SING_IMAGES/umi_tools.sif umi_tools extract --bc-pattern=NNNNNNNNN --stdin=$file --log=$OUT_DIR/UMIextracted/${file%.fastq.gz}_processed.log --stdout $OUT_DIR/UMIextracted/${file%.fastq.gz}_UMI.fastq.gz
done


# ---------------------------
# Read trimming using Trim Galore!
# Use universal Illumina Adapter as IDed by FastQC
# Make new folder named "Trimmed"
# ---------------------------
cd $OUT_DIR/UMIextracted

for file in *.fastq.gz; do
  singularity exec $SING_IMAGES/trim-galore.sif trim_galore -j 6 --fastqc --dont_gzip -o "${OUT_DIR}/Trimmed/" $file
done

python3 /data/jling2/bin/email/sendemail.py -s "${EXP_NAME}: clip_tdp43" -c "clip_tdp43 Trim Galore! complete" -r isinha1@jh.edu

# ---------------------------
# STAR for mapping to genome
# ---------------------------
ml load STAR/2.7.10a samtools/1.10

# Generate index for STAR if needed
if [ -e "${GENOME_DIR}/indices" ]; then
  echo "genome indexing completed previously"
else
  mkdir indices
  STAR --runMode genomeGenerate --runThreadN 6 --genomeDir $GENOME_DIR/indices --genomeFastaFiles $GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile $GENOME_DIR/Homo_sapiens.GRCh38.115.gtf --sjdbOverhang 100
  python3 /data/jling2/bin/email/sendemail.py -s "clip_tdp43" -c "index finished running" -r isinha1@jh.edu
fi

# Align reads
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
cd $OUT_DIR/STARoutput
ml load singularity/3.8.7

for file in *_Aligned.sortedByCoord.out.bam; do \
  singularity exec $SING_IMAGES/umi_tools.sif umi_tools dedup -I $file --output-stats=deduplicated -S $OUT_DIR/UMIdedup/${file%_UMI_trimmed_Aligned.sortedByCoord.out.bam}_AlignSortDedup.bam
done

cd $OUT_DIR/UMIdedup
for file in *_AlignSortDedup.bam; do \
  samtools index -@ 6 $file
done

## Basic tools complete
python3 /data/jling2/bin/email/sendemail.py -s "${EXP_NAME}: clip_tdp43" -c "clip_tdp43 STAR and dedup finished running" -r isinha1@jh.edu

# ---------------------------
# iCount: Quantifying cross-linked sites
# ---------------------------
cd $OUT_DIR/UMIdedup

for file in *_AlignSortDedup.bam; do \
  singularity exec $SING_IMAGES/icount.sif iCount xlsites $file \
  $OUT_DIR/iCount/${file%_AlignSortDedup.bam}_cDNA_unique.bed \
  $OUT_DIR/iCount/${file%_AlignSortDedup.bam}_cDNA_multiple.bed \
  $OUT_DIR/iCount/${file%_AlignSortDedup.bam}_cDNA_skipped.bam \
  --group_by start --quant cDNA
done

# ---------------------------
# iCount: protein-RNA interaction analysis
# https://icount.readthedocs.io/en/latest/tutorial.html
# ---------------------------
# Segment genome if necessary
if [ -e "${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai" ]; then
  if [ -e "${GENOME_DIR}/hs_GRCh38_seg_iCLIP.gtf.gz" ]; then
    echo "genome segmenting for iCount complete"
  else
    singularity exec $SING_IMAGES/icount.sif iCount segment $GENOME_DIR/Homo_sapiens.GRCh38.115.gtf $GENOME_DIR/hs_GRCh38_seg_iCLIP.gtf.gz $GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai;
  fi
else
  samtools faidx $GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa
  singularity exec $SING_IMAGES/icount.sif iCount segment $GENOME_DIR/Homo_sapiens.GRCh38.115.gtf $GENOME_DIR/hs_GRCh38_seg_iCLIP.gtf.gz $GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
fi

# Identifying significantly cross-linked sites
cd $OUT_DIR/iCount
python3 /data/jling2/bin/email/sendemail.py -s "${EXP_NAME}: clip_tdp43" -c "clip_tdp43 start iCount " -r isinha1@jh.edu

for file in *_cDNA_unique.bed; do \
  singularity exec $SING_IMAGES/icount.sif iCount peaks $GENOME_DIR/hs_GRCh38_seg_iCLIP.gtf.gz \
  $file ${file%_cDNA_unique.bed}_peaks.bed \
--scores ${file%_cDNA_unique.bed}_scores.tsv
done

# Identifying clusters of significantly cross-linked sites
for file in *_peaks.bed; do \
  singularity exec $SING_IMAGES/icount.sif iCount clusters $file ${file%_peaks.bed}_clusters.bed
done

# Annotating sites and summary statistics
for file in *_cDNA_unique.bed; do \
  singularity exec $SING_IMAGES/icount.sif annotate $GENOME_DIR/hs_GRCh38_seg_iCLIP.gtf.gz $file ${file%_cDNA_unique.bed}_annotated_sites_biotype.tab
  singularity exec $SING_IMAGES/icount.sif annotate --subtype gene_id $GENOME_DIR/hs_GRCh38_seg_iCLIP.gtf.gz $file ${file%_cDNA_unique.bed}_annotated_sites_genes.tab
  singularity exec $SING_IMAGES/icount.sif summary $GENOME_DIR/hs_GRCh38_seg_iCLIP.gtf.gz $file ${file%_cDNA_unique.bed}_summary.tab $GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
done

# ---------------------------
# Ending the run
# ---------------------------
python3 /data/jling2/bin/email/sendemail.py -s "${EXP_NAME}: clip_tdp43" -c "clip_tdp43 iCount finished running" -r isinha1@jh.edu
exit 0
