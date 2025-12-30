#!/bin/bash

#SBATCH --job-name=eclip_tdp43_Tam
#SBATCH --time=72:00:00
#SBATCH --partition=parallel
#SBATCH --ntasks=6
#SBATCH --mem-per-cpu=100GB

# ---------------------------
# Define variables
# ---------------------------
MAIN_DIR="/data/jling2/irika/clipdata/encode_HepG2_ENCSR187VEQ"
EXPERIMENT_NAME="encode_HepG2_ENCSR187VEQ"
GENOME_DIR="/data/jling2/irika/bin/genomes/GRCh38_v115"
OUT_DIR="${MAIN_DIR}/output/$(date +"%Y%m%d_%H%M%S")"
SING_DIR="$HOME/singularity"
SING_IMAGES="$SING_DIR/singularity_images"
export SINGULARITY_TMPDIR=$SING_DIR/tmp
export SINGULARITY_CACHEDIR=$SING_DIR/singularity_cache
mkdir -p $SING_DIR $SING_IMAGES $SINGULARITY_TMPDIR $SINGULARITY_CACHEDIR $OUT_DIR $MAIN_DIR/data $OUT_DIR/UMIdedup

# ---------------------------
# Load modules
# ---------------------------
ml load sra-tools/3.0.3 singularity/3.8.7 fastqc/0.11.9
python3 /data/jling2/bin/email/sendemail.py -s "clip_tdp43 ${EXPERIMENT_NAME}" -c "${EXPERIMENT_NAME} clip_tdp43 START" -r isinha1@jh.edu

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
# Sort known BAM files
# ---------------------------
# ml load samtools/1.10

# BAM files have been sorted and indexed

# Sort using samtools
# cd $MAIN_DIR/data

# for file in *.bam; do \
#  samtools sort -@ 6 -o ${file%.bam}_Aligned.sortedByCoord.out.bam $file
# done

# for file in *_Aligned.sortedByCoord.out.bam; do \
#   samtools index -@ 6 $file
# done

# ---------------------------
# UMI deduplication
# ---------------------------
cd $MAIN_DIR/data
ml load singularity/3.8.7

for file in *_Aligned.sortedByCoord.out.bam; do \
  singularity exec $SING_IMAGES/umi_tools.sif umi_tools dedup -I $file --output-stats=deduplicated -S $OUT_DIR/UMIdedup/${file%_Aligned.sortedByCoord.out.bam}_AlignSortDedup.bam
done

cd $OUT_DIR/UMIdedup
for file in *_AlignSortDedup.bam; do \
  samtools index -@ 6 $file
done

## Basic tools complete
python3 /data/jling2/bin/email/sendemail.py -s "${EXPERIMENT_NAME} clip_tdp43" -c "${EXPERIMENT_NAME} clip_tdp43 fastq basic tools finished running" -r isinha1@jh.edu

# ---------------------------
# iCount: Quantifying cross-linked sites
# ---------------------------
mkdir $OUT_DIR/iCount
cd $OUT_DIR/UMIdedup

for file in *_AlignSortDedup.bam; do \
  singularity exec $SING_IMAGES/icount.sif iCount xlsites $file \
  ${file%_AlignSortDedup.bam}_cDNA_unique.bed  ${file%_AlignSortDedup.bam}_cDNA_multiple.bed ${file%_AlignSortDedup.bam}_cDNA_skipped.bam \
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
python3 /data/jling2/bin/email/sendemail.py -s "clip_tdp43" -c "clip_tdp43 start iCount " -r isinha1@jh.edu

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
python3 /data/jling2/bin/email/sendemail.py -s "${EXPERIMENT_NAME} CLIP-analysis" -c "${EXPERIMENT_NAME} clip_tdp43 iCount finished running" -r isinha1@jh.edu
exit 0
