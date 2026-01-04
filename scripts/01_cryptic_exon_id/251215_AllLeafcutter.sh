#!/bin/bash

#SBATCH --job-name=eclip_tdp43
#SBATCH --time=72:00:00
#SBATCH --partition=parallel
#SBATCH --ntasks=12
#SBATCH --mem=60G

# ---------------------------
# Define variables
# ---------------------------
MAIN_DIR="/data/jling2/irika/processedRNAseq/OurRuns/2512_NMDi_Leafcutter"
INPUT_DIR="${MAIN_DIR}/input"
OUT_DIR="${MAIN_DIR}/output/$(date +"%Y%m%d_%H%M%S")"
BIN_DIR="/data/jling2/irika/bin/"
mkdir -p $OUT_DIR ${OUT_DIR}/group{1..9}

# ---------------------------
# Converting bams to juncs
# ---------------------------
# # This step is complete. Used samtools.
#
# # Modules
# ml load samtools
# ml load GCC/11.3.0 GCC/12.2.0 GCC/12.3.0 GCC/13.2.0 foss/2022a foss/2022b foss/2023a gfbf/2022b gfbf/2023a gfbf/2023b intel/2019a intel/2020a intel/2023a intel/2023b intel/2024a RegTools/1.0.0-foss-2022a
#
# for file in *.bam; do \
#   echo Converting $i to $i.junc
#   samtools index $i # Index (bam.bai)
#   regtools junctions extract -a 8 -m 50 -M 500000 $i -o $i.junc -s XS
#   echo $i.junc >> NMDi_juncfiles.txt # Add file name to txt file
#
# done

# ---------------------------
# Intron clustering (Python)
# ---------------------------
ml load anaconda gcc/9.3.0 r/4.0.2 r-optparse/1.6.2 GCC/11.3.0 GCC/12.2.0 GCC/12.3.0 GCC/13.2.0 foss/2022a foss/2022b foss/2023a gfbf/2022b gfbf/2023a gfbf/2023b intel/2019a intel/2020a intel/2023a intel/2023b R/4.2.2-foss-2022b leafcutter/0.2.9-foss-2022b-R-4.2.2

cd ${INPUT_DIR}
python ${BIN_DIR}/leafcutter/clustering/leafcutter_cluster_regtools.py -j NMDi_juncfiles.txt -m 50 -o NMDi -l 500000
python /data/jling2/bin/email/sendemail.py -s "leafcutter -  intron clustering" -c "intron clustering complete" -r isinha1@jh.edu

# ---------------------------
# Differential intron excision analysis
# ---------------------------
ml load anaconda gcc/9.3.0 r/4.0.2 r-optparse/1.6.2 GCC/11.3.0 GCC/12.2.0 GCC/12.3.0 GCC/13.2.0 foss/2022a foss/2022b foss/2023a gfbf/2022b gfbf/2023a gfbf/2023b intel/2019a intel/2020a intel/2023a intel/2023b R/4.2.2-foss-2022b leafcutter/0.2.9-foss-2022b-R-4.2.2

for i in {1..9}; do \
  # Set file names
  fName="groups_${i}.txt"

  # Move to output folder
  cd ${OUT_DIR}/group${i}

  # Run LC
  /data/jling2/irika/bin/leafcutter/scripts/leafcutter_ds.R \
    --num_threads 4 ${INPUT_DIR}/NMDi_perind_numers.counts.gz ${MAIN_DIR}/groups_file/${fName} -i 1 -g 1

done

python /data/jling2/bin/email/sendemail.py -s "NMDi LC" -c "splicing leafcutter finished running" -r isinha1@jh.edu

# ---------------------------
# Kill compute
# ---------------------------
exit

