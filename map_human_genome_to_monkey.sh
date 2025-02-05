#!/bin/bash
# Author : Kedar Joshi

# Script: map_human_genome_to_monkey.sh
# Description: Mapping a part of the human genome to the single-cell sequenced adult monkey genome.
# Requirements: Cell Ranger, wget, and Seurat (for downstream analysis in R).

# Step 1: Download the part of the human genome (DUX4 region) in FASTA format
echo "Downloading the DUX4 region of the human genome (chr4) in FASTA format..."
wget -O human_chr4_dux4.fa "https://www.ncbi.nlm.nih.gov/nuccore/NC_000004.12?report=fasta&from=190173774&to=190185911"

# Step 2: Download the GTF file for creating a reference genome
echo "Downloading the GTF file for human chr4..."
wget -O human_chr4.gtf "<insert_gtf_file_url_here>" # Replace with the actual GTF file URL

# Step 3: Download and install cellranger
curl -o cellranger-9.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.0.tar.gz?Expires=1734406683&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=cy7rzUpPPSZLBSHcjO5JjJUAcHJ~9IVOqMh1aJqZN8ucFjPq6pE39TOGFWinbjb7xAa4Q4XJRbsLph9zUnlPu9UFRkcWQafH3aC7XFzWc~fyb-3Ok0ivQfLJ6j-jzGwMXJTXLzNiXMZBTwYPWkzK0Glge2Jsa~kSmHVx-fAunqlyw-8AVq4GjsjRFirYO72jIyqR~H22FjtxNcP4pzdHQ9cAw3zsAadrbxaALy6RargQWgkxeIzSLOPgnLwj17qDK-sS7uNSibIJNWA9245f~rLPIOwc1ekunhPY6zmgHb9MEeS5YC2SSC1bjnTquBcAaW9tVj2182qacajGGpRUpA__"
tar -xzvf cellranger-x.y.z.tar.gz
export PATH=/opt/cellranger-x.y.z:$PATH
which cellranger

# Step 4: Create the reference genome using Cell Ranger's 'mkref'
echo "Creating the reference genome with Cell Ranger..."
cellranger mkref \
    --genome=human_chr4_dux4 \
    --fasta=human_chr4_dux4.fa \
    --genes=human_chr4.gtf

# Step 5: Generate the count matrix using Cell Ranger 'count'
echo "Running Cell Ranger count to generate the count matrix..."
cellranger count \
    --id=Dux_count \
    --transcriptome=human_chr4_dux4 \
    --fastqs=/efs/home/kjoshi \
    --sample=Adult1 \
    --localcores=16 \
    --localmem=128

# Step 6: Use Seurat for visualizing the genes (to be done in R)
echo "Count matrix generated. Proceed to downstream analysis using Seurat in R."
