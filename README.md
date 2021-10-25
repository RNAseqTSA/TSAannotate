# TSAannotate
This code annotates the peptides produced by the core algorithm of TSAFinder

# Rerequisites
# R packages:
matrixStats
data.table
reshape2
biomaRt
# Software:
Blast+ (tested with version 2.9.0) \
bedtools (tested with version 2.29.2)

# Required data:
A reference gtf file: tested with gencode.v22.annotation.gtf
Reference genome (grch38, assembly 22)
Coding reading frame dataset (CCDS, GRCg38.p12 assembly 22, release number 2018614)

# Usage
rnaneoantigen.master_2 is the main code

# Inputs: 
change the inputs and output directories within the rnaneoantigen.master_2 code
Change the peplength argument if needed. Default is 8-11
command line arguments include the sample name

#Output is a spreadsheet of peptides in rows that mapped to the reference genome and their annotations in columns, including genomic location, variants (missense and IN/DELs), gene/exon, and coding reference frame. The coding reading frame is not tested.
