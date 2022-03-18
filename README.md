# TSAannotate
This code annotates the peptides produced by the core algorithm of TSAFinder

# Rerequisites
# R packages:
matrixStats \
data.table \
reshape2 \
biomaRt 
# Software:
Blast+ (tested with version 2.9.0) \
bedtools (tested with version 2.29.2)

# Required data:
A reference gtf file: tested with gencode.v22.annotation.gtf \
Reference genome (grch38, assembly 22) \
Coding reading frame dataset (CCDS, GRCg38.p12 assembly 22, release number 2018614) \

# Usage
rnaneoantigen.master_2.R is the wrapper for the annotate and loader functions.

1. Enter rnaneoantigen.master_2.R and update the wd and wd.in variables. 
  - "wd" is the directory to your supporting files (reference genomes, etc), the temporary files creating during code execution, and the code output. If you are downloading the files from the provided dropbox link above, just replace the "..." with the path to the downloaded dropbox file.
  - "wd.in" is the directory to the output from your TSAcore run. There should be a "sampleid_kmer_neoantigens.txt" and a "sample_id_kmer_neoantmetadata.txt" file. Note that you MUST use the "getmeta" option in your TSAcore run to have the neoantmetadata.txt file ready for TSAannotate.
2. Specify the patient character vector, either via command line or by inputting a character vector of the prefixes. Make sure to uncomment the correct line in the rnaneoantigen.master_2.R file.
3. Change the peplength argument if needed. Default is 8-11 \


**Output** is a spreadsheet of peptides in rows that mapped to the reference genome and their annotations in columns, including genomic location, variants (missense and IN/DELs), gene/exon, and coding reference frame. The coding reading frame is not tested.
