######################## rnaneoantigen.master.R ################################
library(matrixStats)
library(data.table)
library(reshape2) 
library(biomaRt) 

#wd.in is the location of your TSAcore output. Note that it must end in a '/'
wd.in <- ''
#wd is the location of your TSAannotate supporting files (files available for public download on dropbox). Note that it must end in a '/'
wd <- '.../rnaneoantigen/'
bedtools.wd <- NULL 
#extra normal database for filtering peptides
normal.wd <- paste(wd,'references/normalpep_sep/',sep="")
setwd(wd)
load('references/master_input_07272021.RData')
source('TSAannotate/rnaneoantigen.loader_2.R')
source('TSAannotate/rnaneoantigen.annotate.R')

#only uncomment the below two lines if you are NOT loading master_input_07272021
#gtf = gtf.load(wd)
#ccds.chrom = framed.load(wd)

#uncomment the following two lines if you are running R from command line, with the patient IDs as the first argument
#args = commandArgs(trailingOnly=TRUE)
#patients <- toString(args[1])

#patients is the vector of your file name prefixes in your TSAcore output

peplengths <- c(8,9,10,11)
for(patientcounter in 1:length(patients)){
  print(paste('Processing patient',patientcounter,'of',length(patients)))
  patient <- patients[patientcounter]
  ptm2 = proc.time()
  for(pepcounter in 1:length(peplengths)){
    print(paste('Peplength =',peplengths[pepcounter]))
    ptm = proc.time()
    peplength <- peplengths[pepcounter]
    print('Loading and prepping netMHC output and generate contigs')
    netmhc <- netmhc.loader(peplength,patient,wd,wd.in,normal.wd)
    print("Perform nucleotide blast alignment of contigs")
    blastout <- blast(patient,peplength,wd,max.mutations = 2,max.gaps = 1,onlytophits = TRUE,netmhc[[2]])
    print("Extract fasta reference sequences for peptides")
    blastout.fasta <- fastafromblast(patient,wd,blastout,netmhc[[2]],netmhc[[3]],peplength,bedtools.wd)
    print('Add junction-spanning reads')
    blastout <- junctioncaller(blastout.fasta,wd,patient,netmhc[[2]],netmhc[[3]],peplength,bedtools.wd)
    print('Annotate variants')
    blastout.fasta.variant <- variantcaller(patient,wd,blastout,peplength,max.mutations = 2, max.gaps = 1)
    print('Add gene annotations')
    blastout.fasta.variant.gtf <- gtf.blast(blastout.fasta.variant,gtf[[1]],gtf[[2]],gtf[[3]])
    print('Add frame annotations')
    blastout.fasta.variant.gtf.framed <- framed(ccds.chrom,blastout.fasta.variant.gtf)
    save(blastout,blastout.fasta,blastout.fasta.variant,blastout.fasta.variant.gtf,blastout.fasta.variant.gtf.framed, file = paste('RData/master/rnaneoantigen.master',patient,'_',peplength,'.RData',sep=''))
    print('Analyzing this peplength took')
    print(proc.time()-ptm)
  }
  print('Analyzing this patient took')
  print(proc.time()-ptm2)
}


