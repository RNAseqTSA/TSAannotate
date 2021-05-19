######################## rnaneoantigen.master.R ################################
library(matrixStats)
library(reshape2) 
library(biomaRt) 

#wd.in is the location of your TSAcore output. Note that it must end in a '/'
wd.in <- '/home/michaelsharpnack/Documents/PRJNA591860/LT_S01/'
#wd is the location of your TSAannotate supporting files. Note that it must end in a '/'
wd <- '/home/michaelsharpnack/Documents/rnaneoantigen/'
bedtools.wd <- '/home/michaelsharpnack/local/bedtools2/bin/'
#extra normal database for filtering peptides
normal.wd <- '/home/michaelsharpnack/Documents/rnaneoantigen/references/normalpep_sep/'
setwd(wd)
load('RData/master_input.RData')
source('TSAannotate/rnaneoantigen.loader_2.R')
source('TSAannotate/rnaneoantigen.annotate.R')
#wd <- '/wynton/home/jones/sharpnacmi/rnaneoantigen/'
#gtf = gtf.load(wd)
#ccds = framed.load(wd)
#args = commandArgs(trailingOnly=TRUE)
#patient <- toString(args[1])
peplengths <- c(8,9,10,11)
#patients <- c('07H103','10H080','10H118','12H018','luc2','luc4','luc6')
patients <- c('SRR10783095')
for(patientcounter in 1:1){ #length(patients)){
  print(paste('Processing patient',patientcounter,'of',length(patients)))
  patient <- patients[patientcounter]
  ptm2 = proc.time()
  for(pepcounter in 2:length(peplengths)){
    print(paste('Peplength =',peplengths[pepcounter]))
    ptm = proc.time()
    peplength <- peplengths[pepcounter]
    #load(paste('RData/netmhc_load/netmhc_load.',patient,'_',peplength,'.RData',sep=''))
    #load(paste('../laumont_netmhc/rnaneoantigen.master',patient,'_',peplength,'.RData',sep=''))
    print('Loading and prepping netMHC output and generate contigs')
    netmhc <- netmhc.loader(peplength,patient,wd,wd.in,normal.wd)
    print("Perform nucleotide blast alignment of contigs")
    blastout <- blast(patient,peplength,wd,max.mutations = 2,max.gaps = 1,onlytophits = TRUE)
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


