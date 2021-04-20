######################## rnaneoantigen.master.R ################################
#library(data.table,lib.loc='/users/PAS1203/osu1042/R/x86_64-pc-linux-gnu-library/3.5')
library(matrixStats) #,lib.loc='/users/PAS1203/osu1042/R/x86_64-pc-linux-gnu-library/3.5')
library(reshape2) #,lib.loc='/users/PAS1203/osu1042/R/x86_64-pc-linux-gnu-library/3.5')
library(biomaRt) #,lib.loc='/users/PAS1203/osu1042/R/x86_64-pc-linux-gnu-library/3.5')

wd <- '/wynton/home/jones/sharpnacmi/rnaneoantigen/'
setwd(wd)
load('RData/master_input.RData')
source('rnaneoantigen.loader_2.R')
source('rnaneoantigen.annotate.R')
#wd <- '/wynton/home/jones/sharpnacmi/rnaneoantigen/'
#gtf = gtf.load(wd)
#ccds = framed.load(wd)
#args = commandArgs(trailingOnly=TRUE)
#patient <- toString(args[1])
peplengths <- c(8,9,10,11)
#patients <- c('07H103','10H080','10H118','12H018','luc2','luc4','luc6')
patients <- c('SRR10780179','SRR10783095','SRR10783852','SRR10794775')
for(patientcounter in 1:1){ #length(patients)){
  print(paste('Processing patient',patientcounter,'of',length(patients)))
  patient <- patients[patientcounter]
  ptm2 = proc.time()
  for(pepcounter in 1:length(peplengths)){
    print(paste('Peplength =',peplengths[pepcounter]))
    ptm = proc.time()
    peplength <- peplengths[pepcounter]
    #load(paste('RData/netmhc_load/netmhc_load.',patient,'_',peplength,'.RData',sep=''))
    #load(paste('../laumont_netmhc/rnaneoantigen.master',patient,'_',peplength,'.RData',sep=''))
    netmhc <- netmhc.loader(peplength,patient,wd)
    blastout <- blast(patient,peplength,wd,max.mutations = 2,max.gaps = 1,onlytophits = TRUE)
    blastout.fasta <- fastafromblast(patient,wd,blastout,netmhc[[2]],netmhc[[3]],peplength)
    blastout <- junctioncaller(blastout.fasta,wd,patient,netmhc[[2]],netmhc[[3]],peplength)
    blastout.fasta.variant <- variantcaller(patient,wd,blastout,peplength,max.mutations = 2, max.gaps = 1)
    blastout.fasta.variant.gtf <- gtf.blast(blastout.fasta.variant,gtf[[1]],gtf[[2]],gtf[[3]])
    blastout.fasta.variant.gtf.framed <- framed(ccds.chrom,blastout.fasta.variant.gtf)
    save(blastout,blastout.fasta,blastout.fasta.variant,blastout.fasta.variant.gtf,blastout.fasta.variant.gtf.framed, file = paste('RData/master/rnaneoantigen.master',patient,'_',peplength,'.RData',sep=''))
    print('Analyzing this peplength took')
    print(proc.time()-ptm)
  }
  print('Analyzing this patient took')
  print(proc.time()-ptm2)
}


