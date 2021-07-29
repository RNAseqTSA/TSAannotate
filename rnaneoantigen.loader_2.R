################################ rnaneoantigen.loader.R   Load in files  ########################################################

peptide.to.seq <- function(peplength,netmhc_out,pepmet,patient,wd,normal.wd){
  #this function creates a list of contigs for each peptide. Steps:
  #1. find sub-sequence that codes for peptide within rnaseq read
  #2. Center reads on peptide coding sequence (maximum possible contig length)
  #2. Use results from #1 as seed to add on reads
  
  #create pep to codon table
  codon.table <- list()
  codon.table[["F"]] <- c('TTC','TTT')
  codon.table[['L']] <- c('CTA','CTC','CTG','CTT','TTA','TTG')
  codon.table[['I']] <- c('ATA','ATC','ATT')
  codon.table[['M']] <- 'ATG'
  codon.table[['V']] <- c('GTA','GTC','GTG','GTT')
  codon.table[['S']] <- c('TCA','TCC','TCG','TCT','AGT','AGC')
  codon.table[['P']] <- c('CCA','CCC','CCG','CCT')
  codon.table[['T']] <- c('ACA','ACC','ACG','ACT')
  codon.table[['A']] <- c('GCA','GCC','GCG','GCT')
  codon.table[['Y']] <- c('TAC','TAT')
  codon.table[['H']] <- c('CAT','CAC')
  codon.table[['Q']] <- c('CAA','CAG')
  codon.table[['N']] <- c('AAT','AAC')
  codon.table[['K']] <- c('AAA','AAG')
  codon.table[['D']] <- c('GAT','GAC')
  codon.table[['E']] <- c('GAA','GAG')
  codon.table[['C']] <- c('TGT','TGC')
  codon.table[['W']] <- 'TGG'
  codon.table[['R']] <- c('CGA','CGC','CGG','CGT','AGA','AGG')
  codon.table[['G']] <- c('GGA','GGC','GGG','GGT')
  
  bind.names <- rownames(netmhc_out[rowMins(as.matrix(netmhc_out[,1:(dim(netmhc_out)[2]-1)]),na.rm=TRUE) < 50,])
  pepmet <- pepmet[(which(pepmet[[1]] %in% bind.names)),]
  pepmet <- pepmet[order(pepmet[[1]]),]
  pepmet <- pepmet[,c(2:1)]
  
  pep.contig <- vector(mode='character',length(bind.names))
  names(pep.contig) <- bind.names
  peponly.contig <- pep.contig
  contig.pepstart <- pep.contig
  #loop through peptides
  for(y in 1:length(bind.names)){
    #unique rna reads corresponding for peptide y
    reads.temp <- pepmet[[1]][which(pepmet[[2]] == bind.names[y])]
    unique.reads.freq <- table(reads.temp)
    unique.reads <- as.character(unique(reads.temp))
    
    pep.rna.seq <- matrix(NA,nrow=length(unique.reads),ncol=4)  #RNAseq locations and sequences for each peptide
    #figure out read corresponding to each peptide
    for(o in 1:length(unique.reads)){
      n = 1
      for(z in 1:nchar(unique.reads[o])){
        #print(paste('analyzing RNA read position',z))
        if(substr(unique.reads[o],z,z+2) %in% codon.table[[substr(bind.names[y],n,n)]]){
          n.temp <- n+1
          for(m in n:(peplength-1)){
            if(substr(unique.reads[o],z+m*3,z+m*3+2) %in% codon.table[[substr(bind.names[y],n.temp,n.temp)]]){
              if(m < peplength-1){n.temp <- n.temp+1}
            } else {
              #print('subsequent peptides do no match, move down RNA read')
              break
            }
          }
          if(n.temp == peplength){
            pep.rna.seq[o,] <- c(z,substr(unique.reads[o],z,z+peplength*3-1),unique.reads[o],unique.reads.freq[o])
            break
          }
        }
      }
    }
    pep.rna.seq <- pep.rna.seq[order(as.numeric(pep.rna.seq[,1])),]
    if(is.null(dim(pep.rna.seq)) == FALSE){
      read.length.max <- max(nchar(pep.rna.seq[,3]),na.rm=TRUE)
      contig.length <- peplength*3+(read.length.max-peplength*3)*2
      pep.rna.seq.full <- matrix(NA,nrow=dim(pep.rna.seq)[1],ncol=contig.length)
      for(counter in 1:dim(pep.rna.seq.full)[1]){
        read.length <- nchar(pep.rna.seq[counter,3])
        pep.rna.seq.full[counter,(read.length.max-peplength*3-as.numeric(pep.rna.seq[counter,1])+2):(read.length.max-peplength*3-as.numeric(pep.rna.seq[counter,1])+1+read.length)] <- strsplit(pep.rna.seq[counter,3],'')[[1]]
      }
      pep.rna.seq.unique <- list()
      for(counter in 1:dim(pep.rna.seq.full)[2]){
        pep.temp <- cbind(pep.rna.seq[,4],pep.rna.seq.full[,counter])[is.na(pep.rna.seq.full[,counter]) == FALSE,]
        if(is.null(dim(pep.temp))){
          pep.rna.seq.unique[[counter]] <- as.numeric(pep.temp[1])
          names(pep.rna.seq.unique[[counter]]) <- pep.temp[2]
        } else {
          pep.rna.seq.unique[[counter]] <- vector(mode='numeric',length(unique(pep.temp[,2])))
          names(pep.rna.seq.unique[[counter]]) <- unique(pep.temp[,2])
          for(counter2 in 1:length(pep.rna.seq.unique[[counter]])){
            pep.rna.seq.unique[[counter]][counter2] <- sum(as.numeric(pep.temp[pep.temp[,2] == unique(pep.temp[,2])[counter2],1]))
          }
        }
      }
    } else {
      read.length.max <- nchar(pep.rna.seq[3])
      contig.length <- peplength*3+(read.length.max-peplength*3)*2
      pep.rna.seq.full <- vector('character',length=contig.length)
      read.length <- nchar(pep.rna.seq[3])
      pep.rna.seq.full[(read.length.max-peplength*3-as.numeric(pep.rna.seq[1])+2):(read.length.max-peplength*3-as.numeric(pep.rna.seq[1])+1+read.length)] <- strsplit(pep.rna.seq[3],'')[[1]]
      pep.rna.seq.unique <- list()
      for(counter in 1:length(pep.rna.seq.full)){
        pep.temp <- c(pep.rna.seq[4],pep.rna.seq.full[counter])[is.na(pep.rna.seq.full[counter]) == FALSE]
        pep.rna.seq.unique[[counter]] <- as.numeric(pep.temp[1])
        names(pep.rna.seq.unique[[counter]]) <- pep.temp[2]
      }
    }
    #create the contig with winner takes all approach
    contig  <- vector(mode='character',length(pep.rna.seq.unique))
    for(counter in 1:length(pep.rna.seq.unique)){
      contig[counter] <- names(which.max(pep.rna.seq.unique[[counter]]))
    }
    peponly.contig[y] <- paste(contig[(read.length.max-peplength*3+1):(read.length.max)],collapse='')
    contig.pepstart[y] <- read.length.max-peplength*3+1-sum(nchar(contig[1:read.length.max]) == 0)
    pep.contig[y] <- paste(contig,collapse = '')
    
  }
  contig.write <- vector(mode='character',length(pep.contig)*2)
  for(i in 1:length(pep.contig)){
    contig.write[2*i-1] <- paste(">",names(pep.contig)[i],sep='')
    contig.write[2*i] <- pep.contig[i]
  }
  #for(contig.counter in 1:ceiling(length(contig.write)/ceiling(length(contig.write))){
  #  if(contig.counter < ceiling(length(contig.write)/ceiling(length(contig.write))){
  #    write.table(contig.write[((contig.counter-1)*ceiling(length(contig.write)+1):(contig.counter*100)],file=paste(wd,'contig/',patient,peplength,'_',contig.counter,sep=''),quote=FALSE,row.names = FALSE,col.names=FALSE)
  #  } else {
  #    write.table(contig.write[((contig.counter-1)*1+1):length(contig.write)],file=paste(wd,'contig/',patient,peplength,'_',contig.counter,sep=''),quote=FALSE,row.names = FALSE,col.names=FALSE)
  #  }
  #}
  write.table(contig.write,file=paste(wd,'contig/',patient,peplength,sep=''),quote=FALSE,row.names = FALSE,col.names=FALSE)
  pep.contig.save <- cbind(pep.contig,peponly.contig,contig.pepstart)
  rownames(pep.contig.save) <- names(pep.contig)
  
  return(list(contig.write,pep.contig.save))
}


netmhc.loader <- function(peplength,patient,wd,wd.in,normal.wd){
  #col.names <- c('Pos','HLA','Peptide','Core','Of','Gp','Gl','Ip','ll','lcore','Identity','Score','Aff(nM)','%Rank','BindLevel')
  contig.return = TRUE
  #load netmhcpan output
  #netmhc <- read.table(paste(wd.in,patient,"_",peplength,"mer_neoantigens.txt",sep=""),header=FALSE,fill=TRUE,col.names = col.names,sep=',')
  #netmhc <- netmhc[nchar(as.character(netmhc[[3]])) == peplength,]
  #netmhc <- data.frame(netmhc[[3]],netmhc[[2]],netmhc[[13]])
  #colnames(netmhc) <- c('V1','V2','V3')
  #netmhc <- netmhc[duplicated(netmhc) == FALSE,]
  #netmhc <- dcast(netmhc,V1 ~ V2,value.var = 'V3')
  #rownames(netmhc) <- netmhc[[1]]
  #netmhc <- netmhc[,-1]
  #netmhc <- netmhc[-grep('X',rownames(netmhc)),]
 
  netmhc <- fread(paste(wd.in,patient,"_",peplength,"mer_neoantigens.txt",sep=""),header=FALSE)
  netmhc <- netmhc[nchar(as.character(netmhc[[4]])) == peplength,]
  netmhc <- data.frame(netmhc[[4]],netmhc[[2]],netmhc[[5]])
  colnames(netmhc) <- c('V1','V2','V3')
  netmhc <- netmhc[duplicated(netmhc) == FALSE,]
  netmhc <- dcast(netmhc,V1 ~ V2,value.var = 'V3')
  rownames(netmhc) <- netmhc[[1]]
  netmhc <- netmhc[,-1]
     
  #load in peptide counts
  pepmet <- read.table(paste(wd.in,patient,'_',peplength,'mer_neoantmetadata.txt',sep=""))
  pep.count <- table(pepmet[[1]])
  if(length(pep.count) != nrow(netmhc)){print('warning: number of peptides in metadata and neoantigens run through netmhcpan are not equal')}
  pep.count <- pep.count[intersect(names(pep.count),rownames(netmhc))]
  netmhc <- netmhc[intersect(names(pep.count),rownames(netmhc)),]
  netmhc <- cbind(netmhc,as.vector(pep.count))

  if(contig.return == TRUE){
    contig <- peptide.to.seq(peplength,netmhc,pepmet,patient,wd,normal.wd)
    return(list(netmhc,contig[[1]],contig[[2]]))
  } else {
    return(netmhc)
  }
}



