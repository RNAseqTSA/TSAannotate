library(dplyr)
library(reshape2)
`%notin%` <- Negate(`%in%`)

#change to the directory containing the tsa object
X = read.table('~/Documents/tsa',sep='\t',header=TRUE)

X2 <- X %>%
  group_by(Patient.ID, qseqid) %>%
  tally()

X3 <- X2 %>%
  group_by(qseqid) %>%
  tally() %>%
  arrange(-n)

common_TSA = as.matrix(X3[X3$n>=10,"qseqid"])
TSA <- acast(X2[X2$qseqid %in% common_TSA,],Patient.ID~qseqid,value.var="n")
TSA <- rbind(TSA,matrix(NA,sum(unique(X2$Patient.ID) %notin% row.names(TSA)),dim(TSA)[2],dimnames=list(unique(X2$Patient.ID)[unique(X2$Patient.ID) %notin% row.names(TSA)],colnames(TSA))))          
TSA[is.na(TSA)] <- 0

setCover <- function(sets){
  edges = c()
  sets = as.data.frame(sets)
  colnames(sets) = c('node','edge')
  unused = unique(sets$node)
  while(length(unused) > 0){
    sets = sets[sets[,1] %in% unused,]
    #message(length(unique(sets[sets[,2]=="TPSKTRTSTL","Patient.ID"])))
    sets_count <- sets %>%
      group_by(edge) %>%
      tally() %>%
      arrange(-n)
    if(sets_count$n[1]<2){break}
    sets_count = as.data.frame(sets_count)
    edges = c(edges,sets_count$edge[1])
    unused = setdiff(unused,unique(sets$node[sets$edge==sets_count$edge[1]]))
  }
  message('Missing:')
  message(paste0(unused,collapse=', '))
  return(edges)
}

# Optimal TSA (all TSA)
optimal_TSA = setCover(X2[,1:2])

TSAopt <- acast(X2[X2$qseqid %in% optimal_TSA,],Patient.ID~qseqid,value.var="n")
TSAopt <- rbind(TSAopt,matrix(NA,sum(unique(X2$Patient.ID) %notin% row.names(TSAopt)),dim(TSAopt)[2],dimnames=list(unique(X2$Patient.ID)[unique(X2$Patient.ID) %notin% row.names(TSAopt)],colnames(TSAopt))))          
TSAopt[is.na(TSAopt)] <- 0
