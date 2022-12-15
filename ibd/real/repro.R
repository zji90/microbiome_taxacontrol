library(reshape2)
library(lmerTest)
library(parallel)

d <- readRDS('ibd/data/proc/proc.rds')
b <- read.table('ibd/data/dysbio/HMP2_metadata.tsv',as.is=T,sep='\t',header=T)
b <- b[b[,1]%in%colnames(d$rna),]
b <- b[!is.na(b$age),]
source('ibd/model/fit.R')

tb <- b[b$diagnosis=='CD',]
group <- tb$dysbiosis_binary=='Yes'
names(group) <- tb$ID
id <- names(group)
ind <- b$subject[match(id,b$ID)]
age <- b$age[match(id,b$ID)]
ant <- b$antibiotics[match(id,b$ID)]=='Yes'
design <- data.frame(group,age,ant,subject=b$subject[match(id,b$ID)])

  us <- unique(design$subject)
  s1 <- sample(us,length(us)/2)
  s1 <- rownames(design)[design$subject %in% s1]
  s2 <- setdiff(rownames(design),s1)
  
  for (type in c('dnaabu','dna','abu')) {
    suppressMessages(suppressWarnings(k1 <- fit(d$rna[,s1],d$dna[,s1],d$abu[,s1],b$subject[match(s1,b$ID)],d$species,d$pathway,design[s1,],control=type)))
    suppressMessages(suppressWarnings(k2 <- fit(d$rna[,s2],d$dna[,s2],d$abu[,s2],b$subject[match(s2,b$ID)],d$species,d$pathway,design[s2,],control=type)))
    
    saveRDS(list(k1,k2),file=paste0('ibd/real/repro/res/',type,'.rds'))  
  }

