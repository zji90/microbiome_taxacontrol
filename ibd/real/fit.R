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
design <- data.frame(group,age,ant)
suppressMessages(suppressWarnings(k1 <- fit(d$rna[,id],d$dna[,id],d$abu[,id],ind,d$species,d$pathway,design,control='dnaabu')))
suppressMessages(suppressWarnings(k2 <- fit(d$rna[,id],d$dna[,id],d$abu[,id],ind,d$species,d$pathway,design,control='dna')))
suppressMessages(suppressWarnings(k3 <- fit(d$rna[,id],d$dna[,id],d$abu[,id],ind,d$species,d$pathway,design,control='abu')))
suppressMessages(suppressWarnings(k4 <- fit(d$rna[,id],d$dna[,id],d$abu[,id],ind,d$species,d$pathway,design,control='no')))
saveRDS(k1,file='ibd/real/res/CD_dnaabu.rds')
saveRDS(k2,file='ibd/real/res/CD_dna.rds')
saveRDS(k3,file='ibd/real/res/CD_abu.rds')
saveRDS(k4,file='ibd/real/res/CD_no.rds')

