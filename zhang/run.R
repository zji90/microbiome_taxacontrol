library(data.table)
library(reshape2)

af <- unique(sub('\\.mtx_abunds.tsv','',list.files('zhang/data/simu/',pattern='mtx_abunds.tsv')))
rid <- as.numeric(commandArgs(trailingOnly = T))

sf <- af[rid]
rna <- fread(paste0('zhang/data/simu/',sf,'.mtx_abunds.tsv'),data.table=F)
rownames(rna) <- rna[,1]
rna <- as.matrix(rna[,-1])
group <- rna['Phenotype',]
rna <- rna[!rownames(rna)%in%c('Phenotype','SeqDepth'),]
rna <- t(log2(t(rna)/colSums(rna)*1e6+1))

dna <- fread(paste0('zhang/data/simu/',sf,'.mgx_abunds.tsv'),data.table=F)
rownames(dna) <- dna[,1]
dna <- as.matrix(dna[,-1])
dna <- dna[!rownames(dna)%in%c('Phenotype','SeqDepth'),]
dna <- t(log2(t(dna)/colSums(dna)*1e6+1))

abu <- fread(paste0('zhang/data/simu/',sf,'.bug_abunds.tsv'),data.table=F)
rownames(abu) <- abu[,1]
abu <- as.matrix(abu[,-1])
abu <- abu[!rownames(abu)%in%c('Phenotype','SeqDepth'),]
rawabu <- abu
abu[abu==0] <- NA
abu <- log(abu)
abu <- sweep(abu,2,colMeans(abu,na.rm=T),'-')

species <- sub('_.*','',rownames(dna))
pathway <- sub('.*_','',rownames(dna))
rna <- melt(rna)
dna <- melt(dna)
abu <- melt(abu[species,])
rawabu <- melt(rawabu[species,])

id <- which(rawabu[,3] > 0 & dna[,3] > 0 & rna[,3] > 0)
df <- data.frame(species=species[rna[id,1]],pathway=pathway[rna[id,1]],sample=as.character(rna[id,2]),rna=rna[id,3],dna=dna[id,3],abu=abu[id,3],stringsAsFactors = F)

df$feature <- paste0(df$species,'_',df$pathway)
tab <- do.call(rbind,tapply(group[df$sample],list(df$feature),function(i) {
  tab <- table(i)
  c(length(tab),min(tab))
}))
tabcutoff=10
tar <- rownames(tab)[tab[,1]==2&tab[,2]>=tabcutoff]
df <- df[df$feature %in% tar,]

df <- data.frame(df,group=group[df$sample])

sdf <- split(df,df$feature)

pval <- sapply(c('+ dna + abu','+dna','+abu'),function(incw) {
  sapply(names(sdf),function(n) {
    r <- lm(as.formula(paste0('rna ~ group',incw)),data=sdf[[n]])
    summary(r)$coefficients[2,c('Pr(>|t|)')]
  })
})

for (i in 1:ncol(pval))
  pval[,i] <- p.adjust(pval[,i],method='fdr')

saveRDS(pval,file=paste0('zhang/res/',sf,'.rds'))


