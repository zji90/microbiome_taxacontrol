library(reshape2)
library(lmerTest)
library(parallel)
library(ggplot2)
corfunc <- function(abutrans,dnarnatrans) {
  d <- readRDS('ibd/data/proc/proc.rds')
  
  rna=d$rna
  dna=d$dna
  abu=d$abu
  species=d$species
  pathway=d$pathway
  abucut=0;samplecut=0;featurecut=0;tabcutoff=10
  
  rawabu <- abu
  abu[abu==0] <- NA
  #CLR logit raw
  if (abutrans == 'CLR') {
    abu <- log(abu)
    abu <- sweep(abu,2,colMeans(abu,na.rm=T),'-')  
  } else if (abutrans == 'logit') {
    abu <- log(abu/(1-abu))
  }
  
  #log2CPM CPM
  if (dnarnatrans == 'CPM') {
    rna <- 2^rna-1
    dna <- 2^dna-1
  }
  
  rownames(rna) <- rownames(dna) <- 1:nrow(rna)
  rna <- melt(rna)
  dna <- melt(dna)
  abu <- melt(abu[species,])
  rawabu <- melt(rawabu[species,])
  id <- which(rawabu[,3] > abucut & dna[,3] > 0 & rna[,3] > 0)
  df <- data.frame(species=species[rna[id,1]],pathway=pathway[rna[id,1]],sample=as.character(rna[id,2]),rna=rna[id,3],dna=dna[id,3],abu=abu[id,3],stringsAsFactors = F)
  
  df$feature <- paste0(df$species,';',df$pathway)
  
  tab2 <- tapply(df$rna,list(df$feature),function(i) {sum(i > 0)})
  df <- df[df$feature %in% names(tab2)[tab2 >= tabcutoff],]
  
  sdf <- split(df,df$feature)
  
  cor <- mclapply(names(sdf),function(n) {
    ss <- sdf[[n]]
    r1 <- lm(rna~dna,data=ss)
    r2 <- lm(abu~dna,data=ss)
    cor(resid(r1),resid(r2))
  },mc.cores=detectCores())
  names(cor) <- names(sdf)
  rnaabu_dna_cor <- unlist(cor)
  
  cor <- mclapply(names(sdf),function(n) {
    ss <- sdf[[n]]
    r1 <- lm(rna~abu,data=ss)
    r2 <- lm(dna~abu,data=ss)
    cor(resid(r1),resid(r2))
  },mc.cores=detectCores())
  names(cor) <- names(sdf)
  rnadna_abu_cor <- unlist(cor)
  list(rnaabu_dna_cor=rnaabu_dna_cor,rnadna_abu_cor=rnadna_abu_cor)
}

ar <- list()
for (abutrans in c('CLR','logit','raw')) {
  for (dnarnatrans in c('CPM','log2CPM')) {
    print(dnarnatrans)
    ar[[paste0(abutrans,'_',dnarnatrans)]] <- corfunc(abutrans,dnarnatrans)
  }
}
library(reshape2)
library(pheatmap)
library(gplots)
pd <- sapply(ar,function(i) i[[1]])
m <- colnames(pd)
m[m=='CLR_CPM'] <- 'CLR(Taxa)\nCPM(RNA,DNA)'
m[m=='CLR_log2CPM'] <- 'CLR(Taxa)\nlog2CPM(RNA,DNA)'
m[m=='logit_CPM'] <- 'logit(Taxa)\nCPM(RNA,DNA)'
m[m=='logit_log2CPM'] <- 'logit(Taxa)\nlog2CPM(RNA,DNA)'
m[m=='raw_CPM'] <- 'Raw(Taxa)\nCPM(RNA,DNA)'
m[m=='raw_log2CPM'] <- 'Raw(Taxa)\nlog2CPM(RNA,DNA)'
colnames(pd) <- m
dn <- dimnames(pd)
pd <- matrix(cut(pd,c(-1,-0.6,-0.3,0,0.3,0.6,1)),nrow(pd))
dimnames(pd) <- dn
pd <- melt(pd)
ln <- rev(c("(-1,-0.6]","(-0.6,-0.3]","(-0.3,0]",'(0,0.3]','(0.3,0.6]','(0.6,1]'))
pd$value <- factor(pd$value,levels=ln)
colnames(pd)[3] <- 'correlation'
pdf('ibd/real/plot/cor/bar_rnaabu_dna_cor.pdf',width=5,height=3.5)
vll <- rep(c('red','royalblue'),each=3)
names(vll) <- ln
all <- c(1,0.6,0.3,0.3,0.6,1)
names(all) <- ln
ggplot(pd,aes(x=Var2,alpha=correlation,fill=correlation)) + geom_bar() + theme_classic() + scale_fill_manual(values=vll) + scale_alpha_manual(values=all) + coord_flip() + xlab('') + ylab('Number of features')
dev.off()


pd <- sapply(ar,function(i) i[[2]])
m <- colnames(pd)
m[m=='CLR_CPM'] <- 'CLR(Taxa)\nCPM(RNA,DNA)'
m[m=='CLR_log2CPM'] <- 'CLR(Taxa)\nlog2CPM(RNA,DNA)'
m[m=='logit_CPM'] <- 'logit(Taxa)\nCPM(RNA,DNA)'
m[m=='logit_log2CPM'] <- 'logit(Taxa)\nlog2CPM(RNA,DNA)'
m[m=='raw_CPM'] <- 'Raw(Taxa)\nCPM(RNA,DNA)'
m[m=='raw_log2CPM'] <- 'Raw(Taxa)\nlog2CPM(RNA,DNA)'
colnames(pd) <- m
dn <- dimnames(pd)
pd <- matrix(cut(pd,c(-1,-0.6,-0.3,0,0.3,0.6,1)),nrow(pd))
dimnames(pd) <- dn
pd <- melt(pd)
ln <- rev(c("(-1,-0.6]","(-0.6,-0.3]","(-0.3,0]",'(0,0.3]','(0.3,0.6]','(0.6,1]'))
pd$value <- factor(pd$value,levels=ln)
colnames(pd)[3] <- 'correlation'
pdf('ibd/real/plot/cor/bar_rnadna_abu_cor.pdf',width=5,height=3.5)
vll <- rep(c('red','royalblue'),each=3)
names(vll) <- ln
all <- c(1,0.6,0.3,0.3,0.6,1)
names(all) <- ln
ggplot(pd,aes(x=Var2,alpha=correlation,fill=correlation)) + geom_bar() + theme_classic() + scale_fill_manual(values=vll) + scale_alpha_manual(values=all) + coord_flip() + xlab('') + ylab('Number of features')
dev.off()

