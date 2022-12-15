library(reshape2)
library(lmerTest)
library(parallel)

d <- readRDS('ibd/data/proc/proc.rds')

rna=d$rna
dna=d$dna
abu=d$abu
species=d$species
pathway=d$pathway
abucut=0;samplecut=0;featurecut=0;tabcutoff=10

rawabu <- abu
abu[abu==0] <- NA
abu <- log(abu)
abu <- sweep(abu,2,colMeans(abu,na.rm=T),'-')
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
cor <- unlist(cor)
saveRDS(cor,file='ibd/real/plot/cor/rnaabu_dna_cor.rds')

library(ggplot2)
pdf('ibd/real/plot/cor/rnaabu_dna_hist.pdf',width=3,height=3)
ggplot(data.frame(cor),aes(cor)) + geom_histogram(color="grey30",alpha=0.4,fill="royalblue",size=0.3) + theme_classic() + xlab('Correlation') + ylab('Frequency') + theme(legend.position = 'none')
dev.off()

highn <- names(tab2)[tab2 >= 30]
pn <- names(sort(cor[highn],decreasing = T)[1])
ss <- sdf[[pn]]
r1 <- lm(rna~dna,data=ss)
r2 <- lm(abu~dna,data=ss)

pdf(paste0('ibd/real/plot/cor/rnaabu_dna_',pn,'.pdf'),width=3,height=3)
ggplot(data.frame(x=resid(r1),y=resid(r2)),aes(x=x,y=y)) + geom_point(color='orange') + geom_smooth(method='lm',color='royalblue',se=F) + theme_classic() + xlab('RNA | DNA') + ylab('Taxa | DNA') + theme(legend.position = 'none') + ggtitle(paste0('Correlation:',round(cor(resid(r1),resid(r2)),3)))
dev.off()

pn <- names(sort(cor[highn],decreasing = T)[30])
ss <- sdf[[pn]]
r1 <- lm(rna~dna,data=ss)
r2 <- lm(abu~dna,data=ss)

pdf(paste0('ibd/real/plot/cor/rnaabu_dna_',pn,'.pdf'),width=3,height=3)
ggplot(data.frame(x=resid(r1),y=resid(r2)),aes(x=x,y=y)) + geom_point(color='orange') + geom_smooth(method='lm',color='royalblue',se=F) + theme_classic() + xlab('RNA | DNA') + ylab('Taxa | DNA') + theme(legend.position = 'none') + ggtitle(paste0('Correlation:',round(cor(resid(r1),resid(r2)),3)))
dev.off()

pn <- names(sort(cor[highn],decreasing = F)[1])
ss <- sdf[[pn]]
r1 <- lm(rna~dna,data=ss)
r2 <- lm(abu~dna,data=ss)

pdf(paste0('ibd/real/plot/cor/rnaabu_dna_',pn,'.pdf'),width=3,height=3)
ggplot(data.frame(x=resid(r1),y=resid(r2)),aes(x=x,y=y)) + geom_point(color='orange') + geom_smooth(method='lm',color='royalblue',se=F) + theme_classic() + xlab('RNA | DNA') + ylab('Taxa | DNA') + theme(legend.position = 'none') + ggtitle(paste0('Correlation:',round(cor(resid(r1),resid(r2)),3)))
dev.off()

pn <- names(sort(cor[highn],decreasing = F)[30])
ss <- sdf[[pn]]
r1 <- lm(rna~dna,data=ss)
r2 <- lm(abu~dna,data=ss)

pdf(paste0('ibd/real/plot/cor/rnaabu_dna_',pn,'.pdf'),width=3,height=3)
ggplot(data.frame(x=resid(r1),y=resid(r2)),aes(x=x,y=y)) + geom_point(color='orange') + geom_smooth(method='lm',color='royalblue',se=F) + theme_classic() + xlab('RNA | DNA') + ylab('Taxa | DNA') + theme(legend.position = 'none') + ggtitle(paste0('Correlation:',round(cor(resid(r1),resid(r2)),3)))
dev.off()


