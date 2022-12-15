library(reshape2)
library(lmerTest)
library(parallel)
library(ggpubr)
d <- readRDS('data/proc/proc.rds')
b <- read.table('data/dysbio/HMP2_metadata.tsv',as.is=T,sep='\t',header=T)
b <- b[b[,1]%in%colnames(d$rna),]
b <- b[!is.na(b$age),]
source('model/fit.R')

tb <- b[b$diagnosis=='CD',]
group <- tb$dysbiosis_binary=='Yes'
names(group) <- tb$ID
id <- names(group)
ind <- b$subject[match(id,b$ID)]
age <- b$age[match(id,b$ID)]
ant <- b$antibiotics[match(id,b$ID)]=='Yes'
design <- data.frame(group,age,ant)
rna=d$rna[,id];dna=d$dna[,id];abu=d$abu[,id];species=d$species;pathway=d$pathway;includeabu=T;abucut=0;samplecut=0;featurecut=0;tabcutoff=10;includeabu=T

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
df <- data.frame(species=species[rna[id,1]],pathway=pathway[rna[id,1]],sample=as.character(rna[id,2]),rna=rna[id,3],dna=dna[id,3],abu=abu[id,3],ind=ind[rna[id,2]],stringsAsFactors = F)

df$feature <- paste0(df$species,';',df$pathway)

tab <- do.call(rbind,tapply(design[df$sample,1],list(df$feature),function(i) {
  tab <- table(i)
  c(length(tab),min(tab))
}))
tar <- rownames(tab)[tab[,1]==2&tab[,2]>=tabcutoff]
df <- df[df$feature %in% tar,]

df <- data.frame(df,design[df$sample,])

sdf <- split(df,df$feature)

pval <- mclapply(sdf,function(ss) {
  r1 <- lmer(as.formula(paste0('rna ~ ',paste0(colnames(design),collapse = ' + '),' + dna + abu + (1|ind)')),data=ss,REML=F)
  r2 <- lm(as.formula(paste0('rna ~ ',paste0(colnames(design),collapse = ' + '),' + dna + abu')),data=ss)
  anova(r1,r2)["Pr(>Chisq)"][2,1]
},mc.cores = detectCores())
names(pval) <- names(sdf)
pval <- unlist(pval)
fdr_dnaabu <- p.adjust(pval,method='fdr')
names(fdr_dnaabu) <- names(pval)

pval <- mclapply(sdf,function(ss) {
  r1 <- lmer(as.formula(paste0('rna ~ ',paste0(colnames(design),collapse = ' + '),' + dna + (1|ind)')),data=ss,REML=F)
  r2 <- lm(as.formula(paste0('rna ~ ',paste0(colnames(design),collapse = ' + '),' + dna')),data=ss)
  anova(r1,r2)["Pr(>Chisq)"][2,1]
},mc.cores = detectCores())
names(pval) <- names(sdf)
pval <- unlist(pval)
fdr_dna <- p.adjust(pval,method='fdr')
names(fdr_dna) <- names(pval)

pval <- mclapply(sdf,function(ss) {
  r1 <- lmer(as.formula(paste0('rna ~ ',paste0(colnames(design),collapse = ' + '),' + abu + (1|ind)')),data=ss,REML=F)
  r2 <- lm(as.formula(paste0('rna ~ ',paste0(colnames(design),collapse = ' + '),' + abu')),data=ss)
  anova(r1,r2)["Pr(>Chisq)"][2,1]
},mc.cores = detectCores())
names(pval) <- names(sdf)
pval <- unlist(pval)
fdr_abu <- p.adjust(pval,method='fdr')
names(fdr_abu) <- names(pval)

df1 <- readRDS('real/res/CD_dnaabu.rds')
df2 <- readRDS('real/res/CD_dna.rds')
df3 <- readRDS('real/res/CD_abu.rds')
g1 <- rownames(df1)[df1[,2] < 0.05]
g2 <- rownames(df2)[df2[,2] < 0.05]
g3 <- rownames(df3)[df3[,2] < 0.05]

library(ggplot2)
for (n in setdiff(g1,g2)) {
  pval <- signif(c(df1[n,2],df2[n,2]),digits=2)
  ss <- sdf[[n]]
  if (fdr_dna[n] < 0.05) {
    r <- lmer(as.formula(paste0('rna ~ ',paste0(setdiff(colnames(design),'group'),collapse = ' + '),' + dna + (1|ind)')),data=ss)
  } else {
    r <- lm(as.formula(paste0('rna ~ ',paste0(setdiff(colnames(design),'group'),collapse = ' + '),' + dna')),data=ss)
  }
  resid <- resid(r)
  d1 <- data.frame(residual=resid,group=ss$group,type='DNA')
  
  if (fdr_dnaabu[n] < 0.05) {
    r <- lmer(as.formula(paste0('rna ~ abu+',paste0(setdiff(colnames(design),'group'),collapse = ' + '),' + dna + (1|ind)')),data=ss)
  } else {
    r <- lm(as.formula(paste0('rna ~ abu+',paste0(setdiff(colnames(design),'group'),collapse = ' + '),' + dna')),data=ss)
  }
  resid <- resid(r)
  d2 <- data.frame(residual=resid,group=ss$group,type='DNA+Taxa')
  
  df <- rbind(d1,d2)
  df$type <- factor(df$type,levels=c('DNA+Taxa','DNA'))
  df$group <- factor(ifelse(df$group,'dysbiotic','non-dysbiotic'),levels=c('dysbiotic','non-dysbiotic'))
  stat <- data.frame(group1='dysbiotic',group2='non-dysbiotic',p=pval,y.position=max(df$residual),type=c('DNA+Taxa','DNA'))
  pdf(paste0('real/plot/single/dnaabu_dna/',n,'.pdf'),width=3,height=3)
  print(ggplot() + geom_boxplot(data=df,aes(x=type,y=residual,fill=group)) + theme_classic() + xlab('') + ylab('Residual') + theme(legend.title = element_blank()) + scale_fill_manual(values=c('dysbiotic'='orange','non-dysbiotic'='royalblue')) + theme(legend.position = 'none') + stat_pvalue_manual(data=stat,x='type'))
  dev.off()
}


for (n in setdiff(g1,g3)) {
  pval <- signif(c(df1[n,2],df3[n,2]),digits=2)
  ss <- sdf[[n]]
  if (fdr_abu[n] < 0.05) {
    r <- lmer(as.formula(paste0('rna ~ ',paste0(setdiff(colnames(design),'group'),collapse = ' + '),' + abu + (1|ind)')),data=ss)
  } else {
    r <- lm(as.formula(paste0('rna ~ ',paste0(setdiff(colnames(design),'group'),collapse = ' + '),' + abu')),data=ss)
  }
  resid <- resid(r)
  d1 <- data.frame(residual=resid,group=ss$group,type='Taxa')
  
  if (fdr_dnaabu[n] < 0.05) {
    r <- lmer(as.formula(paste0('rna ~ abu+',paste0(setdiff(colnames(design),'group'),collapse = ' + '),' + dna + (1|ind)')),data=ss)
  } else {
    r <- lm(as.formula(paste0('rna ~ abu+',paste0(setdiff(colnames(design),'group'),collapse = ' + '),' + dna')),data=ss)
  }
  resid <- resid(r)
  d2 <- data.frame(residual=resid,group=ss$group,type='DNA+Taxa')
  
  df <- rbind(d1,d2)
  df$type <- factor(df$type,levels=c('DNA+Taxa','Taxa'))
  df$group <- factor(ifelse(df$group,'dysbiotic','non-dysbiotic'),levels=c('dysbiotic','non-dysbiotic'))
  stat <- data.frame(group1='dysbiotic',group2='non-dysbiotic',p=pval,y.position=max(df$residual),type=c('DNA+Taxa','Taxa'))
  
  pdf(paste0('real/plot/single/dnaabu_abu/',sub('/','_',n),'.pdf'),width=3,height=3)
  print(ggplot() + geom_boxplot(data=df,aes(x=type,y=residual,fill=group)) + theme_classic() + xlab('') + ylab('Residual') + theme(legend.title = element_blank()) + scale_fill_manual(values=c('dysbiotic'='orange','non-dysbiotic'='royalblue')) + theme(legend.position = 'none') + stat_pvalue_manual(data=stat,x='type'))
  dev.off()
}

pdf(paste0('real/plot/single/legend.pdf'),width=3,height=3)
print(ggplot(df,aes(x=type,y=residual,fill=group)) + geom_boxplot() + theme_classic() + xlab('') + ylab('Residual') + theme(legend.title = element_blank()) + scale_fill_manual(values=c('dysbiotic'='orange','non-dysbiotic'='royalblue')) + theme(legend.position = 'bottom'))
dev.off()


