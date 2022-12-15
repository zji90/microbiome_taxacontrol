library(reshape2)
d <- readRDS('data/proc/proc.rds')
b <- read.table('data/dysbio/HMP2_metadata.tsv',as.is=T,sep='\t',header=T)
b <- b[b[,1]%in%colnames(d$rna),]
b <- b[!is.na(b$age),]

group <- rbinom(1:nrow(b),size = 1,prob=0.5)
names(group) <- b[,1]
id <- names(group)
ind <- b$subject[match(id,b$ID)]
age <- b$age[match(id,b$ID)]
ant <- b$antibiotics[match(id,b$ID)]=='Yes'
design <- data.frame(group,age,ant)

dna <- d$dna[,id]
rna <- d$rna[,id]
abu <- d$abu[,id]
species <- d$species
pathway <- d$pathway

abucut=0
samplecut=0
featurecut=0
tabcutoff=10

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

rcv1 <- sapply(sdf,function(i) {
  cor(resid(lm(i$rna~i$dna)),resid(lm(i$abu~i$dna)))
})
saveRDS(rcv1,file=paste0('simulation/rcv/rcv_dna.rds'))
rcv2 <- sapply(sdf,function(i) {
  cor(resid(lm(i$rna~i$abu)),resid(lm(i$dna~i$abu)))
})
saveRDS(rcv2,file=paste0('simulation/rcv/rcv_abu.rds'))

tar <- names(which(abs(rcv1) > quantile(abs(rcv1),0.95) | abs(rcv2) > quantile(abs(rcv2),0.95)))
saveRDS(tar,file=paste0('simulation/rcv/tar.rds'))

mp <- df$rna
rs <- cut(mp,quantile(mp,seq(0,1,0.1)),include.lowest=T)
rs <- as.numeric(rs)

odf <- df

for (factor in 1:10) {
    df <- odf
    for (star in tar) {
      bin <- rbinom(n = 1,size=1,prob=0.5)
      tmp <- df[df$feature==star & df$group==bin,'rna']
      df[df$feature==star & df$group==bin,'rna'] <- tmp + sample(mp[rs==factor],length(tmp))
    }
    saveRDS(df,file=paste0('simulation/data/',factor,'.rds'))
}





