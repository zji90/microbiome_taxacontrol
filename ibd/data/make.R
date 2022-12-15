library(data.table)
meta <- fread('data/metadata/hmp2_metadata.csv',data.table = F)
d <- fread('data/metagenome/pathabundances_3.tsv',data.table = F)
rownames(d) <- d[,1]
d <- as.matrix(d[,-1])
d <- d[!grepl('^UNMAPPED|^UNINTEGRATED|unclassified',rownames(d)),]
d <- d[,colSums(d) > 0]
colnames(d) <- sub('_pathabundance_cpm','',colnames(d))
d <- d[grep('\\|',rownames(d)),]
gend <- log2(d+1)

rd <- sub('_CAG_[0-9]*$','',rownames(gend))
gend <- t(sapply(unique(rd),function(i) {
  colMeans(gend[rd==i,,drop=F])
}))


d <- fread('data/metatranscriptome/pathabundances_3.tsv',data.table = F)
rownames(d) <- d[,1]
d <- as.matrix(d[,-1])
d <- d[!grepl('^UNMAPPED|^UNINTEGRATED|unclassified',rownames(d)),]
d <- d[,colSums(d) > 0]
colnames(d) <- sub('_pathabundance_cpm','',colnames(d))
d <- d[grep('\\|',rownames(d)),]
trad <- log2(d+1)

rd <- sub('_CAG_[0-9]*$','',rownames(trad))
trad <- t(sapply(unique(rd),function(i) {
  colMeans(trad[rd==i,,drop=F])
}))

rm('d')
int <- list(intersect(rownames(gend),rownames(trad)),intersect(colnames(gend),colnames(trad)))

gend <- gend[int[[1]],int[[2]]]
trad <- trad[int[[1]],int[[2]]]

p <- fread('data/metagenome/taxonomic_profiles_3.tsv',data.table = F)
rownames(p) <- p[,1]
p <- as.matrix(p[,-1])
colnames(p) <- sub('_profile','',colnames(p))
p <- p[!grepl('^UNKNOWN',rownames(p)),]
p <- p[,colSums(p) > 0]

p <- p[,colnames(trad)]
p <- p[grep('k.*\\|p.*\\|c.*\\|o.*\\|f.*\\|g.*\\|s.*',rownames(p)),]
rd <- sub('_CAG_[0-9]*$','',rownames(p))
p <- t(sapply(unique(rd),function(i) {
  colSums(p[rd==i,,drop=F])
}))
p <- p/100

species <- sub('.*\\|','',rownames(trad))
species <- sub('\\.','|',species)
path <- sub('\\|.*','',rownames(trad))
rownames(p) <- sub('.*\\|g__','g__',rownames(p))

meta <- meta[match(colnames(gend),meta[,2]),]

res <- list(dna=gend,rna=trad,abu=p,species=species,pathway=path,meta=meta)
saveRDS(res,file='data/proc/proc.rds')

