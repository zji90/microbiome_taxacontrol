d1 <- readRDS('real/res/CD_dnaabu.rds')
d2 <- readRDS('real/res/CD_dna.rds')
d3 <- readRDS('real/res/CD_abu.rds')
g1 <- rownames(d1)[d1[,2] < 0.05]
g2 <- rownames(d2)[d2[,2] < 0.05]
g3 <- rownames(d3)[d3[,2] < 0.05]

int <- setdiff(unique(c(g1,g2,g3)),intersect(intersect(g1,g2),g3))
d <- cbind(d1[int,2],d2[int,2],d3[int,2])
rownames(d) <- int
colnames(d) <- c('DNA+Taxa','DNA','Taxa')
path <- sub('.*s__','',rownames(d)[order(d[,1]<0.05,d[,2]<0.05,d[,3]<0.05)])
library(reshape2)
pd <- melt(d)
pd$path <- sub('.*;','',pd$Var1)
pd$species <- sub('.*s__','',sub(';.*','',pd$Var1))
ann <- read.csv('real/ann/ann.csv',as.is=T)

lite <- rep('No pubmed\nsupport',nrow(pd))
lite[pd$path %in% ann[ann$pubmed_number<=30&ann$pubmed_number>0,'pathway']] <- '1-30 pubmed\nsupport'
lite[pd$path %in% ann[ann$pubmed_number>30,'pathway']] <- '> 30 pubmed\nsupport'

pd$lite <- factor(lite,levels=c('> 30 pubmed\nsupport','1-30 pubmed\nsupport','No pubmed\nsupport'))
pd$Var1 <- factor(sub('.*s__','',pd$Var1),levels=path)
pd$sig <- pd$value < 0.05
pd$fdr <- -log10(pd$value)
pd$sig <- factor(as.character(pd$sig),levels=c('TRUE','FALSE'))

library(egg)
library(ggplot2)

g <- ggplot() + geom_point(data=pd,aes(x=Var2,y=Var1,color=sig,size=fdr)) + theme_classic() + scale_color_manual(values=c('TRUE'='red','FALSE'='lightpink')) + xlab('') + ylab('') + guides(color=guide_legend(title="Significance")) + guides(size=guide_legend(title="-log10(FDR)")) + facet_grid(lite~.,scale='free_y',space = "free_y") + theme(strip.background = element_blank(),strip.text.y=element_text(angle=0,size=12),legend.position = 'bottom',legend.box = 'vertical',legend.margin=margin(-2,-2,-2,-2))

library(ggplot2)
pdf('real/plot/heatmap.pdf',height=8,width=9.5)
#grid::grid.draw(gt)
g
dev.off()




