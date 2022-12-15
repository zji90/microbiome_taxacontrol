af <- list.files('simulation/res/')
tar <- readRDS(paste0('ibd/simulation/rcv/tar.rds'))

p <- sapply(af,function(f) {
  r <- readRDS(paste0('ibd/simulation/res/',f))
  sapply(r,function(d) {
    sum(d[,2] < 0.05&!rownames(d)%in%tar)/sum(d[,2] < 0.05)
  })
})

library(reshape2)
library(ggplot2)
rownames(p) <- c('DNA+Taxa','DNA','Taxa')
colnames(p) <- sub('.*_','',sub('.rds','',colnames(p)))
pd <- melt(p)
pd[,2] <- factor(pd[,2],levels=sort(unique(pd[,2])))
pd[is.na(pd[,3]),3] <- 0
pal <- readRDS('ibd/pal/pal.rds')
pdf('ibd/simulation/summary/fdr005.pdf',width=5,height=2)
ggplot(pd,aes(x=Var2,y=value,fill=Var1)) + geom_bar(stat='identity',position = 'dodge') + theme_classic() + geom_hline(yintercept = 0.05,linetype=2) + xlab('Signal Strength') + ylab('FDR') + theme(legend.position = 'right',legend.title = element_blank()) + scale_fill_manual(values=pal)
dev.off()
