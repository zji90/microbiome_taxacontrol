d <- readRDS('zhang/summary/signal.rds')

af <- names(d)
afr <- sapply(af,function(f) {
  r <- readRDS(paste0('zhang/res/',f))
  tar <- read.table(paste0('zhang/data/simu/',sub('.rds','',f),'.mtx_spiked.tsv'))[,1]
  tar <- intersect(tar,rownames(r))
  apply(r,2,function(d) {
    sum(d < 0.05&!rownames(r)%in%tar)/sum(d < 0.05)
  })
})
afr[is.na(afr)] <- 0

library(reshape2)
library(ggplot2)
tt <- sub('_.*','',colnames(afr))
for (st in unique(tt)) {
  p <- afr[,tt==st]
  rownames(p) <- c('DNA+Taxa','DNA','Taxa')
  colnames(p) <- sub('.*_','',sub('.rds','',colnames(p)))
  pd <- melt(p)
  pd[,2] <- as.numeric(sub(':.*','',pd[,2]))*10
  pd[,2] <- factor(pd[,2],levels=sort(unique(pd[,2])))
  pal <- readRDS('ibd/pal/pal.rds')
  pdf(paste0('zhang/fdrplot/',st,'.pdf'),width=5,height=2.5)
  print(ggplot(pd,aes(x=Var2,y=value,fill=Var1)) + geom_bar(stat='identity',position = 'dodge') + theme_classic() + geom_hline(yintercept = 0.05,linetype=2) + xlab('Signal Strength') + ylab('FDR') + theme(legend.position = 'right',legend.title = element_blank()) + scale_fill_manual(values=pal))
  dev.off()  
}

