library(parallel)
r <- sapply(c('dnaabu','dna','abu'),function(type) {
  d <- readRDS(paste0('real/repro/res/',type,'.rds'))
  d1 <- d[[1]]
  d2 <- d[[2]]
  int <- intersect(rownames(d1),rownames(d2))
  d1 <- d1[int,]
  d2 <- d2[int,]
  d1 <- d1[order(d1$pval),]
  d2 <- d2[order(d2$pval),]
  real <- sapply(1:length(int),function(i) {
    length(intersect(rownames(d1)[1:i],rownames(d2)[1:i]))
  })
  real
})

library(ggplot2)
pd=data.frame(p=c(r[,1],r[,2],r[,3]),t=factor(rep(c('DNA+Taxa','DNA','Taxa'),each=nrow(r)),levels=c('DNA+Taxa','DNA','Taxa')),x=rep(1:nrow(r),3))
pd <- pd[pd$x <= 20,]
pd <- pd[nrow(pd):1,]
pdf('real/repro/summary/over.pdf',width=3,height=3.5)
ggplot(data=pd,aes(x=x,y=p,col=t)) + geom_point() + geom_line() + theme_classic() + xlab('Number of features') + ylab('Number of overlaps') + theme(legend.title = element_blank()) +  scale_color_manual(values=readRDS('pal/pal.rds')) + theme(legend.position='bottom')+ theme(legend.position='bottom')
dev.off()

