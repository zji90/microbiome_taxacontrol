library(reshape2)
library(ggplot2)
library(RColorBrewer)

d <- read.csv('ibd/real/ann/ann.csv',as.is=T)
d$cat <- '0'
d$cat[d$pubmed_number <= 10&d$pubmed_number > 0] <- '1-10'
d$cat[d$pubmed_number > 10] <- '> 10'
d$cat <- factor(d$cat,levels=c('> 10','1-10','0'))

pd <- rbind(data.frame(pub=d$cat[d$dnaabu & !d$dna],type='DNA+Taxa differential\nDNA non-differential'),
            data.frame(pub=d$cat[d$dnaabu & !d$abu],type='DNA+Taxa differential\nTaxa non-differential'),
            data.frame(pub=d$cat[!d$dnaabu & d$dna],type='DNA differential\nDNA+Taxa non-differential'),
            data.frame(pub=d$cat[!d$dnaabu & d$abu],type='Taxa differential\nDNA+Taxa non-differential'))

tab <- table(pd[,2],pd[,1])
d <- melt(tab)
colnames(d) <- c('type','cat','count')
d$type <- factor(d$type,levels=rev(c('DNA+Taxa differential\nDNA non-differential','DNA differential\nDNA+Taxa non-differential','DNA+Taxa differential\nTaxa non-differential','Taxa differential\nDNA+Taxa non-differential')))

pdf('ibd/real/plot/pubmed.pdf',width=4.3,height=3.5)
ggplot(d,aes(x=type,fill=cat,y=count)) + geom_bar(stat='identity',width=0.7) + theme_classic() + scale_fill_manual(values=rev(c('grey',brewer.pal(3,'Oranges')[1:2]))) + xlab('') + ylab('Number of features') + guides(fill=guide_legend(title="# Pubmed\nsupport")) + theme(legend.position = 'bottom') + coord_flip()
dev.off()


