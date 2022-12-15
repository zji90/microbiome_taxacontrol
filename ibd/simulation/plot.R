library(ggplot2)
library(reshape2)
fr <- readRDS('ibd/simulation/summary/perf.rds')
pd <- melt(fr,id.vars=c('type','sig'))
colnames(pd) <- c('met','signal','metric','perf')
pd$met[pd$met=='dna'] <- 'DNA'
pd$met[pd$met=='abu'] <- 'Taxa'
pd$met[pd$met=='abudna'] <- 'DNA+Taxa'
pd$met <- factor(pd$met,levels=c('DNA+Taxa','DNA','Taxa'))
pd <- pd[pd$metric!='fdrdiff',]
pd$metric <- factor(c('Statistical power','AUC')[match(as.character(pd$metric),c('power','auc'))],levels=c('Statistical power','AUC'))

pd <- pd[order(match(pd[,1],c('Taxa','DNA','DNA+Taxa'))),]
pd$signal <- factor(pd$signal,levels=1:10)


pdf(paste0('ibd/simulation/summary/auc.pdf'),width=3,height=3)
ggplot() + geom_point(data=pd[pd[,3]=='AUC',],aes(x=signal,y=perf,col=met,group=met)) + geom_line(data=pd[pd[,3]=='AUC',],aes(x=signal,y=perf,col=met,group=met)) + theme_classic() + scale_color_manual(values=readRDS('ibd/pal/pal.rds')) + xlab('Signal Strength') + ylab('AUC') + theme(legend.title = element_blank(),legend.position = 'none')
dev.off()

spd <- pd[pd[,3]=='AUC',]
spd[spd[,1]=='Taxa',4] <- (spd[spd[,1]=='DNA+Taxa',4]-spd[spd[,1]=='Taxa',4])
spd[spd[,1]=='DNA',4] <- (spd[spd[,1]=='DNA+Taxa',4]-spd[spd[,1]=='DNA',4])
spd <- spd[spd[,1]=='DNA',]
spd[,1] <- 'AUC(DNA+Taxa) - AUC(DNA)'

pdf(paste0('ibd/simulation/summary/aucinc.pdf'),width=3,height=3)
ggplot() + geom_point(data=spd,aes(x=signal,y=perf,col=met,group=met)) + geom_line(data=spd,aes(x=signal,y=perf,col=met,group=met)) + theme_classic() + xlab('Signal Strength') + ylab('AUC difference') + theme(legend.margin=margin(t = -0.2, unit='cm'),legend.title = element_blank(),legend.position = 'bottom') + coord_cartesian(ylim=c(-0.1,0.1)) + scale_color_manual(values='black')

pdf(paste0('ibd/simulation/summary/power.pdf'),width=3,height=3)
ggplot() + geom_point(data=pd[pd[,3]=='Statistical power',],aes(x=signal,y=perf,col=met,group=met)) + geom_line(data=pd[pd[,3]=='Statistical power',],aes(x=signal,y=perf,col=met,group=met)) + theme_classic() + scale_color_manual(values=readRDS('ibd/pal/pal.rds')) + xlab('Signal Strength') + ylab('Statistical power') + theme(legend.title = element_blank(),legend.position = 'none')
dev.off()


spd <- pd[pd[,3]=='Statistical power',]
spd[spd[,1]=='Taxa',4] <- (spd[spd[,1]=='DNA+Taxa',4]-spd[spd[,1]=='Taxa',4])
spd[spd[,1]=='DNA',4] <- (spd[spd[,1]=='DNA+Taxa',4]-spd[spd[,1]=='DNA',4])
spd <- spd[spd[,1]=='DNA',]
spd[,1] <- 'power(DNA+Taxa) - power(DNA)'

pdf(paste0('ibd/simulation/summary/powerinc.pdf'),width=3,height=3)
print(ggplot() + geom_point(data=spd,aes(x=signal,y=perf,col=met,group=met)) + geom_line(data=spd[spd[,4]!=Inf,],aes(x=signal,y=perf,col=met,group=met)) + theme_classic() + xlab('Signal Strength') + ylab('Statistical power difference') + theme(legend.margin=margin(t = -0.2, unit='cm'),legend.title = element_blank(),legend.position = 'bottom') + coord_cartesian(ylim=c(-0.1,0.1)) + scale_color_manual(values='black'))
dev.off()

pdf(paste0('ibd/simulation/summary/legend.pdf'),width=3,height=3)
ggplot() + geom_point(data=pd[pd[,3]=='Statistical power',],aes(x=signal,y=perf,col=met,group=met)) + geom_line(data=pd[pd[,3]=='Statistical power',],aes(x=signal,y=perf,col=met,group=met)) + theme_classic() + scale_color_manual(values=readRDS('ibd/pal/pal.rds')) + xlab('Signal Strength') + ylab('AUC') + theme(legend.title = element_blank(),legend.position = 'bottom')
dev.off()


