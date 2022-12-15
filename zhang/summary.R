library(parallel)
perfun <- function(d,tar) {
  fdr <- function(fdrv,gene) {
    diff <- ifelse(grepl(':TRUE',gene),'diff','nodiff')
    sumdifffval <- sum(diff=='diff')
    cs <- cumsum(diff=='diff')
    perf <- cbind(cs/sumdifffval,1-cs/c(1:length(diff)),fdrv)
    
    if (nrow(perf) > 1) {
      for (i in (nrow(perf)):2) {
        if (perf[i-1,2] > perf[i,2]) perf[i-1,2] = perf[i,2]
      }
    }
    colnames(perf) <- c('Sensitivity','Real_FDR','Reported_FDR')
    rbind(c(0,0,0),perf)
  }
  cutoff <- 0.1
  fdrdiff <- function(tmp) {
    bound <- approx(x=tmp[,3],y=tmp[,2],xout=cutoff)$y
    tmp <- rbind(tmp[tmp[,3] < cutoff,2:3],c(bound,cutoff))
    tmp <- unique(tmp)
    diff <- sum(sapply(2:nrow(tmp),function(i) (tmp[i-1,1]+tmp[i,1])*(tmp[i,2]-tmp[i-1,2])/2),na.rm = T)-cutoff*cutoff/2
  }
  
  auc <- function(tmp) {
    bound <- approx(x=tmp[,2],y=tmp[,1],xout=cutoff)$y
    tmp <- rbind(tmp[tmp[,2] < cutoff,1:2],c(bound,cutoff))
    area <- sum(sapply(2:nrow(tmp),function(i) (tmp[i-1,1]+tmp[i,1])*(tmp[i,2]-tmp[i-1,2])/2),na.rm=T)/cutoff
  }
  

  n <- paste0(names(d),':',names(d) %in% tar)
  f <- fdr(d,n)
  c(auc(f),fdrdiff(f),length(intersect(tar,names(d)[d < 0.05]))/length(tar))
}

af <- list.files('zhang/res/',pattern = 'true')
afr <- sapply(af,function(f) {
  r <- readRDS(paste0('zhang/res/',f))
  tar <- read.table(paste0('zhang/data/simu/',sub('.rds','',f),'.mtx_spiked.tsv'))[,1]
  tar <- intersect(tar,rownames(r))
  
  abudna <- perfun(sort(r[,1]),tar)
  dna <- perfun(sort(r[,2]),tar)
  abu <- perfun(sort(r[,3]),tar)
  rbind(abudna,dna,abu)
},simplify = F)

saveRDS(afr,file='zhang/summary/signal.rds')
