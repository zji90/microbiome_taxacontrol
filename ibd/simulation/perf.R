library(parallel)
perfun <- function(d) {
  fdr <- function(fdrv,gene) {
    diff <- ifelse(grepl(':TRUE',gene),'diff','nodiff')
    sumdifffval <- sum(diff=='diff')
    perf <- t(sapply(1:length(diff), function(i) {
      num <- sum(diff[1:i]=='diff')
      c(num/sumdifffval,(i-num)/i,fdrv[i])
    }))
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
  
  tar <- readRDS(paste0('ibd/simulation/rcv/tar.rds'))
  
  n <- paste0(rownames(d),':',rownames(d) %in% tar)
  #of <- fdrfunc(sort(d),as.vector(oraclenull))
  #of <- fdr(of,names(of))
  d$FDR <- p.adjust(d$pval,method='fdr')
  f <- fdr(d[,3],n)
  c(auc(f),fdrdiff(f),length(intersect(tar,rownames(d)[d$FDR < 0.05]))/length(tar))
}

af <- list.files('simulation/res/')
fr <- do.call(rbind,mclapply(af,function(f) {
  r <- readRDS(paste0('ibd/simulation/res/',f))
  k1 <- perfun(r[[1]])
  k2 <- perfun(r[[2]])
  k3 <- perfun(r[[3]])
  rr <- data.frame(rbind(k1,k2,k3))
  colnames(rr) <- c('auc','fdrdiff','power')
  rr$type <- c('abudna','dna','abu')
  rr$sig <- sub('.rds','',sub('.*_','',f))
  #rr$iter <- sub('_.*','',f)
  rr
},mc.cores=detectCores()))
fr$sig <- as.numeric(fr$sig)
saveRDS(fr,file='ibd/simulation/summary/perf.rds')


