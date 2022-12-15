library(parallel)
library(lmerTest)
af <- list.files('simulation/data/')

for (f in af) {
  if (!file.exists(paste0('simulation/res/',f))) {
    print(f)
    df <- readRDS(paste0('simulation/data/',f))
    sdf <- split(df,df$feature)
    
    incw <- '+ dna + abu'
    
    pval <- mclapply(sdf,function(ss) {
      r1 <- lmer(as.formula(paste0('rna ~ group + age + ant',incw,' + (1|ind)')),data=ss,REML=F)
      r2 <- lm(as.formula(paste0('rna ~ group + age + ant',incw)),data=ss)
      anova(r1,r2)["Pr(>Chisq)"][2,1]
    },mc.cores = detectCores())
    names(pval) <- names(sdf)
    pval <- unlist(pval)
    fdr <- p.adjust(pval,method='fdr')
    names(fdr) <- names(pval)
    
    pval <- mclapply(names(sdf),function(n) {
      ss <- sdf[[n]]
      if (fdr[n] > 0.05) {
        r <- lmer(as.formula(paste0('rna ~ group + age + ant',incw,' + (1|ind)')),data=ss)
      } else {
        r <- lm(as.formula(paste0('rna ~ group + age + ant',incw)),data=ss)
      }
      summary(r)$coefficients[2,'Pr(>|t|)']
    })
    
    names(pval) <- names(sdf)
    pval <- unlist(pval)
    res <- data.frame(pval=pval,fdr=p.adjust(pval,method='fdr'))
    k1 <- res[order(res[,1]),]
    
    
    
    incw <- ' + dna'
    
    pval <- mclapply(sdf,function(ss) {
      r1 <- lmer(as.formula(paste0('rna ~ group + age + ant',incw,' + (1|ind)')),data=ss,REML=F)
      r2 <- lm(as.formula(paste0('rna ~ group + age + ant',incw)),data=ss)
      anova(r1,r2)["Pr(>Chisq)"][2,1]
    },mc.cores = detectCores())
    names(pval) <- names(sdf)
    pval <- unlist(pval)
    fdr <- p.adjust(pval,method='fdr')
    names(fdr) <- names(pval)
    
    pval <- mclapply(names(sdf),function(n) {
      ss <- sdf[[n]]
      if (fdr[n] > 0.05) {
        r <- lmer(as.formula(paste0('rna ~ group + age + ant',incw,' + (1|ind)')),data=ss)
      } else {
        r <- lm(as.formula(paste0('rna ~ group + age + ant',incw)),data=ss)
      }
      summary(r)$coefficients[2,'Pr(>|t|)']
    })
    
    names(pval) <- names(sdf)
    pval <- unlist(pval)
    res <- data.frame(pval=pval,fdr=p.adjust(pval,method='fdr'))
    k2 <- res[order(res[,1]),]
    
    
    
    
    incw <- ' + abu'
    
    pval <- mclapply(sdf,function(ss) {
      r1 <- lmer(as.formula(paste0('rna ~ group + age + ant',incw,' + (1|ind)')),data=ss,REML=F)
      r2 <- lm(as.formula(paste0('rna ~ group + age + ant',incw)),data=ss)
      anova(r1,r2)["Pr(>Chisq)"][2,1]
    },mc.cores = detectCores())
    names(pval) <- names(sdf)
    pval <- unlist(pval)
    fdr <- p.adjust(pval,method='fdr')
    names(fdr) <- names(pval)
    
    pval <- mclapply(names(sdf),function(n) {
      ss <- sdf[[n]]
      if (fdr[n] > 0.05) {
        r <- lmer(as.formula(paste0('rna ~ group + age + ant',incw,' + (1|ind)')),data=ss)
      } else {
        r <- lm(as.formula(paste0('rna ~ group + age + ant',incw)),data=ss)
      }
      summary(r)$coefficients[2,'Pr(>|t|)']
    })
    
    names(pval) <- names(sdf)
    pval <- unlist(pval)
    res <- data.frame(pval=pval,fdr=p.adjust(pval,method='fdr'))
    k3 <- res[order(res[,1]),]
    
    saveRDS(list(dnaabu=k1,dna=k2,abu=k3),file=paste0('simulation/res/',f))
  }
}



