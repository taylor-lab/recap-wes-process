#!/usr/bin/env Rscript

##########################################################################################
##########################################################################################
# MSKCC CMO
# Test for evidence of WGD in FACETS outputs 
##########################################################################################
##########################################################################################

'%!in%' <- function(x,y)!('%in%'(x,y))

catverbose <- function(...){
  cat(format(Sys.time(), "%Y%m%d %H:%M:%S |"), ..., "\n")
}

getSDIR <- function(){
  args=commandArgs(trailing=F)
  TAG="--file="
  path_idx=grep(TAG,args)
  SDIR=dirname(substr(args[path_idx],nchar(TAG)+1,nchar(args[path_idx])))
  if(length(SDIR)==0) {
    return(getwd())
  } else {
    return(SDIR)
  }
}

wgd_test <- function(facets,out){
 out_file=as.character(out)
  cat(paste(out,"inside function!\n")) 
  summary <- c()
  flagged <- c()
  samples <- unique(facets$Sample_ID)

  for (s in samples){
    
    cat('\n')
    
    #catverbose(s)
    s.facets <- facets[Sample_ID == s]
    
    catverbose("Loading FACETS Rdata...",s)
    if(file.exists(s.facets$path)) {
      load(s.facets$path)
    } else { next }
    
    facets.fit <- as.data.table(fit$cncf)
    facets.fit[, mcn := tcn.em - lcn.em]
    flags <- out$flags
    wgd <- F
    alt.fit <- F
    
    dipLogR <- out$dipLogR
    if(abs(dipLogR) > 1){
      alt.fit <- T
    }
    
    purity <- fit$purity
    ploidy <- fit$ploidy
    n.bases <- sum(as.numeric(fit$seglen))
    
    if(!is.null(out$flags)){
      alt.fit <- T
      catverbose(flags)
    }
    
    if(n.bases > 0){
      loh <- sum(as.numeric(fit$seglen[which(facets.fit$lcn.em == 0)]))
      f_lo_lcn = loh / n.bases
      f_hi_mcn <- sum(as.numeric(fit$seglen[which((facets.fit$tcn.em - facets.fit$lcn.em) >= 2)])) / n.bases
    } else {
      loh <- f_lo_lcn <- f_hi_mcn <- NA
    }
    
    if(!is.na(f_hi_mcn) & f_hi_mcn > 0.5){ # Major copy number >= 2 across 50% of the genome
      wgd <- T
      alt.fit <- T
    }
    
    if(wgd){
    print(paste(wgd));
      facets.fit[, WGD := "WGD"]
    } else {
      facets.fit[, WGD := "no WGD"]
    }
      f_altered <- 0 # Fraction of genome altered
      f_altered_v2 <- 0 # Fraction of genome altered, excluding diploid regions in WGD cases
    if(wgd){
      cn.neutral <- sum(as.numeric(fit$seglen[which(facets.fit$tcn.em == 4)]))
      cn.neutral.v2 <- sum(as.numeric(fit$seglen[which((facets.fit$tcn.em == 4) | (facets.fit$tcn.em == 2) | (is.na(facets.fit$tcn.em)))]))
      f_altered <- 1 - (cn.neutral / n.bases)
      f_altered_v2 <- 1 - (cn.neutral.v2 / n.bases)
    } else {
      cn.neutral <- sum(as.numeric(fit$seglen[which(facets.fit$tcn.em == 2)]))
      f_altered <- f_altered_v2 <- 1 - (cn.neutral / n.bases)
    }
      #s.out <- c(s.facets$Sample_ID, s.facets$Cancer_Type, dipLogR, purity, ploidy, f_lo_lcn, f_hi_mcn, wgd, f_altered, f_altered_v2, s.facets$path)
     s.out<-c(s.facets$Sample_ID,dipLogR, purity, f_lo_lcn, f_hi_mcn, ploidy,wgd,f_altered,f_altered_v2)
      summary <- rbind(summary,s.out)
  }
  
      colnames(summary) <- c('Sample_ID','dipLogR','Purity','LOH','HI_MCN','Ploidy','WGD','FGA1','FGA2')
  #print(summary)
  #colnames(summary) <- c('Sample_ID', 'Cancer_Type', 'DipLogR', 'Purity', 'Ploidy', 'LOH', 'HI_MCN', 'WGD', 'FGA1', 'FGA2', 'Path')
  write.table(summary,file=out_file, append=F, quote=F, row.names=F, col.names=T, sep="\t")
  
}

  





