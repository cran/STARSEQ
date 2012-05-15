#pleio-map tests;
STAR.test <- function(dat.ped,y2.vec,covar.mat,rv.count.cutoff,maf.vec,maf.cutoff,pct.upper,pct.lower,yub,ylb,no.perm=1000,method='CMC',alternative='two.sided')
  {
    marker <- dat.ped$genotype;
    #maf.cutoff <- 0.05;
    #removing monomorphic sites;
    ind.mono <- which(colSums(marker)==0);

    if(length(ind.mono)>0)
      {
        marker <- as.matrix(marker[,-ind.mono]);
        maf.vec <- maf.vec[-ind.mono];
      }
    
    ind.rare <- which(maf.vec<maf.cutoff);
    marker <- as.matrix(marker[,ind.rare]);
    maf.vec <- maf.vec[ind.rare];

    Nub <- length(which(dat.ped$phenotype>yub));
    Nlb <- length(which(dat.ped$phenotype<ylb));
    
    y1.vec <- dat.ped$phenotype;
    x1.vec <- (rowSums(marker)>0);#other coding is also possible; But ANRV coding is used here;

    res.fit.null <- fit.null.2nd.ex.full.pleiomap(marker=marker,rv.count.cutoff=rv.count.cutoff,covar.mat=covar.mat,y1.vec=y1.vec,y2.vec=y2.vec,yub=yub,ylb=ylb,pct.upper=pct.upper,pct.lower=pct.lower,Nub=Nub,Nlb=Nlb);
    
    r2.vec <- res.fit.null$r.y2;
    
    par.dat <- list(yub.pct=pct.upper,ylb.pct=pct.lower);

    if(method=='CMC') rv.test <- cmc.uniqtl.pop;
    if(method=='WSS') rv.test <- wss.uniqtl.pop.C;
    if(method=='VT') rv.test <- vt.uniqtl.pop.C;
    if(method=='KBAC') rv.test <- kbac.uniqtl.pop.C;
    if(method=='SKAT') rv.test <- skat.uniqtl.simple.C;

    dat.ped.new <- list(genotype=marker,phenotype=r2.vec);
    
    res.test <- rv.test(dat.ped.new,par.dat,maf.vec,maf.cutoff,no.perm,alternative);

    res <- c(res.fit.null,res.test);
    #return(res.fit.null);
    return(res);
  }
 
