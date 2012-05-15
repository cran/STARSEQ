#source('datagen.uniqtl.R');source('wss.uniqtl.R');source('score.stat.R');pars <- list(maf=rep(0.005,10),ind.noncausal=1:5,yub=qnorm(0.7),ylb=qnorm(0.3),Nub=300,Nlb=300,beta1=0.25,scenario='beta1',maf.cutoff=0.1);res.datagen <- datagen.uniqtl(pars);dat.ped <- res.datagen$dat.sim;par.dat=res.datagen$par.dat;maf.vec <- pars$maf;maf.cutoff <- pars$maf.cutoff;res.wss <- wss.uniqtl(dat.ped,par.dat,maf.vec,maf.cutoff,1000,alternative='twosided');print(res.wss);
wss.uniqtl.pop.C <- function(dat.ped,par.dat,maf.vec,maf.cutoff,no.perm=1000,alternative='greater')
{
  yub.pct <- par.dat$yub.pct;
  ylb.pct <- par.dat$ylb.pct;
  if(length(yub.pct)==0) yub.pct <- 0.75;
  if(length(ylb.pct)==0) ylb.pct <- 0.25;
  
  
  y <- dat.ped$phenotype;
  yub <- as.numeric(quantile(y,yub.pct));
  ylb <- as.numeric(quantile(y,ylb.pct));
  marker <- dat.ped$genotype;

  ind.rare <- which(maf.vec<maf.cutoff);
  marker <- marker[,ind.rare];
  
  res.C <- .C("genericPerm2",markerPt=as.double(marker),
              yPt=as.double(y),
              yUbPt=as.double(yub),
              yLbPt=as.double(ylb),
              noPersonPt=as.integer(length(y)),
              noMarkerPt=as.integer(ncol(marker)),
              noPermPt=as.integer(no.perm),
              pValuePt=as.double(0),
              testPt=as.integer(3),
              markerMaxPt=as.double(2),
              statisticPt=as.double(rep(0,no.perm+1)));
  res <- list(p.value=res.C$pValuePt,
              statistic=res.C$statisticPt[1],
              statperm=res.C$statisticPt[-1]);  
  return(res);
  
}
