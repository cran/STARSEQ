#source('include.uniqtl.R');pars <- list(maf=rep(0.005,20),ind.noncausal=1:10,yub=qnorm(0.7),ylb=qnorm(0.3),Nub=300,Nlb=300,beta1=-0.5,scenario='beta1',maf.cutoff=0.1);res.datagen <- datagen.uniqtl(pars);dat.ped <- res.datagen$dat.sim;par.dat=res.datagen$par.dat;maf.vec <- pars$maf;maf.cutoff <- pars$maf.cutoff;res.vt <- vt.uniqtl(dat.ped,par.dat,maf.vec,maf.cutoff);print(res.vt);vt.proscore(dat.ped,par.dat,maf.vec,maf.cutoff);

vt.uniqtl.pop.C <- function(dat.ped,par.dat,maf.vec,maf.cutoff,no.perm=1000,alternative=c('two.sided','greater','less'))
{

  if(length(alternative)>1) alternative <- alternative[1];
  
  vt.proscore.simple <- function(dat.ped,par.dat,maf.vec,maf.cutoff,no.perm=1000,alternative='greater')
    {
      yub.pct <- par.dat$yub.pct;
      ylb.pct <- par.dat$ylb.pct;

      if(length(yub.pct)==0) yub.pct <- 0.75;
      if(length(ylb.pct)==0) ylb.pct <- 0.25;
      sampling.status <- as.numeric((dat.ped$phenotype)>0);
      ind.rare <- which(maf.vec<maf.cutoff);
      dat.mat <- cbind(dat.ped$genotype,dat.ped$phenotype);
      dat.binary <- cbind(dat.ped$genotype,sampling.status);
      ind.rare <- which(maf.vec<maf.cutoff);
      marker <- dat.mat[,-c(ncol(dat.mat))];
      marker <- marker[,ind.rare];
      y.dat <- dat.mat[,ncol(dat.mat)];
      ##yub <- quantile(y.dat,yub.pct);
      ##ylb <- quantile(y.dat,ylb.pct);
      y.binary <- sampling.status;
      dat.mat <- cbind(marker,y.dat);
      dat.binary <- cbind(marker,y.binary);
      ind.rm <- which(colSums(marker)==0);
      if(length(ind.rm)>0) {
        marker <- marker[,-ind.rm];
        dat.mat <- cbind(marker,y.dat);
      }
      marker <- as.matrix(marker);
      rv.count <- colSums(as.matrix(marker));
      
      if(length(rv.count)<1) return(1);
      rv.ix <- sort(unique(rv.count));
      
      if(length(rv.ix)<1) return(1);
      x.mat <- matrix(NA,nrow=nrow(marker),ncol=length(rv.ix));
      
      for(ii in 1:length(rv.ix))
        {
          ind.ii <- which(rv.count<=(rv.ix[ii]));
          x.mat[,ii] <- rowSums(as.matrix(marker[,ind.ii]));
        }
      
      par.nuis <- NA;
      yub <- 0;ylb <- 0;

      res.C <- .C("genericPerm2",markerPt=as.double(x.mat),
                  yPt=as.double(y.dat),
                  yUbPt=as.double(yub),
                  yLbPt=as.double(ylb),
                  noPersonPt=as.integer(length(y.dat)),
                  noMarkerPt=as.integer(ncol(x.mat)),
                  noPermPt=as.integer(no.perm),
                  pValuePt=as.double(0),
                  testPt=as.integer(2),
                  markerMaxPt=as.double(max(x.mat)),
                  statisticPt=as.double(rep(0,no.perm+1)));

      res <- list(p.value=res.C$pValuePt,
                  statistic=res.C$statisticPt[1],
                  statperm=res.C$statisticPt[-1]);      
      
      return(res);
    }
  
  return(vt.proscore.simple(dat.ped,par.dat,maf.vec,maf.cutoff,no.perm=1000,alternative=alternative));
}



