#test code;
#source('datagen.uniqtl.R');source('kbac.uniqtl.R');source('score.stat.R');pars <- list(maf=rep(0.005,15),ind.noncausal=1:5,yub=qnorm(0.7),ylb=qnorm(0.3),Nub=300,Nlb=300,beta1=-0.5,scenario='beta1',maf.cutoff=0.1);res.datagen <- datagen.uniqtl(pars);dat.ped <- res.datagen$dat.sim;par.dat=res.datagen$par.dat;maf.vec <- pars$maf;maf.cutoff <- pars$maf.cutoff;res.kbac <- kbac.uniqtl(dat.ped,par.dat,maf.vec,maf.cutoff,1000,'less');print(res.kbac);


kbac.uniqtl.pop.C <- function(dat.ped,par.dat,maf.vec,maf.cutoff,no.perm=1000,alternative='greater')
{

    ind.rare <- which(maf.vec<maf.cutoff);

    marker <- dat.ped$genotype;
    marker <- marker[,ind.rare];

    y <- dat.ped$phenotype;

    yub.pct <- par.dat$yub.pct;
    ylb.pct <- par.dat$ylb.pct;
    
    if(length(yub.pct)==0) yub.pct <- 0.75;
    if(length(ylb.pct)==0) ylb.pct <- 0.25;
    
    yub <- as.numeric(quantile(y,yub.pct));
    ylb <- as.numeric(quantile(y,ylb.pct));
    y.binary.greater <- rep(NA,length(y));
    y.binary.less <- rep(NA,length(y));

    y.binary.greater[which(y>=yub)] <- 1;
    y.binary.greater[which(y<=yub)] <- 0;

    y.binary.less[which(y>=ylb)] <- 0;
    y.binary.less[which(y<=ylb)] <- 1;

    dat.binary.greater <- cbind(marker,y.binary.greater);
    dat.binary.less <- cbind(marker,y.binary.less);
    dat.mat<- cbind(marker,y);
    dat.rec <- recode.uniqtl(dat.mat);

    dat.rec.binary.greater <- cbind(dat.rec[,1],y.binary.greater);
    dat.rec.binary.less <- cbind(dat.rec[,1],y.binary.less);

    m.dat <- dat.rec[,1];

    x.dat.greater <- kbac.weight(dat.rec.binary.greater);
    x.dat.less <- kbac.weight(dat.rec.binary.less)

    par.nuis <- NA;

    kbac.dat.greater <- std.score.stat.new(x.dat.greater,y,yub,ylb,par.nuis);
    kbac.dat.less <- std.score.stat.new(x.dat.less,y,yub,ylb,par.nuis);

    kbac.perm.greater <- vector(length=no.perm);
    kbac.perm.less <- vector(length=no.perm);

    res.C <- .C("genericPerm2",markerPt=as.double(dat.rec[,1]-1),
                yPt=as.double(y),
                yUbPt=as.double(yub),
                yLbPt=as.double(ylb),
                noPersonPt=as.integer(length(y)),
                noMarkerPt=as.integer(1),
                noPermPt=as.integer(no.perm),
                pValuePt=as.double(0),
                testPt=as.integer(1),
                markerMaxPt=as.double(max(dat.rec[,1])),
                statisticPt=as.double(rep(0,no.perm+1)));
    
    res <- list(p.value=res.C$pValuePt,
                statistic=res.C$statisticPt[1],
                statperm=res.C$statisticPt[-1]);

    return(res);
}



recode.uniqtl <- function(dat)
{

    marker <- as.matrix(dat[,-ncol(dat)]);
    marker <- (marker>0);
    marker.str <- 0;

    for(ii in 1:nrow(marker))
    {
        marker.str[ii] <- paste(as.integer(marker[ii,]),collapse='',sep='');
    }
  marker.pool <- paste(rep(0,ncol(marker)),collapse='',sep='');
  marker.rec <- 1;
  for(ii in 1:nrow(marker))
  {
    rec <- which(marker.str[ii]==marker.pool);
    if(length(rec)<1) marker.pool <- c(marker.pool,marker.str[ii]);
    rec <- which(marker.str[ii]==marker.pool);
    marker.rec[ii] <- rec;
  }
  return(cbind(marker.rec,dat[,ncol(dat)]))
}




kbac.weight <- function(dat)
{
  k <- max(dat[,1]);
  case.count <- 0;
  control.count <- 0;

  for(ii in 1:k)
  {
    case.count[ii] <- sum(dat[,1]==ii & dat[,2]==1);
    control.count[ii] <- sum(dat[,1]==ii & dat[,2]==0);
  }
  nocol <- ncol(dat);
  norow <- nrow(dat);

  weight <- phyper(case.count,sum(case.count),sum(control.count),case.count+control.count);

  if(length(which((case.count+control.count)==0))>0) {
    weight[which((case.count+control.count)==0)] <- 1/2;
  }
  weight <- weight;
  weight[1] <- 0;
  #return the weighted genotype coding;
  x <- weight[dat[,1]];
  return(x);

}
