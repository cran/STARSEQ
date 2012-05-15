wss.lin.simple.C <- function(dat.ped,par.dat,maf.vec,maf.cutoff,no.perm=1000,alternative=c('two.sided','greater','less'))
{
  if(length(alternative)>1) alternative <- alternative[1];
  yub <- par.dat$yub;ylb <- par.dat$ylb;
  y.binary.greater <- (dat.ped$phenotype>yub);
  y.binary.less <- (dat.ped$phenotype<ylb);
  y.dat <- dat.ped$phenotype;
  marker <- dat.ped$genotype;
  ix.rm <- which(colSums(marker)==0);
  #print(maf.vec);
  if(length(ix.rm)>0)
    {
      marker <- marker[,-ix.rm];
      maf.vec <- maf.vec[-ix.rm];
    }
 
  ind.rare <- which(maf.vec<maf.cutoff);

  marker <- marker[,ind.rare];
  par.nuis <- NA;
  q <- (colSums(marker)+1)/(2*nrow(marker)+2);
  w <- 1/sqrt(nrow(marker)*q*(1-q));
  #print(w);
  x.dat <- marker%*%(w);
  N <- length(y.dat);
  wss.dat <- as.numeric(sqrt(N)*cor(x.dat,y.dat));
  wss.perm <- 0;
  score.stat.vec <- 0;
  for(ii in 1:ncol(marker))
    {
      score.stat.vec[ii] <- sqrt(N)*cor(marker[,ii],y.dat)
    }
  if(no.perm>0)
    {
      res.C <- .C("genericPerm1",xPt=as.double(x.dat),
         yPt=as.double(y.dat),
         noPerson=as.integer(length(y.dat)),
         noPermPt=as.integer(no.perm),
         pValuePt=as.double(0),
         statisticPt=as.double(0));
      statistic <- res.C$statisticPt;
      p.value <- res.C$pValuePt;
    }
  
  if(no.perm==0)
    {
      if(alternative=='greater') {
        statistic <- wss.dat;
        p.value <- 1-pnorm(statistic);
      }
      if(alternative=='less') {
        statistic <- wss.dat;
        p.value <- pnorm(statistic);
      }
      if(alternative=='two.sided')
        {
          statistic <- wss.dat^2;
          p.value <- 1-pchisq(statistic,df=1);
        }
    }
  res <- list(p.value=p.value,
              statistic=statistic,
              weight=w,
              score.stat.vec=score.stat.vec);
  return(res);
}
lse.wss.misshsq <- function(marker,y,weight)
  {
    x <- as.vector(marker%*%weight);
    var.x <- as.numeric(var(x));
    res.lm <- lm(y ~ x);
    beta1.est <- res.lm$coefficients[2];
    sigma.gae.sq <- beta1.est^2*var.x;
    return(list(beta1.est=beta1.est,
                hsq.est=sigma.gae.sq));
  }
