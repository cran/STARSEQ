cmc.uniqtl.pop <- function(dat.ped,par.dat,maf.vec,maf.cutoff,no.perm=1000,alternative='greater')
{
  cmc.weight <- function(dat,yub,ylb)
    {
      marker <- as.matrix(dat[,-ncol(dat)]);
      x <- (rowSums(marker)>0);
      return(x);
    }
  
  yub <- NA;ylb <- NA;
  
  ind.rare <- which(maf.vec<maf.cutoff);
  
  ind.rare <- which(maf.vec<maf.cutoff);
  marker <- dat.ped$genotype;
  marker <- marker[,ind.rare];
  
  dat.mat <- cbind(marker,dat.ped$phenotype);
  
  x.dat <- cmc.weight(dat.mat,yub,ylb);
  y.dat <- dat.mat[,ncol(dat.mat)];
  par.nuis <- list();
  stat.dat <- std.score.stat.new(x.dat,y.dat,yub,ylb,par.nuis)
  p.value.greater <- 1-pnorm(stat.dat);
  p.value.less <- pnorm(stat.dat);
  p.value.two.sided <- 1-pchisq(stat.dat^2,df=1);  
  if(alternative=='greater') p.value <- p.value.greater;
  if(alternative=='less') p.value <- p.value.less;
  if(alternative=='two.sided') p.value <- p.value.two.sided;
  
  res <- list(p.value=p.value,
              p.value.greater=p.value.greater,
              p.value.less=p.value.less,
              p.value.two.sided=p.value.two.sided);
  
  return(res);
  
}


