skat.uniqtl.simple.C <- function(dat.ped,par.dat,maf,maf.cutoff,no.perm=1000,alternative='two.sided',out.type="C")
  {
    marker <- dat.ped$genotype;
    ind.rm <- which(colSums(marker)==0);

    if(length(ind.rm)>0) {
      marker <- marker[,-ind.rm];
      maf <- maf[-ind.rm];
    }
    
    ind.rare <- which(maf<maf.cutoff);
    marker <- marker[,ind.rare];

    y <- dat.ped$phenotype;
    marker <- as.matrix(marker);

    stat.perm <- vector(length=no.perm);
    if(no.perm==0) {
      p.value.twosided <- res.skat$p.value;
      statistic <- res.skat$Q;
    }

    W <- diag(dbeta(colSums(marker)/nrow(marker)/2,1,25));
    res.skat <- mySKAT(marker,y,W);
    stat.dat <- res.skat$statistic;
    yub <- 0;ylb <- 0;

    if(no.perm>0)
      {
        markerPt <- marker%*%W%*%t(marker);
        res.C <- .C("genericPerm2",
                    markerPt=as.double(markerPt),
                    yPt=as.double(y),
                    yUbPt=as.double(yub),
                    yLbPt=as.double(ylb),
                    noPersonPt=as.integer(length(y)),
                    noMarkerPt=as.integer(length(y)),
                    noPermPt=as.integer(no.perm),
                    pValuePt=as.double(0),
                    testPt=as.integer(4),
                    markerMaxPt=as.double(2),
                    statisticPt=as.double(rep(0,no.perm+1)));
        statistic <- res.skat$Q
        p.value.twosided <- res.C$pValuePt;
      }
    
    score.stat.vec <- 0;
    N <- nrow(marker);
    for(ii in 1:ncol(marker))
      {
        score.stat.vec[ii] <- sqrt(N)*cor(marker[,ii],y)
      }
    res <- list(p.value=p.value.twosided,
                p.value.greater=NA,
                p.value.less=NA,
                p.value.two.sided=p.value.twosided,
                statistic=statistic,
                score.stat.vec=score.stat.vec);
    return(res);
  }


mySKAT <- function(marker,y,W)
  {
    statistic <- t(y-mean(y))%*%(marker%*%W%*%t(marker))%*%(y-mean(y));
    N <- nrow(marker);
    X.T.times.X <-cov(marker)*(N-1);
    var.Y <- as.numeric(var(y)*(N-1)/N);
    lambda <- eigen(var.Y*(W%*%X.T.times.X))$values;
    p.value <- davies(statistic,lambda=lambda)$Qq;
    return(list(statistic=statistic,
                p.value=p.value));
  }
