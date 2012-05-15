library('numDeriv');






std.score.stat.new <- function(x,y,yub,ylb,par.nuis)
{
    x.mean <- mean(x);
    y.mean <- mean(y);
    N <- length(x);
    numer <- sum((x-x.mean)*y);
    denom <- sqrt(1/N*sum((x-x.mean)^2)*sum((y-y.mean)^2));
    return(numer/denom);
}


std.score.proscore.stat.mat <- function(x.mat,y,yub,ylb,par.nuis)
  {
    N <- ncol(x.mat);
    stat.vec <- vector(length=N);
    for(ii in 1:N)
      {
        x <- x.mat[,ii];
        stat.vec[ii] <- std.score.stat.new(x,y,yub,ylb,par.nuis);
      }
    return(stat.vec);
  }
