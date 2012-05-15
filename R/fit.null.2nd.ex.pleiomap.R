#when common variants are also used in the model;
#they are listed in the covariates;
#now will provide a variable threshold version of the test;
fit.null.2nd.ex.full.pleiomap <- function(marker,rv.count.cutoff,covar.mat,y1.vec,y2.vec,yub,ylb,pct.upper,pct.lower,Nub,Nlb)
  {
    rv.count <- colSums(marker);

    ind.rare <- which(rv.count<rv.count.cutoff);

    #collapse very rare variants:
    if(length(ind.rare)>0) {
      ind.common <- (1:ncol(marker))[-ind.rare];
      marker.rare <- marker[,ind.rare];
      x.rare <- rowSums(marker.rare)>0;
      if(length(ind.common)>0) {
        marker.common <- marker[,ind.common];
        marker.new <- cbind(x.rare,marker.common);
      }
      if(length(ind.common)==0) {
        marker.new <- as.matrix(x.rare);
      }
    }

    #set up initial parameter:
    beta.10 <- 0;beta.11.vec <- rep(0,ncol(marker.new));beta.20 <- 0;
    beta.2y1 <- 0;log.sigma.1 <- 0;log.sigma.2 <- 0;
    
    res.fit.null.ex <- list();
    pars.init <- c(beta.10,beta.11.vec,beta.20,beta.2y1,log.sigma.1,log.sigma.2);
     
    if(length(covar.mat)!=0)
      {
        beta.1z <- rep(0,ncol(covar.mat));
        beta.2z <- rep(0,ncol(covar.mat));
        pars.init <- c(beta.10,beta.11.vec,beta.20,beta.2y1,log.sigma.1,log.sigma.2,beta.1z,beta.2z);
      }

    #print(length(y2.vec))
    res.fit.null.ex <- nlminb(pars.init,minus.L0.2nd.ex.full.pleiomap,
                              marker=marker.new,
                              covar.mat=covar.mat,
                              y1.vec=y1.vec,
                              y2.vec=y2.vec,
                              yub=yub,
                              ylb=ylb,
                              pct.upper=pct.upper,
                              pct.lower=pct.lower,
                              Nub=Nub,
                              Nlb=Nlb);
    
    pars.est <- res.fit.null.ex$par;

    beta.10.est <- pars.est[1];beta.11.vec.est <- pars.est[2:(1+ncol(marker.new))];
    beta.20.est <- pars.est[2+ncol(marker.new)];
    beta.2y1.est <- pars.est[3+ncol(marker.new)];
    sigma.1.est <- exp(pars.est[4+ncol(marker.new)]);
    sigma.2.est <- exp(pars.est[5+ncol(marker.new)]);
    
    beta.1z.est <- vector(length=0);
    beta.2z.est <- vector(length=0);
    
    
    if(length(covar.mat)!=0) {
      beta.1z.est <- pars.est[(5+ncol(marker.new)+1):(5+ncol(marker.new)+ncol(covar.mat))];
      beta.2z.est <- pars.est[(5+ncol(marker.new)+ncol(covar.mat)+1):(5+ncol(marker.new)+2*ncol(covar.mat))];
    }
    #get residual:

    r.y1 <- y1.vec-beta.10.est-marker.new%*%beta.11.vec.est;
    r.y2 <- y2.vec-beta.20.est-beta.2y1.est*r.y1;

    if(length(covar.mat)>0) {
      r.y2 <- r.y2-covar.mat%*%beta.2z.est;
    }
    
#    r.y2 <- qqnorm(r.y2,plot.it=F)$x;#quantile normalize the residual;

    res <- list(beta.10.est=beta.10.est,
                beta.11.vec.est=beta.11.vec.est,
                beta.20.est=beta.20.est,
                beta.2y1.est=beta.2y1.est,
                sigma.1.est=sigma.1.est,
                sigma.2.est=sigma.2.est,
                beta.1z.est=beta.1z.est,
                beta.2z.est=beta.2z.est,
                r.y1=r.y1,
                r.y2=r.y2);

                                        #    print(c(res$beta.11.est,res$beta.21.est));
    return(res);    
  }

minus.L0.2nd.ex.full.pleiomap <- function(pars,marker,covar.mat,y1.vec,y2.vec,yub,ylb,pct.upper,pct.lower,Nub,Nlb)
{
  #pars.init <- c(beta.10,beta.11,beta.20,beta.2y1,log.sigma.1,log.sigma.2,beta.1z,beta.2z);
  beta.10 <- pars[1];beta.11.vec <- pars[(1+1):(ncol(marker)+1)];

  beta.20 <- pars[ncol(marker)+2];beta.2y1 <- pars[ncol(marker)+3];

  rho <- beta.2y1;
  
  log.sigma.1 <- pars[ncol(marker)+4];
  log.sigma.2 <- pars[ncol(marker)+5];

  sigma.1 <- exp(log.sigma.1);sigma.2 <- exp(log.sigma.2);sigma.1.sq <- sigma.1*sigma.1;sigma.2.sq <- sigma.2*sigma.2;

  prob.y1.upper <- pct.upper;
  prob.y1.lower <- pct.lower;
  
  if(length(covar.mat)>0) {
    beta.1z <- pars[(ncol(marker)+5+1):(ncol(marker)+5+ncol(covar.mat))];
    beta.2z <- pars[(ncol(covar.mat)+ncol(marker)+5+1):(ncol(marker)+5+2*ncol(covar.mat))];
    y1.vec <- y1.vec-covar.mat%*%beta.1z;
    y2.vec <- y2.vec-covar.mat%*%beta.2z;}
  
  #print(length(y2.vec));
  #compute likelihoods:

  r.y1 <- y1.vec-beta.10-marker%*%beta.11.vec;
  
  prob.y1.y2.given.x <- dnorm(y2.vec,beta.20+sigma.2/sigma.1*rho*r.y1,sigma.2*sqrt(1-rho^2))*dnorm(y1.vec,beta.10+marker%*%beta.11.vec,sigma.1);
  prob.spl.given.y1 <- (y1.vec>yub)*Nub/prob.y1.upper+(y1.vec<ylb)*Nlb/prob.y1.lower;

  prob.spl.upper <- Nub/pct.upper;
  prob.spl.lower <- Nlb/pct.lower;
  
  prob.y1.upper.given.x <- 1-pnorm(yub,beta.10+marker%*%beta.11.vec,sigma.1);
  prob.y1.lower.given.x <- pnorm(ylb,beta.10+marker%*%beta.11.vec,sigma.1);
 
  numer <- sum(log(prob.spl.given.y1*prob.y1.y2.given.x));
  denom <- sum(log(prob.spl.upper*prob.y1.upper.given.x+
                   prob.spl.lower*prob.y1.lower.given.x));
  return(-numer+denom);  
}
