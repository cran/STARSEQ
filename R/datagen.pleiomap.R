mat2ped <- function(dat.mat)
{

    genotype <- dat.mat[,1:(ncol(dat.mat)-2)];
    phenotype<- dat.mat[,ncol(dat.mat)-1];
    family.id <- paste('FAM',1:nrow(dat.mat),sep='');
    individual.id <- paste('individiual',1:nrow(dat.mat),sep='');
    sex <- rep('unknown',nrow(dat.mat));
    father.id <- paste('FATHER',1:nrow(dat.mat),sep='');
    mother.id <- paste('MOTHER',1:nrow(dat.mat),sep='');
    dat.ped <- list(family.id=family.id,
                    individual.id=individual.id,
                    father.id=father.id,
                    mother.id=mother.id,
                    genotype=genotype,
                    phenotype=phenotype,
                    sex=sex)
    return(dat.ped);
}
datagen.pleiomap <- function(pars)
{
    scenario <- pars$scenario;
    if(scenario=='extreme-trait')
      {
        maf.vec <- pars$maf.vec;
        beta.10 <- 0;
        beta.11 <- pars$beta.11;
        beta.20 <- 0;
        beta.21.vec <- pars$beta.21.vec;
        sigma.mat <- pars$sigma.mat;
        sigma.11 <- sigma.mat[1,1];
        sigma.22 <- sigma.mat[2,2];
        sigma.12 <- sigma.mat[1,2];
        rho <- sigma.12/sqrt(sigma.11*sigma.22);
        sigma.1 <- sqrt(sigma.11);
        sigma.2 <- sqrt(sigma.22);
        pool.size <- pars$pool.size;
        no.sample <- pars$sample.size;
        pct.upper <- pars$pct.upper;
        pct.lower <- pars$pct.lower;
        pct.causal <- pars$pct.causal;
        ##        ind.noncausal <- pars$ind.noncausal;
        ##        ind.causal <- (1:length(maf.vec))[-ind.noncausal];
        ##        ind.causal.1 <- integer(0);ind.causal.2 <- integer(0);
        ##        ind.causal.1 <- sample(ind.causal,as.integer(pct.causal*length(ind.causal)));
        ##        ind.causal.2 <- sample(ind.causal,as.integer(pct.causal*length(ind.causal)));
        ##Print(dim(hapmat));
        ind.causal.1 <- pars$ind.causal.1;
        ind.causal.2 <- pars$ind.causal.2;
        
        case.count <- 0;ctrl.count <- 0;
        dat.mat.pool <- matrix(nrow=pool.size,ncol=length(maf.vec)+2);      
        count <- 0;
        while(count<pool.size)
          {
            hap1 <- (runif(length(maf.vec))<maf.vec);
            hap2 <- (runif(length(maf.vec))<maf.vec);
            geno <- hap1+hap2;
            x.1 <- 0;x.2 <- 0;
            if(length(ind.causal.1)>0 && length(ind.causal.2)>0) {
              x.1 <- sum((beta.11*geno)[ind.causal.1]);
              x.2 <- sum((beta.21.vec*geno)[ind.causal.2]);
            }
            y.1 <- rnorm(1,x.1,sqrt(sigma.11));
            y.2 <- rnorm(1,beta.20+x.2+sigma.2/sigma.1*rho*(y.1-x.1),sqrt(sigma.22*(1-rho^2)));    
            count <- count+1;
            dat.quant.ins <- c(geno,y.1,y.2);#insert into dat.ori;
            dat.mat.pool[count,] <- dat.quant.ins;
          }
        y1.vec <- dat.mat.pool[,ncol(dat.mat.pool)-1];
        yub <- quantile(y1.vec,pct.upper);
        ylb <- quantile(y1.vec,pct.lower);
        dat.mat.quant <- dat.mat.pool[c(which(y1.vec>yub),which(y1.vec<ylb)),];
        dat.ped.quant <- mat2ped(dat.mat.quant);
        total.no.sites <- length(maf.vec);
        sumstat <- list(yub=yub,
                        ylb=ylb,
                        total.no.causal.sites.1=length(ind.causal.1),
                        total.no.causal.sites.2=length(ind.causal.2),
                        total.causal.var.freq.1=sum(maf.vec[ind.causal.1]),
                        total.causal.var.freq.2=sum(maf.vec[ind.causal.2]),
                        total.no.sites=total.no.sites,
                        total.no.pleiotropic.sites=length(set.intersect(ind.causal.1,ind.causal.2)),
                        total.pleiotropic.var.freq=sum(maf.vec[set.intersect(ind.causal.1,ind.causal.2)]),
                        ind.causal.1=ind.causal.1,
                        ind.causal.2=ind.causal.2,
                        maf.vec=maf.vec,
                        hsq.1=beta.11^2*sum(maf.vec[ind.causal.1]*(1-maf.vec[ind.causal.1])),
                        hsq.2=sum((beta.21.vec[ind.causal.2])^2*maf.vec[ind.causal.2]*(1-maf.vec[ind.causal.2])),
                        par.sim=pars);

        res <- list(dat.ped=dat.ped.quant,
                    y.2=dat.mat.quant[,ncol(dat.mat.quant)],
                    sumstat=sumstat);
        return(res);
      }
  }
