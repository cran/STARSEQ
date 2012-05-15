meta.p <- function(p.value.vec,N.vec,no.perm.vec,alternative=c('two.sided','greater','less'))
  {
    if(length(alternative)>1)
      alternative <- 'two.sided';
    for(ii in 1:length(no.perm.vec))
      {
        if(p.value.vec[ii]==0 & no.perm.vec[ii]>0)
          p.value.vec[ii] <- 1/(no.perm.vec[ii]);
      }
    
    z.score.vec <- qnorm(1-p.value.vec);
    statistic <- sum(sqrt(N.vec)*z.score.vec)/sqrt(sum(N.vec));
    if(alternative=='greater')
      p.value <- 1-pnorm(statistic);
    if(alternative=='less')
      p.value <- pnorm(statistic);
    if(alternative=='two.sided')
      {
        p.value <- 1-pchisq(statistic^2,df=1);
        statistic <- statistic^2;
      }

    return(list(p.value=p.value,
                statistic=statistic));

  }
