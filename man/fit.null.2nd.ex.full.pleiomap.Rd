\name{fit.null.2nd.ex.full.pleiomap}
\alias{fit.null.2nd.ex.full.pleiomap}

\title{
  Model Fitting Function for STARSEQ

}
\description{
This implements the STARSEQ model for correcting the secondary trait
distribution for selected samples. This function is not intended for
taking user input. 
}
\usage{
fit.null.2nd.ex.full.pleiomap(marker, rv.count.cutoff, covar.mat, y1.vec, y2.vec, yub, ylb, pct.upper, pct.lower, Nub, Nlb)
}
\arguments{
  \item{marker}{
    A Matrix for marker genotypes.
  }
  \item{rv.count.cutoff}{
    An upper bound for collapsing variants.
  }
  \item{covar.mat}{
    Matrix for covariate
  }
  \item{y1.vec}{
    Primary trait.
  }
  \item{y2.vec}{
    Secondary trait.
  }
  \item{yub}{
    Upper trait threshold for selective sampling. 
  }
  \item{ylb}{
    Lower trait threshold for selective sampling. 
  }
  \item{pct.upper}{
    Upper percentage for selective sampling.
  }
  \item{pct.lower}{
    Lower percentage for selective sampling. 
  }
  \item{Nub}{
    Number of samples from upper extremes.
  }
  \item{Nlb}{
    Number of samples from lower extremes.
  }
}

\value{
  \item{beta.10.est}{Estimates of intercept for the primary trait}
  \item{beta.11.vec.est}{beta.11.vec.est is the estimates for primary
    trait effects}
  \item{beta.20.est}{beta.20.est is the estimated intercept for
    secondary trait}
  \item{ beta.2y1.est}{beta.2y1.est is the estimated secondary trait effect }
  \item{sigma.1.est}{ sigma.1.est is the estimated primary trait
    standard deviation}
  \item{sigma.2.est}{ sigma.2.est is the estimated secondary trait
    standard deviation} 
  \item{beta.1z.est}{ beta.1z.est is the estimated covariate primary
    trait effect}
  \item{beta.2z.est}{beta.2z.est is the estimated covariate secondary
    trait effect} 
  \item{r.y1}{ r.y1 is the primary trait residuals}
  \item{r.y2}{r.y2 is the secondary trait residuals}
}

\author{
  Dajiang Liu
}
