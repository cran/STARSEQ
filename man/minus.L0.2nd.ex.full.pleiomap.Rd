\name{minus.L0.2nd.ex.full.pleiomap}
\alias{minus.L0.2nd.ex.full.pleiomap}
\title{
  Minus Log-likelihoood Function for STAR
}
\description{
  Minus Log-likelihoood Function for STAR
}
\usage{
minus.L0.2nd.ex.full.pleiomap(pars, marker, covar.mat, y1.vec, y2.vec, yub, ylb, pct.upper, pct.lower, Nub, Nlb)
}
\arguments{
  \item{pars}{
    Vector of model paramters
}
  \item{marker}{
    Genotype matrix
  }
  \item{covar.mat}{
    Matrix of covariates
  }
  \item{y1.vec}{
    Vector of primary traits
  }
  \item{y2.vec}{
    Vector of secondary traits
  }
  \item{yub}{
    Upper extreme primary trait cutoff for selective sampling
  }
  \item{ylb}{
    Lower extreme primary trait cutoff for selective sampling
  }
  \item{pct.upper}{
    Proportion of samples from the upper extreme
  }
  \item{pct.lower}{
    Proportion of samples from the lower extreme
  }
  \item{Nub}{
    Number of samples from the upper extreme of the primary traits
  }
  \item{Nlb}{
    Number of samples from lower extreme of the primary traits
  }
}

\value{
  Returns the value for minus log likelihood, which will be used for
  maximum likelihood estimates. This function is not intended to take
  user inputs. 
}
\author{
  Dajiang Liu
}
