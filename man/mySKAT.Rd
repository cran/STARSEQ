\name{mySKAT}
\alias{mySKAT}

\title{
  Implementation of SKAT Test
}
\description{
This implements the SKAT test according to Wu et al. 
}
\usage{
mySKAT(marker, y, W)
}
\arguments{
  \item{marker}{
    Marker genotypes
}
  \item{y}{
    Vector of phenotypes
  }
  \item{W}{
    Weights calculated based the kernel used. 
  }
}
\value{
  \item{p.value}{p values from SKAT calculated based upon mixture
    chi-square distribution.} 
}

\author{
  Dajiang Liu
}


