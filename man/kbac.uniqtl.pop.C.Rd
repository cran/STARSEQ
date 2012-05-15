\name{kbac.uniqtl.pop.C}
\alias{kbac.uniqtl.pop.C}
\title{ KBAC Test for Population-based Studies of Quantitative Trait
}
\description{
This function implements the extended kernel-based adaptive cluster method.
}
\usage{
kbac.uniqtl.pop.C(dat.ped, par.dat, maf.vec, maf.cutoff, no.perm = 1000, alternative = "greater")
}
\arguments{
  \item{dat.ped}{
    A list of ped files.
  }
  \item{par.dat}{
    A list of parameters for ascertainment. The default in an empty list.
  }
  \item{maf.vec}{
    User specified minor allele frequency vector
  }
  \item{maf.cutoff}{
    Upper minor allele frequency cutoff for rare variant analysis
  }
  \item{no.perm}{
    The number of permutations. Set to 0 is recommended for CMC test. 
  }
  \item{alternative}{
    Alternative hypothesis, default choice is two.sided. Other options
    include greater or less.
  }
}
\value{
  \item{p.value}{P-value as determined by the alternative hypothesis tested}
  \item{statistic}{Statistic value for the CMC test}
}
\author{
  Dajiang Liu
}
