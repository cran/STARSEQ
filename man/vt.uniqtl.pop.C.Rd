\name{vt.uniqtl.pop.C}
\alias{vt.uniqtl.pop.C}
\title{VT Test for Population-basd Studies of Quantitative Trait
}
\description{
This function implements the variable threshold test 
}
\usage{
vt.uniqtl.pop.C(dat.ped, par.dat, maf.vec, maf.cutoff, no.perm = 1000, alternative = c("two.sided","greater","less"))
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
    The number of permutations. Set to 1000 is default for VT
  test. Adaptive permutatoin is implemented 
}
  \item{alternative}{
    Alternative hypothesis, default choice is two.sided. Other options
  include greater or less.
}
}

\value{
  \item{p.value}{P-value as determined by the alternative hypothesis tested}
  \item{statistic}{Statistic value for the VT test}
}
\author{
Dajiang Liu
}
