\name{skat.uniqtl.simple.C}
\alias{skat.uniqtl.simple.C}
\title{SKAT Test for Population-basd Studies of Quantitative Trait
}
\description{
This function implements the sequence kernel association test. 
}
\usage{
skat.uniqtl.simple.C(dat.ped, par.dat, maf, maf.cutoff, no.perm =
1000, alternative = "two.sided" , out.type="C")
}
\arguments{
  \item{dat.ped}{
    A list of ped files.
}
  \item{par.dat}{
    A list of parameters for ascertainment. The default in an empty list.
}
  \item{maf}{
    User specified minor allele frequency vector
}
  \item{maf.cutoff}{
    Upper minor allele frequency cutoff for rare variant analysis
}
  \item{no.perm}{
    The number of permutations. Set to 1000 is default for SKAT
    test. Adaptive permutatoin is implemented 
  }
  \item{alternative}{
    Alternative hypothesis, default choice is two.sided. Other options
    include greater or less.
  }
  \item{out.type}{C for continuous trait}
}
\value{
  \item{p.value}{P-value as determined by the alternative hypothesis tested}
  \item{statistic}{Statistic value for the SKAT test}
}

\author{
Dajiang Liu
}
