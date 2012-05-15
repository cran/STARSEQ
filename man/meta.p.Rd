\name{meta.p}
\alias{meta.p}

\title{
P-value Based Meta-Analysis
}
\description{
Meta-Analysis Combining P-values
}
\usage{
meta.p(p.value.vec, N.vec, no.perm.vec, alternative = c("two.sided", "greater", "less"))
}

\arguments{
  \item{p.value.vec}{
    A vector of p-values from each study;
  }
  \item{N.vec}{
    Sample size from each study
  }
  \item{no.perm.vec}{
    If p-values are permutation based, the number of permutations has to
    be specified
  }
  \item{alternative}{
    Alternative hypothesis, which can be two.sided, greater or less. The
    default choice is two.sided.
  }
}
\value{
  \item{p.value}{P value of meta analysis}
  \item{statistic}{Meta-analysis test statisitcs. The statistic is
    obtained by combining Z-score statistics, which are transformed from
    p-values and weighted by the square root of the sample size. 
  }
} 
\author{
  Dajiang Liu
}
