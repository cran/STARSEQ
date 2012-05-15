\name{STAR.test}
\alias{STAR.test}
\title{Secondary Trait Association Analysis with Rare Variants
}
\description{
  Core function for STAR
}
\usage{
STAR.test(dat.ped, y2.vec, covar.mat, rv.count.cutoff, maf.vec, maf.cutoff, pct.upper, pct.lower, yub, ylb, no.perm = 1000, method = "CMC", alternative = "two.sided")
}
\arguments{
  \item{dat.ped}{
    The list of pedigree files. Components include individual.id,
    family.id, genotype, phenotype, father.id, mother.id, and sex. 
  }
  \item{y2.vec}{
    Vector of secondary trait. 
  }
  \item{covar.mat}{
    Covariate matrix.
  }
  \item{rv.count.cutoff}{
    The upper bound for the very rare variants that will be collapsed in
    the analysis. 
  }
  \item{maf.vec}{
    Mandatory minor allele frequencies. 
  }
  \item{maf.cutoff}{
    The minor allele frequency threshold which is used to determine the
    set of variants that are analyzed. 
  }
  \item{pct.upper}{
    Upper percentage of selection.
  }
  \item{pct.lower}{
    Lower percentage of selection. 
  }
  \item{yub}{
    Upper primary trait threshold. 
  }
  \item{ylb}{
    Lower primary trait threshold. 
  }
  \item{no.perm}{
    Number of permutations to be performed. Default choice is
    1000. Adaptive permutation is implemented. For many tests, the
    asymptotic p-values can be obtained,
    however, they may not be accurate. 
  }
  \item{method}{ Supported methods include CMC, WSS (extended by DY Lin
    and ZZ Tang and by Liu and Leal), KBAC, SKAT, and VT. 
  }
  \item{alternative}{ Three possible alternative hypothesis two.sided,
    greater and less. The default choice is two.sided. For SKAT, only
    two sided alternative hypothesis can be tested
  }
}

\value{
  \item{p.value}{P values for the specified test}
  \item{statistic}{Statistic values}
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
  Dajiang J. Liu & Suzanne M. Leal
}