\name{rsem.emmusig}
\alias{rsem.emmusig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust mean and covariance matrix using Huber-type weight
}
\description{
Robust mean and covariance matrix using Huber-type weight.
}
\usage{
rsem.emmusig(xpattern, varphi=.1, max.it=1000, st='i')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xpattern}{
Missing data pattern output from \code{rsem.pattern}. 
}
  \item{varphi}{
  Proportion of data to be down-weighted. Default is 0.1.
}
  \item{max.it}{Maximum number of iterations for EM. Default is 1000}
  \item{st}{Starting values for EM algorithm. The default is 0 for mean and I for covariance. Alternative, the starting values can be estimated according to MCD.}
}
\details{
Estimate mean and covariance matrix using the expectation robust (ER) algorithm.
}
\value{

  \item{err}{Error code. 0: good. 1: maximum iterations are exceeded.}
  \item{mu}{Mean vector}
  \item{sigma}{Covariance matrix}
  \item{weight}{weight used in robust mean and covariance estimation.}
}
\references{
Yuan, K.-H., & Zhang, Z. (2012). Robust Structural Equation Modeling with Missing Data and Auxiliary Variables. Psychometrika, 77(4), 803-826.
}
\author{
Zhiyong Zhang and Ke-Hai Yuan
}

\seealso{
\code{\link{rsem.emmusig}}
}

