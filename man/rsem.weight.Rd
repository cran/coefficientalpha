\name{rsem.weight}
\alias{rsem.weight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate weight for each subject
}
\description{
Calculate weight for each subject in estimating the mean and covariance matrix.
}
\usage{
rsem.weight(x, varphi, mu0, sig0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data}
  \item{varphi}{
Downweight rate.
}
  \item{mu0}{
Robust mean
}
  \item{sig0}{
Robust covariance matrix.
}
}
\value{
\item{w1}{Weight for robust mean estimates}
\item{w2}{Weight for robust covariance estimates}
}
\references{
Yuan, K.-H., & Zhang, Z. (2012). Robust Structural Equation Modeling with Missing Data and Auxiliary Variables. Psychometrika, 77(4), 803-826.
}
\author{
Zhiyong Zhang and Ke-Hai Yuan
}

