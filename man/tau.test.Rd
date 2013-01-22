\name{tau.test}
\alias{tau.test}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test the tau-equivalence and the homogeneity of items
}
\description{
Calculate alpha or omega for a given data set.}
\usage{
tau.test(y, varphi = 0.1, complete = FALSE, drop)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Data
}
  \item{varphi}{
Downweight rate
}

  \item{complete}{
Calculate alpha/omega only based on complete data (listwise deletion)}

  \item{drop}{
The row number of cases to be dropped from the analysis.
}
}

\author{
Zhiyong Zhang and Ke-Hai Yuan
}
\references{
Zhang, Z. & Yuan, K.-H. (2014). Robust Coefficients Alpha and Omega and their Confidence Intervals: Methods and Software.
}
\examples{
data(example)
tau.test(example)
}