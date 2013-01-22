\name{alpha}
\alias{alpha}
\alias{omega}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust Cronbach's alpha and McDonald's omega
}
\description{
Calculate alpha or omega for a given data set.}
\usage{
alpha(y, varphi = 0.1, se = FALSE, complete = FALSE, drop)
omega(y, varphi = 0.1, se = FALSE, complete = FALSE, drop)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Data
}
  \item{varphi}{
Downweight rate
}
  \item{se}{
Whether to estimate standard error. It might take significant more time if se is requested with many items.
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

alpha(example, varphi=.01)
omega(example, varphi=.01)
}
