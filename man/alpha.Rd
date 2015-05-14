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
alpha(y, varphi = 0.1, se = FALSE, test = TRUE, 
complete =FALSE, auxiliary = NULL, drop, silent = TRUE)
omega(y, varphi = 0.1, se = FALSE, test = TRUE, 
complete =FALSE, auxiliary = NULL, drop, silent = TRUE)
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
\item{test}{
Whether to test the assumption of alpha (tau equivalent test) or omega (homogeneous test).
}
  \item{complete}{
Calculate alpha/omega only based on complete data (listwise deletion)}
\item{auxiliary}{
Provide a matrix or data frame of auxiliary variables for missing data analysis.
}
\item{drop}{
The row number of cases to be dropped from the analysis.
}
\item{silent}{
Whether to print information of the analysis.
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
