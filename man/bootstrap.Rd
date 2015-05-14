\name{bootstrap}
\alias{bootstrap}
\title{
Bootstrap se and CI for alpha and omega
}
\description{
Bootstrap se and CI for alpha and omega.}
\usage{
bootstrap(y, type="omega", alpha=.95, nboot=1000, ci="bc", plot=FALSE, 
varphi=0, complete=FALSE, auxiliary=NULL, silent=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Data
}
  \item{type}{
\code{omega}: coefficient omega. Otherwise, coefficient alpha.
}
  \item{alpha}{
Confidence level.
}
  \item{nboot}{
Number of bootstrap samples to use
}
  \item{ci}{
\code{bc}: Bias-corrected CI. Otherwise, the percentile CI is used.
}
  \item{plot}{
Whether to plot the bootstrap density.
}
  \item{varphi}{
Downweight rate
}
  \item{complete}{
Calculate alpha/omega only based on complete data (listwise deletion)}
\item{auxiliary}{
Provide a matrix or data frame of auxiliary variables for missing data analysis.
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

boot.test<-bootstrap(example, type='alpha', nboot=10, plot=TRUE)
}
