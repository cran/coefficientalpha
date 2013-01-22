\name{summary.alpha}
\alias{summary}
\alias{summary.alpha}
\alias{summary.omega}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print alpha/omega and its confidence interval.
}
\description{
Print alpha/omega and its confidence interval.}
\usage{
\method{summary}{alpha}(object, prob = 0.95,...)
\method{summary}{omega}(object, prob = 0.95,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Results from the function \code{\link{alpha}} or \code{\link{omega}}.
}
  \item{prob}{
Alpha level for confidence interval. By default, the 95% CI is calculated.
}
  \item{...}{Option can be passed to the summary function.}
}
\author{
Zhiyong Zhang and Ke-Hai Yuan
}
\references{
Zhang, Z. & Yuan, K.-H. (2014). Robust Coefficients Alpha and Omega and their Confidence Intervals: Methods and Software.
}
\examples{
data(example)

res<-alpha(example, varphi=.01)
## diagnostic plot
plot(res, type='d')

## alpha with varphi=.01 & standard error
res<-alpha(example, varphi=.01, se=TRUE)
## confidence interval
summary(res)
}
