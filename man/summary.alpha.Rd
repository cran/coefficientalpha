\name{summary.alpha}
\alias{summary}
\alias{summary.alpha}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print alpha and its confidence interval.
}
\description{
Print alpha and its confidence interval.}
\usage{
\method{summary}{alpha}(object, prob = 0.95,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Results from the function \code{\link{cronbach}}.
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
Zhang, Z. & Yuan, K.-H. (2013). Robust Cronbach's alpha.
}
\examples{
data(holzinger)
spatial<-holzinger[, c('visual', 'cubes', 'paper', 'lozenge', 'paperrev','flagssub')]
verbal<-holzinger[, c('general', 'paragrap', 'sentence', 'wordc', 'wordm')]
speed<-holzinger[, c('add', 'code', 'counting', 'straight')]
memory<-holzinger[, c('wordr', 'numberr', 'figurer', 'object', 'numberf', 'figurew')]

alpha.spatial<-cronbach(spatial)
## diagnostic plot
plot(alpha.spatial, type='d')

## alpha with varphi=.01 & standard error
alpha.spatial<-cronbach(spatial, varphi=.01, se=TRUE)
## confidence interval
summary(alpha.spatial)
}
