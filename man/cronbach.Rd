\name{cronbach}
\alias{cronbach}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust Cronbach's alpha
}
\description{
Calculate alpha for a given data set.}
\usage{
cronbach(y, varphi = 0.1, se = FALSE, complete = FALSE)
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
Whether to estimate standard error. It will take significant more time if se is requested.
}
  \item{complete}{
Calculate alpha only based on complete data (lisewise deletion)}
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

alpha.spatial<-cronbach(spatial, varphi=.01)
}