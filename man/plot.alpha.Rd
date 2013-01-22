\name{plot.alpha}
\alias{plot.alpha}
\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Alpha related plot
}
\description{
Generate plot related to alpha. Three plot can be generated. (1) The weight plot will plot the weight associated with alpha calculation. (2) Profile plot will plot the cases with smallest weights and the average value. (3) The diagnostic plot plots the alpha according to different values of the tuning parameter.}
\usage{
\method{plot}{alpha}(x, type="weight", profile=5, interval=0.01,center=TRUE, scale=FALSE, w1=FALSE, numbered=FALSE, pos="topright",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Results from the function \code{\link{cronbach}}.
}
  \item{type}{
Three types of plots can be generated. \code{type="weight"} generates the plot of weight of each case.
\code{type="profile"} generates a profile plot for the smallest weight. 
\code{type="weight"} generates the diagnostic plot for \code{varphi} with an interval of .01. 
}
  \item{profile}{
Number of cases used on the profile plot. At most 10.}
  \item{interval}{
The interval used in the diagnostic plot. The default is 0.01.}
\item{center}{
 Whether to center the data in the profile plot.
}
\item{scale}{
 Whether to scale the data using variance parameters in the profile plot.
}
\item{w1}{Whether to plot the weight for means}
\item{numbered}{Whether to number the profile plot}
  \item{pos}{Position of legend. If pos=NULL, no legend is plotted.}
\item{...}{Options can be passed to the \code{\link{plot}} function.}
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

## weight plot
plot(alpha.spatial)
# or
plot(alpha.spatial, type='w')

## profile plot
plot(alpha.spatial, type='p')


alpha.verbal<-cronbach(verbal, varphi=.02, se=TRUE)
## confidence interval
summary(alpha.verbal)

## weight plot
plot(alpha.verbal)

## profile plot
plot(alpha.verbal, type='p', profile=6)

alpha.speed<-cronbach(speed, varphi=.01, se=TRUE)
## confidence interval
summary(alpha.speed)

## weight plot
plot(alpha.speed)

## profile plot
plot(alpha.speed, type='p', pos='topleft')

alpha.memory<-cronbach(memory, varphi=.04, se=TRUE)
## confidence interval
summary(alpha.memory)

## weight plot
plot(alpha.memory)

## profile plot
plot(alpha.memory, type='p', profile=10, pos='bottomright')
 plot(alpha.memory, type='p', profile=10, scale=TRUE, pos='bottomright')
}
