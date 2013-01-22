\name{plot.alpha}
\alias{plot.alpha}
\alias{plot.omega}
\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Alpha/Omega related plot
}
\description{
Generate plot related to alpha or omega. Three plot can be generated. (1) The weight plot will plot the weight associated with alpha calculation. (2) Profile plot will plot the cases with smallest weights and the average value. (3) The diagnostic plot plots the alpha according to different values of the tuning parameter.}
\usage{
\method{plot}{alpha}(x, type="weight", profile=5, interval=0.01,
center=TRUE, scale=FALSE, w1=FALSE, numbered=FALSE, pos="topright",...)
\method{plot}{omega}(x, type="weight", profile=5, interval=0.01,
center=TRUE, scale=FALSE, w1=FALSE, numbered=FALSE, pos="topright",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Results from the function \code{\link{alpha}} or \code{\link{omega}}.
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

## weight plot
plot(res)
# or
plot(res, type='w')

## profile plot
plot(res, type='p')
}
