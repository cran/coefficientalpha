\name{rsem.pattern}
\alias{rsem.pattern}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtaining missing data patterns
}
\description{
This function obtains the missing data patterns and the number of cases in each patterns. It also tells the number of observed variables and their indices for each pattern.
}
\usage{
rsem.pattern(x, print=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix as data
}
  \item{print}{
Whether to print the missing data pattern. The default is FALSE.
}
}
\details{
The missing data pattern matrix has 2+p columns. The first column is the number cases in that pattern. The second column is the number of observed variables. The last p columns are a matrix with 1 denoting observed data and 0 denoting missing data.

In addition, a matrix of 0/1 is also used to indicate missing data. 1 means missing and 0 means observed. 
}
\value{
  \item{x}{Data ordered according to missing data pattern}
  \item{misinfo}{Missing data pattern matrix}
  \item{mispat}{Missing data pattern in better readable form.}
  \item{y}{The original data.}
}
\references{
Yuan, K.-H., & Zhang, Z. (2012). Robust Structural Equation Modeling with Missing Data and Auxiliary Variables. Psychometrika, 77(4), 803-826.
}
\author{
Zhiyong Zhang and Ke-Hai Yuan
}

