\name{rsem.Ascov}
\alias{rsem.Ascov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sandwich-type covariance matrix 
}
\description{
Returns the sandwich type covariance matrix. This function is not intended to use seperately from the \code{rsem.emmusig} function.
}
\usage{
rsem.Ascov(xpattern, musig, varphi=.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xpattern}{
Missing data pattern output from \code{rsem.pattern}. 
}
  \item{musig}{
Robust mean and covariance matrix from \code{rsem.emmusig}
}

  \item{varphi}{
  Proportion of data to be down-weighted. Default is 0.1.
}

}
\details{
Data should be a matrix. To change a data frame to a matrix, using data.matrix(x).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Abeta}{A matrix}
  \item{Bbeta}{B matrix}
  \item{Gamma}{Sandwich type covariance matrix}
}


\seealso{
\code{\link{rsem.emmusig}}
}


