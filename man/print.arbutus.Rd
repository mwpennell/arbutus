% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare.summStat.R
\name{print.arbutus}
\alias{print.arbutus}
\title{Print arbutus object}
\usage{
\method{print}{arbutus}(x, ...)
}
\arguments{
\item{x}{an \code{arbutus} object}

\item{...}{additional arguments to be passed to print}
}
\value{
prints p.values to screen
}
\description{
Generic function for printing 'arbutus' objects
from \code{\link{compare_pic_stat}}
}
\examples{
data(finch)
phy <- finch$phy
dat <- finch$data[,"wingL"]

## use arbutus() to assess model adequacy
## pretend as if tree has already been rescaled
p <- arbutus(phy, data=dat, nsim=10)

## print output
p

}
