% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lalonde.R
\docType{data}
\name{lalonde}
\alias{lalonde}
\title{Lalonde data with external controls}
\format{
An object of class \code{list} of length 7.
}
\usage{
lalonde
}
\description{
The data are drawn from a paper by Robert Lalonde, "Evaluating the Econometric Evaluations of Training Programs," American Economic Review, Vol. 76, pp. 604-620, and published by Dehejia and Wabha on \url{http://users.nber.org/~rdehejia/data/nswdata2.html}.
}
\details{
Each entry is a data.frame with name "cps", "cps2", "cps3", "psid", "psid2", "psid3" and "rct". "rct" corresponds to the dataset from the randomized experiment with 185 treated and 260 control units (with earnings in 1974 available). The other datasets append the corresponding external controls to the dataset given by "dw". See \url{http://users.nber.org/~rdehejia/data/nswdata2.html} for details.
}
\references{
\insertRef{lalonde1986evaluating}{ovalue}

\insertRef{dehejia1999causal}{ovalue}
}
\keyword{datasets}
