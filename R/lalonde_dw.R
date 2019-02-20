#' Lalonde data with external controls
#'
#' The data are drawn from a paper by Robert Lalonde, "Evaluating the Econometric Evaluations of Training Programs," American Economic Review, Vol. 76, pp. 604-620, and published by Dehejia and Wabha on \url{http://users.nber.org/~rdehejia/data/nswdata2.html}.
#'
#' Each entry is a data.frame with name "cps", "cps2", "cps3", "psid", "psid2", "psid3" and "dw". "dw" corresponds to the dataset from the RCT with 185 treated and 260 control units (with earnings in 1974 available). The other datasets append the corresponding external controls to the dataset given by "dw". See \url{http://users.nber.org/~rdehejia/data/nswdata2.html} for details. 
#'
#' @references
#' \insertRef{lalonde1986evaluating}{ovalue}
#'
#' \insertRef{dehejia1999causal}{ovalue}
#' 
"lalonde_dw"
