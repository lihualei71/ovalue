## Log upper-tail inequalities of U statistics
ustat_hoeffding_upper <- function(mu, x, n, k){
    n <- floor(n / k)
    -n * h1(pmin(mu, x), mu)
}

ustat_bentkus_upper <- function(mu, x, n, k){
    n <- floor(n / k)
    log(pbinom(ceiling(n * x), n, mu, lower.tail = TRUE)) + 1
}

ustat_maurer_upper <- function(mu, x, n, k){
    lambda_fun <- function(lambda){
        Glambda <- exp(lambda) - lambda - 1
        if (is.nan(Glambda)){
            Glambda <- Inf
        }
        lambda * (x - mu / (k * Glambda / lambda + 1))
    }
    left <- 1e-5
    right <- k
    while (TRUE){
        obj <- optimize(lambda_fun, c(1e-5, k))
        if (obj$minimum >= k){
            right <- 10 * right
            left <- right
            obj <- optimize(lambda_fun, c(left, right))
        } else {
            break
        }
    }
    logprob <- n * obj$objective / k
    logprob
}

#' Upper confidence bound of mean of U statistics of order 2
#' via Hoeffding-Bentkus-Maurer inequalities
#'
#' @param U the value of the U statistic
#' @param n sample size
#' @param alpha confidence level
#' 
HBM_U_upper <- function(U, n, alpha){
    m <- floor(n / 2)
    tailprob <- function(EU){
        hoeffding_EU <- ustat_hoeffding_upper(EU, U, n, 2)
        bentkus_EU <- ustat_bentkus_upper(EU, U, n, 2)
        maurer_EU <- ustat_maurer_upper(EU, U, n, 2)
        min(hoeffding_EU, bentkus_EU, maurer_EU) - log(alpha)
    }
    if (tailprob(1 - 1e-4) > 0){
        1
    } else {
        uniroot(tailprob, c(U + 1e-10, 1 - 1e-4))$root
    }
}

## Exact DiR O-value
DiR_ovalue_exact <- function(score1, score0, delta){
    n1 <- length(score1)
    n0 <- length(score0)
    n <- n1 + n0
    delta <- delta / 4

    score1 <- pmin(score1 + 1e-10 * runif(n1), 1)
    score0 <- pmin(score0 + 1e-10 * runif(n0), 1)

    ustat_left <- 2 * sum(fast_rank(score0, score1)) / n / (n - 1) # for p*
    EU_left_upper <- HBM_U_upper(ustat_left, n, delta)
    ustat_right <- 2 * sum(fast_rank(score1, score0)) / n / (n - 1) # for 1 - p*
    EU_right_upper <- HBM_U_upper(ustat_right, n, delta)

    ovalue_fun <- function(pi, type){
        if (type == "ATE"){
            tmp <- pmin(EU_right_upper - pi * (1 - pi), pi * (1 - pi))
            0.5 + tmp - sqrt((1 - 2 * pi)^2 / 4 + tmp^2)
        } else if (type == "ATT"){
            EU_right_upper / (EU_right_upper + pi^2)
        } else if (type == "ATC"){
            EU_left_upper / (EU_left_upper + (1 - pi)^2)
        } 
    }
    return(ovalue_fun)
}

## Inexact DiR O-value
DiR_ovalue_inexact <- function(score1, score0){
    n1 <- length(score1)
    n0 <- length(score0)
    pstar <- sum(fast_rank(score0, score1)) / n1 / n0

    ovalue_fun <- function(pi, type){
        if (type == "ATE"){
            tmp <- pi * (1 - pi) * abs(1 - 2 * pstar)
            0.5 - tmp - sqrt((1 - 2 * pi)^2 / 4 + tmp^2)
        } else if (type == "ATT"){
            numer <- 2 * (1 - pi) * (1 - pstar)
            numer / (numer + pi)
        } else if (type == "ATC"){
            numer <- 2 * pi * pstar
            numer / (numer + 1 - pi)
        } 
    }
    return(ovalue_fun)
}
