## Log upper-tail inequalities of variance
hoeffding_var <- function(sigma2, x, n){
    m <- floor(n / 2)
    -m * h1(pmin(sigma2, x), sigma2)
}

bentkus_var <- function(sigma2, x, n){
    m <- floor(n / 2)
    log(pbinom(ceiling(m * x), m, sigma2, lower.tail = TRUE)) + 1
}

maurer_pontil_var <- function(sigma2, x, n){
    -(n - 1) / 2 / sigma2 * pmax(sigma2 - x, 0)^2
}

maurer_var <- function(sigma2, x, n){
    a <- n / (n - 1)
    lambda_fun <- function(lambda){
        Glambda <- exp(lambda) - lambda - 1
        if (is.nan(Glambda)){
            Glambda <- Inf
        }
        lambda * (x - sigma2 / (a * Glambda / lambda + 1))
    }
    left <- 1e-5
    right <- 2
    while (TRUE){
        obj <- optimize(lambda_fun, c(1e-5, 2))
        if (obj$minimum >= 2){
            right <- 10 * right
            left <- right
            obj <- optimize(lambda_fun, c(left, right))
        } else {
            break
        }
    }
    logprob <- n * obj$objective
    logprob
}

## Upper confidence bound for the variance
sigma_upper_bound <- function(x, alpha, n){
    sigma2hat <- pmax(var(x), 1e-10)
    tailprob <- function(sigma2){
        hoeffding_bound <- hoeffding_var(sigma2, sigma2hat, n)
        bentkus_bound <- bentkus_var(sigma2, sigma2hat, n)
        maurer_bound <- maurer_var(sigma2, sigma2hat, n)
        min(hoeffding_bound, bentkus_bound, maurer_bound) - log(alpha)
    }
    if (tailprob(0.25) > 0){
        0.5
    } else {
        sqrt(uniroot(tailprob, c(sigma2hat, 0.25))$root)
    }
}

## WSR upper confidence bound 
WSR_upper <- function(x, alpha){
    n <- length(x)
    muhat <- (cumsum(x) + 0.5) / (1 + 1:n)
    sigma2hat <- (cumsum((x - muhat)^2) + 0.25) / (1 + 1:n)
    sigma2hat <- c(0.25, sigma2hat[-n])
    lambda <- pmin(sqrt(2 * log(1 / alpha) / n / sigma2hat), 1)

    Kn <- function(mu){
        max(cumsum(log(1 - lambda * (x - mu)))) + log(alpha)
    }
    if (Kn(1) < 0){
        1
    } else {
        uniroot(Kn, c(0, 1),
                tol = .Machine$double.eps^0.5)$root
    }
}

## WSR lower confidence bound 
WSR_lower <- function(x, alpha){
    n <- length(x)
    muhat <- (cumsum(x) + 0.5) / (1 + 1:n)
    sigma2hat <- (cumsum((x - muhat)^2) + 0.25) / (1 + 1:n)
    sigma2hat <- c(0.25, sigma2hat[-n])
    lambda <- pmin(sqrt(2 * log(1 / alpha) / n / sigma2hat), 1)

    Kn <- function(mu){
        max(cumsum(log(1 + lambda * (x - mu)))) + log(alpha)
    }
    if (Kn(0) < 0){
        0
    } else {
        uniroot(Kn, c(0, 1),
                tol = .Machine$double.eps^0.5)$root
    }
}

## Exact DiM O-value
DiM_ovalue_exact <- function(score1, score0, delta){
    n1 <- length(score1)
    n0 <- length(score0)
    delta <- delta / 4
    mu1_lower <- WSR_lower(score1, delta)
    mu0_upper <- WSR_upper(score0, delta)
    mu_diff_lower <- max(0, mu1_lower - mu0_upper)    

    ATE_sigma1_upper <- sigma_upper_bound(score1, delta, n1)
    ATE_sigma0_upper <- sigma_upper_bound(score0, delta, n0)
    ATTC_sigma1_upper <- sigma_upper_bound(score1, 2 * delta, n1)
    ATTC_sigma0_upper <- sigma_upper_bound(score0, 2 * delta, n0)
    
    ATE_T_lower <- mu_diff_lower / c(ATE_sigma0_upper, ATE_sigma1_upper)
    ATTC_T_lower <- mu_diff_lower / c(ATTC_sigma0_upper, ATTC_sigma1_upper)
    
    ovalue_fun <- function(pi, type){
        if (type == "ATE"){
            numer <- pi * (1 - pi)
            denom <- 1 + pmax(pi * ATE_T_lower[1], (1 - pi) * ATE_T_lower[2])^2
            0.5 - sqrt(0.25 - numer / denom)
        } else if (type == "ATT"){
            numer <- 1 - pi
            denom <- 1 + pi * ATTC_T_lower[1]^2
            numer / denom
        } else if (type == "ATC"){
            numer <- pi
            denom <- 1 + (1 - pi) * ATTC_T_lower[2]^2
            numer / denom
        } 
        
    }
    return(ovalue_fun)
}

## Inexact DiM O-value
DiM_ovalue_inexact <- function(score1, score0){
    mu1 <- mean(score1)
    mu0 <- mean(score0)
    mu_diff <- abs(mu1 - mu0)
    sigma1 <- sd(score1)
    sigma0 <- sd(score0)
    T1 <- mu_diff / sigma1
    T0 <- mu_diff / sigma0

    ovalue_fun <- function(pi, type){
        if (type == "ATE"){
            numer <- pi * (1 - pi)
            denom <- 1 + pmax(pi * T0, (1 - pi) * T1)^2
            0.5 - sqrt(0.25 - numer / denom)
        } else if (type == "ATT"){
            numer <- 1 - pi
            denom <- 1 + pi * T0^2
            numer / denom
        } else if (type == "ATC"){
            numer <- pi
            denom <- 1 + (1 - pi) * T1^2
            numer / denom
        }
    }
    return(ovalue_fun)    
}
