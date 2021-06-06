## Normalized Vapnik upper confidence bound (obsolete)
normalized_vapnik_tail_up <- function(n, alpha, eta){
    c1 <- log(1 / 4 / pnorm(sqrt(2), lower.tail = FALSE))
    c2 <- 5 * sqrt(2 * pi * exp(1)) * (2 * pnorm(1) - 1)
    tailprob <- function(x){
        params <- expand.grid(gamma = seq(0.001, 0.5, 0.001),
                              np = seq(0.5, 3, 0.1))
        gamma <- params$gamma
        np <- ceiling(n^(params$np))
        kappa <- eta + x^2 / 2 + x * sqrt(x^2 / 4 + eta)
        kappa <- eta + (gamma + n / np) / (1 + n / np) * sqrt(kappa)
        fac1 <- 1 - exp(-np * x^2 / 2 * gamma^2 / (1 + gamma^2 * x^2 / 36 / eta))
        fac2 <- 1 - (sqrt(1 + eta) - sqrt(eta))^2 / (np * x^2 * gamma^2)
        log_denom <- log(pmax(0, fac1, fac2))

        g2 <- n / (1 + n / np) * (1 - 1 / (n + np)) * x^2 / 2 * (1 - gamma)^2 / (1 + (1 - gamma)^2 * x^2 / 36 / kappa)
        log_Delta <- log(2) + log(1 + n / np) + params$np * log(n)
        log_prob_greene_wellner <- safe_min(log_Delta - g2 - log_denom)

        tmp <- sqrt(n * (1 + eta) / 2) * (1 - gamma) * x
        gauss <- pnorm(tmp, lower.tail = FALSE)
        extra_term <- log(4 * n) - log_denom
        log_prob_bentkus_dzindzalieta <- safe_min(c1 + log(gauss) + extra_term)
        log_prob_pinelis <- safe_min(log(gauss + dnorm(tmp) / (9 + tmp^2) * c2) + extra_term)
        log_prob_hoeffding <- safe_min(-tmp^2 / 2 + extra_term)

        log_prob <- min(log_prob_greene_wellner,
                        log_prob_bentkus_dzindzalieta,
                        log_prob_pinelis,
                        log_prob_hoeffding)
        log_prob - log(alpha)
    }
    uniroot(tailprob, c(1 / n, 1))$root
}

normalized_vapnik_CI_up <- function(n, muhat, alpha){
    bound <- sapply(muhat, function(mu){
        thr <- normalized_vapnik_tail_up(n, alpha, 0)
        mu + thr^2 / 2 + thr * sqrt(thr^2 / 4 + mu)
    })
    pmin(bound, 1)
}

## Improved Devroye upper confidence bound (obsolete)
improved_devroye <- function(n, alpha){
    tailprob <- function(x){
        params <- expand.grid(gamma = seq(0.001, 0.5, 0.001),
                              np = seq(0.5, 3, 0.1))
        gamma <- params$gamma
        np <- ceiling(n^(params$np))

        tmp <- floor(gamma * x * np) / np
        g3 <- pmax(2 * np * (gamma * x)^2,
                   2 * np * tmp^2 +
                   log(8 * pi) / 2 +
                   log(np) +
                   log(tmp) +
                   log((1 - 2 * tmp) / (1 + 2 * tmp)) / 2)
        log_denom <- log(1 - exp(-g3))

        g4 <- 2 * n / (1 + n / np) / (1 + 1 / n) * x^2 * (1 - gamma)^2
        log_Delta <- log(1 + (n + 1) / np) + params$np * log(n)

        log_prob <- safe_min(log_Delta - g4 - log_denom)
        log_prob - log(alpha)
    }
    uniroot(tailprob, c(0, 1))$root
}

## Two obsolete CE O-values
CE_ovalue_exact_nvapnik <- function(score1, score0, delta,
                                    kfrac = 0.5,
                                    turn = 5){
    n1 <- length(score1)
    n0 <- length(score0)
    n <- n1 + n0

    cutoff_list <- c(score1, score0)
    err <- (fast_rank(score1, cutoff_list) +
            fast_rank(1 - score0, 1 - cutoff_list)) / n
    err <- min(err)
    err_upper <- normalized_vapnik_CI_up(n, err, delta)
    
    ovalue_fun <- function(pi, type){
        if (type == "ATE"){
            err_upper
        } else if (type == "ATT"){
            NA
        } else if (type == "ATC"){
            NA
        } 
    }
    return(ovalue_fun)
}

CE_ovalue_exact_devroye <- function(score1, score0, delta,
                                    kfrac = 0.5,
                                    turn = 5){
    n1 <- length(score1)
    n0 <- length(score0)
    n <- n1 + n0

    cutoff_list <- c(score1, score0)
    err <- (fast_rank(score1, cutoff_list) +
            fast_rank(1 - score0, 1 - cutoff_list)) / n
    err <- min(err)
    err_upper <- err + improved_devroye(n, delta)
    
    ovalue_fun <- function(pi, type){
        if (type == "ATE"){
            err_upper
        } else if (type == "ATT"){
            NA
        } else if (type == "ATC"){
            NA
        } 
    }
    return(ovalue_fun)
}

## Exact CE O-value
CE_ovalue_exact <- function(score1, score0, delta,
                            kfrac = 0.5){
    n1 <- length(score1)
    n0 <- length(score0)
    n <- n1 + n0
    delta <- delta / 2

    score1 <- pmin(score1 + 1e-10 * runif(n1), 1)
    score0 <- pmin(score0 + 1e-10 * runif(n0), 1)

    class_err1 <- hybrid_upper(score1, delta, kfrac)
    class_err0 <- hybrid_upper(1 - score0, delta, kfrac)
    cutoff_list <- c(score1, score0)
    
    ovalue_fun <- function(pi, type){
        if (type == "ATE"){
            err <- pi * class_err1(cutoff_list) + (1 - pi) * class_err0(1 - cutoff_list)
            min(err)
        } else if (type == "ATT"){
            NA
        } else if (type == "ATC"){
            NA
        } 
    }
    return(ovalue_fun)
}

## Inexact CE O-value
CE_ovalue_inexact <- function(score1, score0){
    cutoff_list <- c(score1, score0)    
    class_err1 <- fast_ecdf(score1, cutoff_list)
    class_err0 <- fast_ecdf(1 - score0, 1 - cutoff_list)
    
    ovalue_fun <- function(pi, type){
        if (type == "ATE"){
            err <- pi * class_err1 + (1 - pi) * class_err0
            min(err)
        } else if (type == "ATT"){
            NA
        } else if (type == "ATC"){
            NA
        } 
    }
    return(ovalue_fun)
}
