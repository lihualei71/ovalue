## Dempster lower confidence band of F(x)
dempster <- function(n, a, b){
    fac1 <- lgamma(n + 1) + log(a)
    tmp1 <- 0:(floor(n * (1 - b)))
    tmp2 <- tmp1 / n * (1 - a) / (1 - b) + a
    inds <- which(tmp2 > 0 & tmp2 < 1)
    tmp1 <- tmp1[inds]
    tmp2 <- tmp2[inds]
    fac2 <- -lgamma(n - tmp1 + 1) - lgamma(tmp1 + 1) + (tmp1 - 1) * log(tmp2) + (n - tmp1) * log(1 - tmp2)
    sum(exp(fac1 + fac2))
}

dempster_lower_beta1 <- function(n, alpha, beta0){
    if (beta0 > 1){
        return(0)
    }
    find_beta1 <- function(beta1){
        a <- (beta0 + beta1) / (1 + beta1)
        b <- beta0
        dempster(n, a, b) - alpha
    }
    lower <- -beta0 + 1e-6
    upper <- 1
    if (find_beta1(lower) < 0){
        return(lower)
    }
    while (upper < 1e3){
        if (find_beta1(upper) > 0){
            lower <- upper
            upper <- 10 * upper
        } else {
            break
        }
    }
    if (find_beta1(upper) < 0){
        uniroot(find_beta1, c(lower, upper))$root
    } else {
        upper
    }
}

dempster_lower <- function(x, alpha, turn = 5){
    n <- length(x)
    beta0 <- turn / n
    beta1 <- dempster_lower_beta1(n, alpha, beta0)
    function(y){
        ecdf <- fast_ecdf(x, y)
        pmax((ecdf - beta0) / (1 + beta1), 0)
    }
}

## DKWM lower confidence band of F(x)
DKWM_lower <- function(x, alpha){
    n <- length(x)
    correct <- sqrt(log(1 / alpha) / 2 / n)
    function(y){
        pmax(fast_ecdf(x, y) - correct, 0)
    }
}

## Hybrid upper confidence band of F(x)
hybrid_lower <- function(x, alpha, turn = 5){
    alpha <- alpha / 2
    dempster_fun <- dempster_lower(x, alpha, turn)
    DKWM_fun <- DKWM_lower(x, alpha)
    function(y){
        pmax(dempster_fun(y), DKWM_fun(y))
    }
}

## dempster upper confidence band of F(x)
dempster_upper_beta1 <- function(n, alpha, beta0){
    find_beta1 <- function(beta1){
        a <- beta0 / (1 - beta1)
        b <- beta0 + beta1
        dempster(n, a, b) - alpha
    }
    lower <- -beta0 + 1e-6
    upper <- 1 - beta0 - 1e-6
    if (find_beta1(lower) < 0){
        return(lower)
    }
    if (find_beta1(upper) > 0){
        return(upper)
    }
    uniroot(find_beta1, c(lower, upper))$root
}

dempster_upper <- function(x, alpha, turn = 5){
    n <- length(x)
    beta0 <- turn / n
    beta1 <- dempster_upper_beta1(n, alpha, beta0)
    function(y){
        ecdf <- fast_ecdf(x, y)
        pmin((ecdf + beta0) / (1 - beta1), 1)
    }
}

## Simes upper confidence band of F(x)
simes_upper <- function(x, alpha, kfrac = 0.5){
    x <- sort(x)
    n <- length(x)
    k <- ceiling(kfrac * n)
    fac1 <- log(alpha) / k - mean(log((n - k + 1):n))
    fac2 <- zoo::rollmean(log(1:n), k)
    bseq <- c(1 - exp(rev(fac2 + fac1)), rep(1, k))
    function(y){
        inds <- find_posit_vec(y, x, "right", FALSE)
        bseq[inds]
    }
}

## DKWM upper confidence band of F(x)
DKWM_upper <- function(x, alpha){
    n <- length(x)
    correct <- sqrt(log(1 / alpha) / 2 / n)
    function(y){
        pmin(fast_ecdf(x, y) + correct, 1)
    }
}

## Hybrid upper confidence band of F(x)
hybrid_upper <- function(x, alpha, kfrac = 0.5){
    alpha <- alpha / 2
    ## dempster_fun <- dempster_upper(x, alpha, turn)
    DKWM_fun <- DKWM_upper(x, alpha)    
    simes_fun <- simes_upper(x, alpha, kfrac)
    function(y){
        pmin(DKWM_fun(y), simes_fun(y))
    }
}

## Exact DiT O-value
DiT_ovalue_exact <- function(score1, score0, delta,
                             kfrac = 0.5,
                             turn = 5){
    n1 <- length(score1)
    n0 <- length(score0)
    delta <- delta / 8

    score1_left <- pmin(score1 + 1e-10 * runif(n1), 1)
    score1_right <- pmin(1 - score1 + 1e-10 * runif(n1), 1)
    score0_left <- pmin(score0 + 1e-10 * runif(n0), 1)
    score0_right <- pmin(1 - score0 + 1e-10 * runif(n0), 1)

    ATE_F1_left_upper <- hybrid_upper(score1_left, delta, kfrac)
    ATE_F0_left_upper <- hybrid_upper(score0_left, delta, kfrac)
    ATE_F1_left_lower <- hybrid_lower(score1_left, delta, turn)
    ATE_F0_left_lower <- hybrid_lower(score0_left, delta, turn)
    ATE_F1_right_upper <- hybrid_upper(score1_right, delta, kfrac)
    ATE_F0_right_upper <- hybrid_upper(score0_right, delta, kfrac)
    ATE_F1_right_lower <- hybrid_lower(score1_right, delta, turn)
    ATE_F0_right_lower <- hybrid_lower(score0_right, delta, turn)    

    ATTC_F1_left_upper <- hybrid_upper(score1_left, delta * 2, kfrac)
    ATTC_F0_left_upper <- hybrid_upper(score0_left, delta * 2, kfrac)
    ATTC_F1_left_lower <- hybrid_lower(score1_left, delta * 2, turn)
    ATTC_F0_left_lower <- hybrid_lower(score0_left, delta * 2, turn)
    ATTC_F1_right_upper <- hybrid_upper(score1_right, delta * 2, kfrac)
    ATTC_F0_right_upper <- hybrid_upper(score0_right, delta * 2, kfrac)
    ATTC_F1_right_lower <- hybrid_lower(score1_right, delta * 2, turn)
    ATTC_F0_right_lower <- hybrid_lower(score0_right, delta * 2, turn)    

    ovalue_fun <- function(pi, type){
        x_left <- c(score1_left, score0_left)
        x_right <- c(score1_right, score0_right)        
        if (type == "ATE"){
            nu0_left <- max(ATE_F0_left_lower(x_left) / ATE_F1_left_upper(x_left))
            nu1_left <- max(ATE_F1_left_lower(x_left) / ATE_F0_left_upper(x_left))
            nu0_right <- max(ATE_F0_right_lower(x_right) / ATE_F1_right_upper(x_right))
            nu1_right <- max(ATE_F1_right_lower(x_right) / ATE_F0_right_upper(x_right))
            nu0 <- max(nu0_left, nu0_right)
            nu1 <- max(nu1_left, nu1_right)
            pmin(pi / (pi + (1 - pi) * nu0),
                 (1 - pi) / (1 - pi + pi * nu1))
        } else {
            nu0_left <- max(ATTC_F0_left_lower(x_left) / ATTC_F1_left_upper(x_left))
            nu1_left <- max(ATTC_F1_left_lower(x_left) / ATTC_F0_left_upper(x_left))
            nu0_right <- max(ATTC_F0_right_lower(x_right) / ATTC_F1_right_upper(x_right))
            nu1_right <- max(ATTC_F1_right_lower(x_right) / ATTC_F0_right_upper(x_right))
            nu0 <- max(nu0_left, nu0_right)
            nu1 <- max(nu1_left, nu1_right)
            if (type == "ATT"){
                (1 - pi) / (1 - pi + pi * nu1)
            } else if (type == "ATC"){
                pi / (pi + (1 - pi) * nu0)
            }
        } 
    }
    return(ovalue_fun)
}
