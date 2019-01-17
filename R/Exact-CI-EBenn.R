## link function between T(gamma) and eta(gamma)
link_Tgamma <- function(Tfun){
    function(gamma){
        denom <- Tfun(gamma)^2 + (sqrt(gamma) + 1 / sqrt(gamma))^2
        1 / 2 - sqrt(1 / 4 - 1 / denom)
    }
}

## h function and inverse of h
h <- function(y){
    (1 + y) * log(1 + y) - y
}

hinv <- function(z){
    ## h(y) <= y^2; h(y) >= y + 2 for any y >= e^2 - 1
    interval <- c(sqrt(z), max(z - 2, exp(2) - 1))
    uniroot(function(y){h(y) - z}, interval)$root
}
hinv <- Vectorize(hinv)

## Solve a univariate equation with multiple roots.
uniroot2 <- function(fun, val, grid, type = "lower"){
    if (type == "lower"){
        tmp <- which(fun(grid) >= val)
        if (length(tmp) == 0) return(NA)
        ind <- min(tmp)
        return(grid[ind])
    } else if (type == "upper"){
        tmp <- which(fun(grid) <= val)
        if (length(tmp) == 0) return(NA)
        ind <- max(tmp)
        return(grid[ind])        
    }
}

## Empirical Bennett's bound
eta_EBenn <- function(T, score, delta){
    T <- clean_format_treat(T)    
    score1 <- score[T]
    score0 <- score[!T]
    m1 <- sum(T)
    m0 <- sum(!T)
    delta_mu1 <- delta / 4
    delta_sigma1 <- delta / 4    
    delta_mu0 <- delta / 4
    delta_sigma0 <- delta / 4
    
    mu1 <- mean(score1)
    mu0 <- mean(score0)
    sigma1 <- sd(score1)
    sigma0 <- sd(score0)

    temp1 <- log(1 / delta_sigma1) / 2 / (m1 - 1)
    sigma1_hat <- sqrt(sigma1^2 + temp1) + sqrt(temp1)
    temp0 <- log(1 / delta_sigma0) / 2 / (m0 - 1)
    sigma0_hat <- sqrt(sigma0^2 + temp0) + sqrt(temp0)

    ## Lemma A.1. implies that
    ##    mu <= RHS of (27) <= mu + sigma^2 h^{-1}(log(1/delta)/m/sigma^2)
    temp1 <- log(1 / delta_mu1) / m1 / sigma1_hat^2
    GL <- function(mu){
        mu + sigma1_hat^2 / (1 - mu) * hinv((1 - mu)^2 * temp1)
    }    
    mu1_hat_upper <- mu1
    mu1_hat_lower <- mu1 - sigma1_hat^2 * hinv(temp1)
    mu1_grid <- seq(mu1_hat_lower, mu1_hat_upper,
                    length.out = 1000)
    mu1_hat <- uniroot2(GL, mu1, mu1_grid, type = "lower")
    if (is.na(mu1_hat)){
        mu1_hat <- mu1_hat_lower
    }

    ## Lemma A.1. implies that
    ##    mu >= RHS of (28) >= mu - sigma^2 h^{-1}(log(1/delta)/m/sigma^2)    
    temp0 <- log(1 / delta_mu0) / m0 / sigma0_hat^2
    GR <- function(mu){
        mu - sigma0_hat^2 / mu * hinv(mu^2 * temp0)
    }    
    mu0_hat_lower <- mu0
    mu0_hat_upper <- mu0 + sigma0_hat^2 * hinv(temp0)
    mu0_grid <- seq(mu0_hat_lower, mu0_hat_upper,
                    length.out = 1000)
    mu0_hat <- uniroot2(GR, mu0, mu0_grid, type = "upper")
    if (is.na(mu0_hat)){
        mu0_hat <- mu0_hat_upper
    }
    
    if (mu1_hat <= mu0_hat){
        TBenn <- function(gamma){0}
    } else {
        TBenn <- function(gamma){
            term1 <- pmax(sqrt(gamma) / sigma0_hat,
                          1 / sqrt(gamma) / sigma1_hat)
            term2 <- mu1_hat - mu0_hat
            term1 * term2
        }
    }
    eta_fun <- link_Tgamma(TBenn)
    return(eta_fun)
}
