eta_ROC <- function(T, score, delta){
    T <- clean_format_treat(T)
    score1 <- score[T]
    score0 <- score[!T]
    m1 <- sum(T)
    m0 <- sum(!T)
    delta1 <- delta / 2
    delta0 <- delta / 2

    F1 <- ecdf(score1)
    F0 <- ecdf(score0)
    adj1 <- sqrt(log(2 / delta1) / 2 / m1)
    adj0 <- sqrt(log(2 / delta0) / 2 / m0)    
    F1plus <- function(x){min(F1(x) + adj1, 1)}
    F1minus <- function(x){max(F1(x) - adj1, 0)}
    F0plus <- function(x){min(F0(x) + adj0, 1)}
    F0minus <- function(x){max(F0(x) - adj0, 0)}    
    nu_minus1 <- suppressWarnings(optimize(
        function(x){F1plus(x) / F0minus(x)},
        c(0, 1))$objective)
    nu_minus2 <- suppressWarnings(optimize(
        function(x){(1 - F1minus(x)) / (1 - F0plus(x))},
        c(0, 1))$objective)
    nu_minus <- min(nu_minus1, nu_minus2)

    nu_plus1 <- suppressWarnings(optimize(
        function(x){F1minus(x) / F0plus(x)},
        c(0, 1), maximum = TRUE)$objective)
    nu_plus2 <- suppressWarnings(optimize(
        function(x){(1 - F1plus(x)) / (1 - F0minus(x))},
        c(0, 1), maximum = TRUE)$objective)
    nu_plus <- max(nu_plus1, nu_plus2)

    eta_fun <- function(gamma, type){
        eta_ATC <- 1 - 1 / (1 + gamma * nu_minus)
        eta_ATT <- 1 / (1 + gamma * nu_plus)
        eta_ATE <- pmin(eta_ATC, eta_ATT)
        eta <- switch(type,
                      ATE = eta_ATE,
                      ATT = eta_ATT,
                      ATC = eta_ATC)
        return(eta)
    }
    return(eta_fun)
}
