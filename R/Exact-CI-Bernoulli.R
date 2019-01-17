gamma_CI <- function(T, delta){
    if (nlevels(as.factor(T)) != 2){
        stop("The treatment 'T' must be binary")
    }
    if (!is.logical(T)){
        T <- as.logical(T)
    }
    m1 <- sum(T)
    m <- length(T)
    
    pi_ci <- exactci::exactbinomCI(m1, m, conf.level = 1 - delta)
    pi_ci <- as.numeric(pi_ci)
    gamma_ci <- pi_ci / (1 - pi_ci)
    return(gamma_ci)
}
