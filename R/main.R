#' @export
ovalue <- function(T = NULL, X = NULL,
                   formula = NULL, data = NULL,
                   alpha = 0.05,
                   scorefun = "rf",
                   trainprop = 0.5,
                   nreps = 50,
                   verbose = TRUE,
                   ...){
    eta_gen_fun <- eta_hybrid
    scorefun <- clean_format_scorefun(scorefun)

    if (is.null(T)){
        if (is.null(formula)){
            stop("Either 'T' or 'formula' should be specified")
        }
        if (is.null(data)){
            stop("'data' is not specified")
        }
        Tname <- as.character(formula[[2]])
        T <- data[[Tname]]
        X <- model.matrix(formula, data = data)
    }
    if (is.null(X)){
        stop("'X' is not specified")
    }
    n <- length(T)
    ntrain <- ceiling(n * trainprop)
    delta_gamma <- alpha / 10
    delta_others <- alpha - delta_gamma
    gamma_ci <- gamma_CI(T, delta_gamma)
    gamma_grid <- seq(gamma_ci[1], gamma_ci[2],
                      length.out = 1000)

    eta_list <- rep(NA, nreps)
    nfails <- 0
    if (verbose){
        cat("Fitting the scores\n")
        pb <- txtProgressBar(min = 0, max = nreps, style = 3, width = 50)
    }
    for (i in 1:nreps){
        trainid <- sample(n, ntrain)
        if (sum(T[trainid]) %in% c(0, n)){
            nfails <- nfails + 1
        }
        score <- scorefun(T, X, trainid)
        Ttest <- T[-trainid]
        eta_fun <- eta_gen_fun(Ttest, score, delta_others)
        eta_list[i] <- max(eta_fun(gamma_grid))
        if (verbose){
            setTxtProgressBar(pb, i)
        }
    }

    if (verbose){    
        cat("\n")
    }
    if (nfails > 0){
        warning(paste0(nfails, " replicates involve only one class in the training set."))
    }
    
    return(list(eta = median(eta_list),
                etalist = eta_list))
}
