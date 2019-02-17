#' O-values of Observational Studies
#'
#' \code{ovalue} is a framework to calculate O-values of observational studies allowing for arbitrary classifiers to calculate 1-d scores.
#'
#' @details \code{ovalue} provides a bunch of user-friendly wrappers and a flexible framework that allows for arbitrary classifiers. For the former, \code{scorefun} can be a valid string, including
#' * "logistic" for logistic regression,
#' * "lasso" for L1-penalized logistic regression,
#' * "gam" for generalized additive model,
#' * "gbm" for generalized boosting machine,
#' * "rf" for random forest.
#' 
#' For the latter, \code{scorefun} can be a function object whose inputs must include
#' * "T" for treatment vector, must be a logical vector or a factor/vector encoded by 0 and 1,
#' * "X" for covariates, must be a vector/matrix/data.frame,
#' * "trainid" for the index of training samples, must be a logical vector or a vector of integers.
#'
#' \code{ovalue} supports two types of data inputs: (1) \code{T} and \code{X} or (2) \code{formula} and \code{data}. One of the pair has to be specified.  
#' @md
#'
#' @param T a logical vector or a factor/vector encoded by 0 and 1. Treatment assignment
#' @param X a matrix or a data.frame. Covariate
#' @param formula a formula object. See Details
#' @param data a data.frame. See Details
#' @param alpha numeric. Confidence level
#' @param type a vector of strings. Types of overlap condition to be considered. Currently support "ATE", "ATT" and "ATC"
#' @param scorefun a string or a function. See Details
#' @param trainprop numeric. Proportion of training samples
#' @param nreps an integer. Number of times for data splitting
#' @param verbose logical. Indicate whether relevant information is outputted to the console
#' @param return_list logical. Indicate whether the O-values for each split are returned.
#' @param ... Other arguments passed into \code{scorefun}
#'
#' @return
#' \item{ATE/ATT/ATC}{ median of O-values for all splitted data for ATE/ATT/ATC.}
#' \item{etalist}{ optional (only returned when \code{return_list = TRUE}). List of O-values for each splitted data. Each entry corresponds to a type of overlap condition.}
#' 
#' @examples
#' \donttest{# Generate data from a logistic model
#' set.seed(1)
#' n <- 1000
#' p <- 50
#' X <- matrix(stats::rnorm(n * p), n, p)
#' beta <- rep(1 / sqrt(p), p)
#' probs <- 1 / (1 + exp(-X %*% beta))
#' T <- stats::runif(n) <= probs
#' data <- data.frame(T = T, X)
#'
#' # Calculate the O-value with inputs \code{T} and \code{X}
#' set.seed(1)
#' ovalue(T, X, scorefun = "gbm")
#'
#' # Calculate the O-value with inputs \code{formula} and \code{data}
#' set.seed(1)
#' ovalue(formula = T ~ ., data = data, scorefun = "gbm")
#' }
#' 
#' @export
ovalue <- function(T = NULL, X = NULL,
                   formula = NULL, data = NULL,
                   alpha = 0.05,
                   type = c("ATE", "ATT", "ATC"),
                   scorefuns = c("rf", "gbm"),
                   sw = rep(1, length(scorefuns)),
                   methods = c("ROC", "EBenn"),
                   mw = rep(1, length(scorefuns)),
                   trainprop = 0.5,
                   nreps = 50,
                   verbose = TRUE,
                   return_list = FALSE,
                   ...){
    eta_gen_fun <- eta_hybrid(methods, mw)
    scorefuns <- lapply(scorefuns, clean_format_scorefun)

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
    delta_others <- (alpha - delta_gamma) * sw / sum(sw)
    gamma_ci <- gamma_CI(T, delta_gamma)
    gamma_grid <- seq(gamma_ci[1], gamma_ci[2],
                      length.out = 1000)

    eta_list <- list()
    for (tp in type){
        for (j in 1:length(scorefuns)){
            eta_list[[tp]][[j]] <- rep(NA, nreps)
        }
    }
    nfails <- 0
    if (verbose){
        cat("Computing O-values for each data splits\n")
        pb <- txtProgressBar(min = 0, max = nreps, style = 3, width = 50)
    }
    for (i in 1:nreps){
        trainid <- sample(n, ntrain)
        if (sum(T[trainid]) %in% c(0, n)){
            nfails <- nfails + 1
        }
        for (j in 1:length(scorefuns)){
            scorefun <- scorefuns[[j]]
            score <- scorefun(T, X, trainid)
            Ttest <- T[-trainid]
            eta_fun <- eta_gen_fun(Ttest, score, delta_others[j])
            for (tp in type){
                eta_list[[tp]][[j]][i] <- max(eta_fun(gamma_grid, tp))
            }
        }
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

    eta <- lapply(eta_list, function(x){
        temp <- lapply(x, median)
        min(unlist(temp))
    })
    if (return_list){
        return(c(eta, list(etalist = eta_list)))
    } else {
        return(eta)
    }
}
