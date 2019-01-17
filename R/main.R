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
#' @param scorefun a string or a function. See Details
#' @param trainprop numeric. Proportion of training samples
#' @param nreps an integer. Number of times for data splitting
#' @param verbose logical. Indicate whether relevant information is outputted to the console
#' @param ... Other arguments passed into \code{scorefun}
#'
#' @return
#' \item{eta}{ median of O-values for all splitted data.}
#' \item{etalist}{ list of O-values for all splitted data.}
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
