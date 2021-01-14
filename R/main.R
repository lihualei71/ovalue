#' O-values of Observational Studies
#'
#' \code{ovalue} is a framework to calculate O-values of observational studies allowing for arbitrary classifiers to calculate 1-d scores.
#'
#' @details \code{ovalue} provides a bunch of user-friendly wrappers and a flexible framework that allows for arbitrary classifiers. It also supports using multiple classifiers with each classifier being assigned a fraction of confidence budget. The argument \code{scorefuns} gives either one or multiple classifiers and the argument \code{sw} gives their weights. Specifically, the algorithm will assign \code{sw[i] / sum(sw)} fraction of confidence budget to the i-th method.
#'
#' Each element of \code{scorefuns} can be a valid string, including
#' * "logistic" for logistic regression,
#' * "lasso" for L1-penalized logistic regression,
#' * "gam" for generalized additive model,
#' * "gbm" for generalized boosting machine,
#' * "rf" for random forest,
#'
#' or a function object whose inputs must include
#' * "T" for treatment vector, must be a logical vector or a factor/vector encoded by 0 and 1,
#' * "X" for covariates, must be a vector/matrix/data.frame,
#' * "trainid" for the index of training samples, must be a logical vector or a vector of integers.
#' * "testid" for the index of testing samples, must be a logical vector or a vector of integers. "testid" is allowed to overlap with "trainid".
#' The default setting is \code{scorefuns = c("rf", "gbm"), sw = c(1, 1)}.
#'
#' \code{ovalue} supports two types of data inputs: (1) \code{T} and \code{X} or (2) \code{formula} and \code{data}. One of the pair has to be specified.
#'
#' Similar to the classifiers, \code{ovalue} provides a bunch of testing methods and a flexible framework that allows for user-specified external testing methods. It also supports hybrid version of multiple testing methods with each method assigned a fraction of confidence budget. The argument \code{methods} gives either one or multiple testing methods and the argument \code{mw} gives their weights. Specifically, the algorithm will assign \code{mw[i] / sum(mw)} fraction of confidence budget to the i-th test.
#'
#' Each element of \code{methods} can be a valid string, including
#' * "ROC" for ROC bound,
#' * "EBenn" for \eqn{$\chi^2$} bound based on Empirical Bennett inequality,
#'
#' or a function object whose inputs must include
#' * "T" for treatment vector, must be a logical vector or a factor/vector encoded by 0 and 1,
#' * "score" for the univariate scores, must be a vector with the same length as "T",
#' * "delta" for the confidence level, must be a real number in \eqn{[0, 1]}.
#' The default setting is \code{scorefuns = c("ROC", "EBenn"), sw = c(1, 1)}.
#'
#' Derandomization step is crucial to reduce the external randomness from data splitting. \code{ovalue} calculate O-values for \code{nreps} data splits and report the \code{drq}-th quantile. If the confidence level for each split is \eqn{$\beta$}, it is guaranteed that the coverage is at least \eqn{$\beta / drq$}. To guarantee the coverage in worst case, \code{ovalue} corrects \code{alpha} to \code{alpha} * \code{drq}, if \code{drcorrect = TRUE} by default. In practice, users might awant to avoid this level of correction because each O-value is conservative and the overall coverage will still be guaranteed even without the correction.
#'
#' @md
#'
#' @param T logical vector or a factor/vector encoded by 0 and 1. Treatment assignment
#' @param X matrix or data.frame. Covariate
#' @param formula formula object. See Details
#' @param data data.frame. See Details
#' @param alpha numeric. Confidence level
#' @param type vector of strings. Types of overlap condition to be considered. Currently support "ATE", "ATT" and "ATC"
#' @param scorefuns vector of strings or functions. See Details
#' @param sw vector of non-negative numbers. See Details
#' @param methods vector of strings or functions. See Details
#' @param mw vector of non-negative numbers. See Details
#' @param datasplit logical. Indicate whether data splitting is performed.
#' @param trainprop numeric. Proportion of training samples
#' @param nreps an integer. Number of times for data splitting
#' @param drq numeric. The quantile at which the O-values are reported. See Details
#' @param drcorrect logical. Indicate whether the confidence level needs to be corrected for derandomization. See Details
#' @param verbose logical. Indicate whether relevant information is outputted to the console
#' @param return_list logical. Indicate whether the O-values for each split are returned.
#' @param ... Other arguments passed into \code{scorefuns}
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
#' probs <- 1 / (1 + exp(-X \%*\% beta))
#' T <- stats::runif(n) <= probs
#' data <- data.frame(T = T, X)
#'
#' # Calculate the O-value with inputs \code{T} and \code{X}
#' set.seed(1)
#' ovalue(T, X, scorefuns = "gbm")
#'
#' # Calculate the O-value with inputs \code{formula} and \code{data}
#' set.seed(1)
#' ovalue(formula = T ~ ., data = data, scorefuns = "gbm")
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
                   mw = rep(1, length(methods)),
                   datasplit = FALSE,
                   trainprop = 0.5,
                   nreps = 50,
                   drq = 0.5,
                   drcorrect = FALSE,
                   verbose = TRUE,
                   return_list = FALSE,
                   ...){
    if (!datasplit){
        #cat("Warning: there is no theoretical guarantee on Type-I error control without data splitting. \n")
        nreps <- 1
    }
    if (drcorrect){
        alpha <- alpha * drq
    }

    oldw <- getOption("warn")
    options(warn = -1)
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
    gamma_ci<- gamma_CI(T, delta_gamma)
    eta0 <- gamma_ci$pi[2]
    gamma_grid <- seq(gamma_ci$gamma[1], gamma_ci$gamma[2],
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
        if (datasplit){
            trainid <- sample(n, ntrain)
            testid <- setdiff(1:n, trainid)
        } else {
            trainid <- testid <- 1:n
        }
        if (sum(T[trainid]) == 0 || sum(T[testid]) == 0){
            nfails <- nfails + 1
            next
        }
        for (j in 1:length(scorefuns)){
            scorefun <- scorefuns[[j]]
            score <- scorefun(T, X, trainid, testid, ...)
            Ttest <- T[testid]
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
        options(warn = oldw)
        warning(paste0(nfails, " replicates involve only one class in the training or the testing set."))
    }

    eta <- lapply(eta_list, function(x){
        temp <- lapply(x, function(y){
            quantile(y, drq, na.rm = TRUE)
        })
        min(unlist(temp), eta0)
    })

    if (return_list){
        return(c(eta, list(etalist = eta_list)))
    } else {
        return(eta)
    }
}
