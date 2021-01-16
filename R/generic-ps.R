#' Propensity Score Estimation
#'
#' \code{propensityScore} is a generic wrapper to estimate propensity scores for binary treatment with or without cross-fitting.
#'
#' @details \code{propensityScore} estimates propensity score for binary treatment. It supports both in-sample fitting, where the model is used on the whole dataset and the propensity scores are fitting based on the model, and cross fitting, where the data is splitted into \code{nfolds} subsets and the propensity score for each subset is estimated based on the model fitted on other subsets. In practice, it is crucial to apply cross fitting because the in-sample fits may be pulled towards 0 and 1 because of the non-negligible finite-sample errors.
#'
#' Users can avoid using cross fitting by \code{nfolds = 1}.  Otherwise users can specify the splitting scheme by one of the following argument:
#' * \code{foldid}: a list of indices that form a partition of 1:n. If \code{foldid} is not NULL, then the arguments \code{nfolds} and \code{foldsizes} will be ignored,
#' * \code{foldsizes}: a vector of positive integers that sum up to n. \code{propensityScore} will generate a uniformly random partition with sizes given by the elements of \code{foldsizes},
#' * \code{nfolds}: a positive integer. \code{propensityScore} will generate a uniformly random partition with \code{nfolds} folds and equal sizes (up to 1) within each fold.
#'
#'
#'
#' \code{propensityScore} supports a bunch of algorithms and provides a flexible framework that allows for arbitrary classifiers. The algorithm to use is specified by the argument \code{algo}. \code{algo} can be a valid string, including
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
#' The default setting is \code{algo = "gbm"}.
#'
#' \code{propensityScore} supports two types of data inputs: (1) \code{T} and \code{X} or (2) \code{formula} and \code{data}. One of the pair has to be specified.
#'
#' @md
#'
#' @param T logical vector or factor/vector encoded by 0 and 1. Treatment assignment
#' @param X matrix or data.frame. Covariate
#' @param formula formula object. See Details
#' @param data data.frame. See Details
#' @param nfolds positive integer. Number of folds for cross fitting. See Details
#' @param foldsizes vector of positive integers. See Details
#' @param foldid list of vectors. See Details
#' @param algo vector of strings or functions. See Details
#' @param ... Other arguments passed into \code{algo}
#'
#' @return
#' A vector of propensity scores
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
#' # Calculate propensity scores using GBM with 5-folds cross fitting
#' propensityScore(T, X, algo = "gbm", nfolds = 5)
#'
#' # Calculate propensity scores using GBM without cross fitting
#' propensityScore(T, X, algo = "gbm", nfolds = 1)
#' }
#'
#' @export
propensityScore <- function(T = NULL, X = NULL,
                            formula = NULL, data = NULL,
                            nfolds = 1,
                            foldsizes = NULL,
                            foldid = NULL,
                            algo = c("gbm", "rf", "gam",
                                "lasso", "logistic"),
                            ...){
    algo <- algo[1]
    if (is.character(algo)){
        algo <- switch(
            algo,
            gbm = score_Boosting,
            rf = score_RandomForest,
            gam = score_GAM,
            lasso = score_LassoLogisticReg,
            logistic = score_LogisticReg,
            NULL)
    }
    if (!is.function(algo)){
        stop("algo must be a function or a valid string")
    }
    if (any(!c("T", "X", "trainid", "testid") %in% formalArgs(algo))){
        stop("'algo' must be a valid string or a function with inputs (at least) 'T', 'X', 'trainid' and 'testid'")
    }

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

    if (nfolds == 1){
        return(algo(T, X, 1:n, 1:n, ...))
    }

    if (!is.null(foldid)){
        if (!is.list(foldid)){
            stop("'foldid' must be a list.")
        }
        if (!all(sort(unlist(foldid)) == 1:n)){
            stop("the elements in 'foldid' must be a partition of 1 to n.")
        }
        nfolds <- length(foldid)
    } else {
        if (is.null(foldsizes)){
            quotient <- floor(n / nfolds)
            residual <- n - nfolds * quotient
            foldsizes <- c(rep(quotient + 1, residual),
                           rep(quotient, nfolds - residual))
        }
        if (sum(foldsizes) != n){
            stop("'foldsizes' must sum up to n.")
        }
        nfolds <- length(foldsizes)
        tmpinds <- sample(n, n)
        tmplocs <- cumsum(foldsizes)
        foldid <- lapply(1:nfolds, function(i){
            start <- ifelse(i == 1, 1, tmplocs[i - 1] + 1)
            end <- tmplocs[i]
            tmpinds[start:end]
        })
    }

    ps <- rep(NA, n)
    for (i in 1:nfolds){
        trainid <- unlist(foldid[-i])
        testid <- foldid[[i]]
        ps[testid] <- algo(T, X, trainid, testid, ...)
    }
    return(ps)
}
