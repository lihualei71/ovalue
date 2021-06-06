clean_format_treat <- function(T){
    if (is.logical(T)){
        T <- as.numeric(T)
    } else if (!is.numeric(T) && !is.factor(T)){
        stop("T must be of type numeric/logical/factor")
    } 
    if (nlevels(as.factor(T) != 2)){
        stop("The treatment 'T' must be binary")
    }
    if (!0 %in% T || !1 %in% T){
        stop("T must be encoded as 0 (for the controls) and 1 (for the treated)")
    }
    return(T)
}

clean_format_scorefun <- function(scorefun){
    if (is.character(scorefun)){
        scorefun <- switch(
            scorefun,
            rf = score_RandomForest,
            logistic = score_LogisticReg,
            lasso = score_LassoLogisticReg,
            gbm = score_Boosting,
            gam = score_GAM,
            NA)
    }
    if (!is.function(scorefun)){
        stop("'scorefun' must be a valid string or a function")
    }
    if (any(!c("T", "X", "trainid", "testid") %in% formalArgs(scorefun))){
        stop("'scorefun' must be a valid string or a function with inputs (at least) 'T', 'X', 'trainid' and 'testid'")
    }
    return(scorefun)
}

ovalue_method_map <- function(method,
                              exact = TRUE){
    if (is.function(method)){
        method
    } else if (exact){
        switch(method,
               DiM = DiM_ovalue_exact,
               DiT = DiT_ovalue_exact,
               DiR = DiR_ovalue_exact,
               CE = CE_ovalue_exact)
    } else if (!exact){
        switch(method,
               DiM = DiM_ovalue_inexact,
               DiT = function(score1, score0){
                   DiT_ovalue_exact(score1, score0, 0.05)
               },
               DiR = DiR_ovalue_inexact,
               CE = CE_ovalue_inexact)
    }
}

ovalue_method_hybrid <- function(methods, weights,
                                 exact = TRUE){
    if (is.function(methods)){
        methods <- list(methods)
    }
    ovalue_funs <- lapply(methods,
                          function(method){
                              ovalue_method_map(method, exact)
                          })
    if (exact){
        function(score1, score0, delta){    
            deltas <- delta * weights / sum(weights)
            ovalue_fun_list <- mapply(
                function(ovalue_fun, delta){
                    ovalue_fun(score1, score0, delta)
                }, ovalue_fun = ovalue_funs, delta = deltas)
            ovalue_fun <- function(pi, type){
                ovalues <- lapply(ovalue_fun_list, function(ovalue_fun){
                    ovalue_fun(pi, type)
                })
                do.call(pmin, ovalues)
            }
            return(ovalue_fun)
        }
    } else {
        function(score1, score0){    
            ovalue_fun <- function(pi, type){
                ovalues <- lapply(ovalue_funs, function(ovalue_fun){
                    ovalue_fun <- ovalue_fun(score1, score0)
                    ovalue_fun(pi, type)
                })
                do.call(pmin, ovalues)
            }
            return(ovalue_fun)
        }
    }
}

safe_min <- function(x){
    if (any(is.na(x))){
        -Inf
    } else {
        min(x)
    }
}

h1 <- function(y, mu){
    y * log(y / mu) + (1 - y) * log((1 - y) / (1 - mu))
}

find_posit_vec <- function(vec, target, dir, decreasing = TRUE){
    if (decreasing){
        posit <- rank(-c(vec, target))[1:length(vec)] - rank(-vec)
    } else {
        posit <- rank(c(vec, target))[1:length(vec)] - rank(vec)
    }
    if (dir == "right"){
        posit <- posit + 1
    }
    return(posit)
}

fast_ecdf <- function(x, testx){
    find_posit_vec(testx, x, "left", FALSE) / length(x)
}

fast_rank <- function(x, testx){
    find_posit_vec(testx, x, "left", FALSE)
}
