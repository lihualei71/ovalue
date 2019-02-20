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

