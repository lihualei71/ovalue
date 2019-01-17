score_RandomForest <- function(T, X, trainid, ...){
    T <- clean_format_treat(T)
    if (class(X) == c("numeric")){
        X <- as.matrix(X)
    }
    if (is.logical(trainid)){
        trainid <- which(trainid)
    }
    
    Xtrain <- X[trainid, , drop = FALSE]
    Ttrain <- as.factor(T[trainid])
    if (any(as.numeric(table(Ttrain)) == 0)){
        stop("The training set only includes data from one class. Check if the original data is highly unbalanced or change the random seed for data splitting.")
    }
    Xtest <- X[-trainid, , drop = FALSE]
    
    mod <- randomForest::randomForest(Xtrain, Ttrain)
    score <- predict(mod, newx = Xtest, type = "prob")
    ind <- which(colnames(score) == "1")
    return(as.numeric(score[, ind]))
}
