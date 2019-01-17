score_Boosting <- function(T, X, trainid,
                           n.trees = 100,
                           ...){
    T <- clean_format_treat(T)
    if (class(X) == c("numeric")){
        X <- as.matrix(X)
    }
    if (is.logical(trainid)){
        trainid <- which(trainid)
    }

    data <- data.frame(T = T, X)
    Ttrain <- as.factor(T[trainid])
    if (any(as.numeric(table(Ttrain)) == 0)){
        stop("The training set only includes data from one class. Check if the original data is highly unbalanced or change the random seed for data splitting.")
    }
    
    mod <- gbm::gbm(T ~ ., distribution = "bernoulli",
                    data = data, n.trees = n.trees)
    score <- predict(mod, newdata = data[-trainid, ], type = "response", n.trees = n.trees)
    return(as.numeric(score))
}
