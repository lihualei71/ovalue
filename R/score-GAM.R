score_GAM <- function(T, X, trainid,
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
    
    mod <- gam::gam(T ~ ., family = "binomial", data = data[trainid, ])
    score <- predict(mod, newdata = data[-trainid, ], type = "response")
    return(as.numeric(score))
}
