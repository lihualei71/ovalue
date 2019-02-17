eta_hybrid <- function(T, score, delta){
    delta_ROC <- delta / 2
    delta_EBenn <- delta / 2
    eta_ROC_fun <- eta_ROC(T, score, delta_ROC)
    eta_EBenn_fun <- eta_EBenn(T, score, delta_EBenn)
    eta_fun <- function(gamma, type){
        pmin(eta_ROC_fun(gamma, type),
             eta_EBenn_fun(gamma, type))
    }
    return(eta_fun)
}
