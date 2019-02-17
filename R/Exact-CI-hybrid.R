eta_methods_map <- function(method){
    switch(method,
           ROC = eta_ROC,
           EBenn = eta_EBenn)
}

eta_hybrid <- function(methods, weights){
    methods <- lapply(methods, eta_methods_map)
    function(T, score, delta){    
        deltas <- delta * weights / sum(weights)
        eta_fun_list <- mapply(
            function(etafun, delta){
                etafun(T, score, delta)
            }, etafun = methods, delta = deltas)
        eta_fun <- function(gamma, type){
            etas <- lapply(eta_fun_list, function(etafun){
                etafun(gamma, type)
            })
            do.call(pmin, etas)
        }
        return(eta_fun)
    }
}
