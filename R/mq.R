##' Calculates a moving quantile for a given band
##'
##' Given data take a moving band and calculate a quantile based on that band.
##' @title Moving quantile for a given band
##' @param x a vector
##' @param alpha a quantile
##' @param k band structure
##' @param ... additional parameters
##' @return a matrix
##' @author Vaidotas Zemlys
##' @export
##' @importFrom midasr mls
mq <- function(x, alpha, k, d = NULL, ... ) {    
    X <- mls(x, k, 1)
    if(is.function(d)) return(x)
    else {
        if(!is.null(d)) {
            if(!is.numeric(d)) warning("Discount factors must be a numeric vector. Coercing to numeric")
            d <- as.numeric(d)
            if(length(d)!=ncol(X))stop("The number of discount factors must be the same as the number of lags")
            X <- sweep(X,2,d,"*")        
        }        
        r <- apply(X, 1, quantile, probs = alpha, na.rm = TRUE)
        if(length(alpha)==1) {
            r <- matrix(r, nrow = length(r))
        }else {
             r <- t(r)
         }
        r
    }
}

