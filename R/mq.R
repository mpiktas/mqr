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
mq <- function(x, alpha, k, ... ) {    
    X <- mls(x, k, 1)
    r <- apply(X, 1, quantile, probs = alpha, na.rm = TRUE)
    if(length(alpha)==1) {
        r <- matrix(r, nrow = length(r))
    }else {
        r <- t(r)
    }
    r
}
