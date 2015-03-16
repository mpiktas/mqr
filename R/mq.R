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
    X <- mls(x, k)
    apply(X, 1, quantile, probs = alpha)
}
