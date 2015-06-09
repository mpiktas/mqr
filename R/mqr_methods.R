##' @export
##' @method print mqr
print.mqr <- function(x, digits=max(3,getOption("digits")-3),...) {
    cat("Moving quantile regression model\n")
    cat(" model:", deparse(formula(x)),"\n")
    print(coef(x),digits = digits, ...)
    invisible(x)
}

##' Moving quantile regression model deviance
##'
##' Returns the deviance of a fitted MIDAS regression object
##' @param object a \code{\link{midas_r}} object
##' @param ... currently nothing
##' @return The sum of squared residuals
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @rdname deviance.mqr
##' @method deviance mqr
##' @export
deviance.mqr <- function(object,...) {
    sum(residuals(object)^2,na.rm=TRUE)
}

##' @import sandwich
##' @export
##' @method estfun midas_r
estfun.midas_r <- function(x,...) {
    XD <- numDeriv::grad(x$rhs,coef(x))
    rval <- as.vector(residuals(x))*XD
    colnames(rval) <- names(coef(x))
    rval
}
