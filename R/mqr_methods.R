##' @export
##' @method print mqr
print.mqr <- function(x, digits=max(3,getOption("digits")-3),...) {
    cat("Moving quantile regression model\n")
    cat(" model:", deparse(formula(x)),"\n")
    print(coef(x),digits = digits, ...)
    invisible(x)
}

