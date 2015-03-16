##' Performs moving quantile regression
##'
##' Performs moving quantile regression
##' @title Moving quantile regression
##' @param formula a formula specifying moving quantile regression
##' @param data data data
##' @param ... additional parameters
##' @return an mqr object
##' @author Vaidotas Zemlys
##' @export
mqr <- function(formula, data, ...) {
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf)
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    X <- model.matrix(mt, mf)

    if(attr(mt,"intercept")==1) {
        fit <- lm(y~.,data=data.frame(cbind(y,X[,-1]),check.names=FALSE))
    } else {
        fit <- lm(y~.-1,data=data.frame(cbind(y,X),check.names=FALSE))
    }
    
    out <- list(coefficients = coef(fit),
                model <- cbind(y,X),
                terms = mt,
                call = cl)                
}
