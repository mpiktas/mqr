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
##' Returns the deviance of a fitted moving quantile regression object
##' @param object a \code{\link{mqr}} object
##' @param ... currently nothing
##' @return The sum of squared residuals
##' @author Virmantas Kvedaras, Vaidotas Zemlys
##' @rdname deviance.mqr
##' @method deviance mqr
##' @export
deviance.mqr <- function(object,...) {
    sum(residuals(object)^2)
}

##' @import sandwich
##' @export
##' @method estfun mqr
estfun.mqr <- function(x,...) {
    XD <- numDeriv::jacobian(x$rhs,coef(x))
    rval <- as.vector(residuals(x))*XD
    colnames(rval) <- names(coef(x))
    rval
}

##' @export
##' @method vcov mqr
vcov.mqr <- function(object, ...) {
    ##Code stolen from stats:::vcov.nls
    sm <- summary(object)
    sm$cov.unscaled * sm$sigma^2
}

##' @export
##' @method bread mqr
bread.mqr <- function(x,...) {
    sx <- summary(x, vcov.=NULL)
    return(sx$cov.unscaled * as.vector(sum(sx$df[1:2])))
}

##' @export
##' @method summary mqr
summary.mqr <- function(object, vcov.=NULL, df=NULL, prewhite=TRUE, ...) {
    r <- as.vector(residuals(object))
    param <- coef(object)
    pnames <- names(param)
    n <- length(r)
    p <- length(param)
    rdf <- n - p
    resvar <- deviance(object)/rdf
    XD <- numDeriv::jacobian(object$rhs, coef(object))
    R <- qr.R(qr(XD))
    XDtXDinv <- chol2inv(R)
    dimnames(XDtXDinv) <- list(pnames,pnames)
        
    se <- sqrt(diag(XDtXDinv)*resvar)

    f <- as.vector(object$fitted.values)
    mss <- if (attr(object$terms, "intercept")) 
        sum((f - mean(f))^2)
    else sum(f^2)
    rss <- sum(r^2)

    n <- length(r)
    p <- length(coef(object))
    rdf <- n-p
    df.int <- if (attr(object$terms, "intercept")) 1L
    else 0L

    r_squared <- mss/(mss + rss)
    adj_r_squared <- 1 - (1 - r_squared) * ((n - df.int)/rdf)
    
    if(!is.null(vcov.)) {
        set <- try(sqrt(diag(vcov.(object,prewhite=prewhite,...))))
        if(class(set)=="try-error") {
            warning("Unable to compute robust standard errors, using non-robust ones. This is an indication of problems with optimisation, please try other starting values or change optimisation method")
        } else se <- set
    }
    tval <- param/se

    #Code stolen from coeftest function
    if(is.null(df)) {
        df <- rdf
    }
    if (is.finite(df) && df > 0) {
        pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
        cnames <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
        mthd <- "t"
    }
    else {
        pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
        cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
        mthd <- "z"
    }
    
    param <- cbind(param,se,tval,pval)
    dimnames(param) <- list(pnames, c("Estimate", "Std. Error", 
        "t value", "Pr(>|t|)"))
    ans <- list(formula=formula(object$terms),residuals=r,sigma=sqrt(resvar),
                df=c(p,rdf), cov.unscaled=XDtXDinv, call=object$call,
                coefficients=param,
                r_squared = r_squared, adj_r_squared = adj_r_squared)
    class(ans) <- "summary.mqr"
    ans
}

##' @export
##' @method print summary.mqr
print.summary.mqr <- function(x, digits=max(3, getOption("digits") - 3 ), signif.stars = getOption("show.signif.stars"), ...) {
    cat("\n Formula", deparse(formula(x)),"\n")
    df <- x$df
    rdf <- df[2L]
    cat("\n Parameters:\n")
    printCoefmat(x$coefficients,digits=digits,signif.stars=signif.stars,...)
    cat("\n Residual standard error:", format(signif(x$sigma,digits)), "on", rdf , "degrees of freedom\n")
 #   cat(" Multiple R-squared: ", formatC(x$r_squared, digits = digits))
 #   cat(",\tAdjusted R-squared: ", formatC(x$adj_r_squared, digits = digits),"\n")
    invisible(x)
}
