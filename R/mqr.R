##' Performs moving quantile regression
##'
##' Performs moving quantile regression
##' @title Moving quantile regression
##' @param formula a formula specifying moving quantile regression
##' @param data data data
##' @param start the starting values for discount function
##' @param ... additional parameters
##' @return an mqr object
##' @author Vaidotas Zemlys
##' @export
mqr <- function(formula, data, start=NULL, Ofunction="lm", ...) {
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf)
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    X <- model.matrix(mt, mf)

    out <- prepmqr(y, X, mt, Zenv, cl, start, Ofunction)
    class(out) <- mqr
    mqr.fit(out)
}

prepmqr <- function(y, X, mt, Zenv, cl, start, Ofunction) {
    
    terms.lhs <- as.list(attr(mt,"variables"))[-2:-1]
    term.labels <- attr(mt,"term.labels") 

    rfd <- vector("list",length(terms.lhs))

    dterm <- function(fr){
        qlev <- eval(fr[[2]], Zenv)
        lagstruct <- eval(fr[[3]], Zenv)
        mf <- fr[c(-5, -2)]
        mf[[1]] <- fr[[5]]
        dcf <- function(p) {
            mf[[2]] <- p
            eval(mf, Zenv)
        }
        return(list(term_name = as.character(fr[[2]]),
                    qlev = qlev,
                    lag_structure = lagstruct,
                    ncol = 1,
                    discount_factor_name = as.character(fr[[5]]),
                    discount_factor = dcf,
                    start = c(0, 0)
                    ))        
    }
    
    for(i in 1:length(rfd)) {
        fr <- terms.lhs[[i]]
        fun <- as.character(fr)[1]
        rfd[[i]] <- dterm(fr,fun)
    }

    if(any(sapply(rfd, "[[", "discount_factor_name") != "")) {
        Ofunction <- "lm"
    }

    build_indices <- function(ci) {
        inds <- cbind(c(1,ci[-length(ci)]+1),ci)
        inds <- apply(inds,1,function(x)list(x[1]:x[2]))
        inds <- lapply(inds,function(x)x[[1]])
        inds
    }

    xinds <- build_indices(cumsum(sapply(rfd, "[[", "ncol")))
    dinds <- {}
    tinds <- {}
        
    
    rfd <- mapply(function(l, x, d, t) c(l, list(xindex = x, dindex = d, tindex = t)),
                      rfd, xinds, dinds, tinds, SIMPLIFY = FALSE)
    
    quantile_discounted<- function(X, p) {
        do.call("cbind",
                lapply(rfd,
                       function(l) {
                           if(l$discount_factor_name == "") return(X[, l$xindex])
                           dcf <- l$discount_factor(p[l$dindex])
                           mq(X[, l$xindex], l$qlev, l$lag_structure, dcf)                   
                       })
                )
    }

    term_coef <- function(p) {
        p[sapply(rfd, "[[", "tindex")]
    }

    rhs <- function(p, ...) {
        Xd <- quantile_discounted(X, p)
        coefs <- term_coef(p)
        X %*% coefs
    }
    
    fn0 <- function(p, ...) {
        r <- y - rhs(p)
        sum(r^2)
    }
    
    out <- list(coefficients = coef(fit),
                model <- cbind(y,X),
                terms = mt,
                call = cl)                
}


mqr.fit <- function(object) {
    if (Ofunction == "lm") {
        if(attr(mt,"intercept") == 1) {
            fit <- lm(y~.,data=data.frame(cbind(y,X[,-1]),check.names=FALSE))
        } else {
              fit <- lm(y~.-1,data=data.frame(cbind(y,X),check.names=FALSE))
          }
        
    }
}
