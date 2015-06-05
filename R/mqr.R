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

    Zenv <- new.env(parent=environment(formula))

    cl <- match.call()
    args <- list(...)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf)
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    X <- model.matrix(mt, mf)

    out <- prepmqr(y, X, mt, Zenv, cl, args, start, Ofunction)
    class(out) <- mqr
    mqr.fit(out)
}

prepmqr <- function(y, X, mt, Zenv, cl, args, start, Ofunction) {
    
    terms.lhs <- as.list(attr(mt,"variables"))[-2:-1]
    term.labels <- attr(mt,"term.labels") 

    rfd <- vector("list",length(terms.lhs))

    dterm <- function(fr) {
        fun <- as.character(fr[[1]])
        tname <- deparse(fr)
        qlev <- dcf <- start <- NULL
        ncol <- 1
        lagstruct <- 0        

        if(fun == "mqr") {
            qlev <- eval(fr[[2]], Zenv)
            lagstruct <- eval(fr[[3]], Zenv)
            tname <- as.character(fr[[2]])
            if(length(fr)>=5) {
                discount <- eval(fr[[5]],Zenv)
                if(is.function(discount)) {
                    mf <- fr[c(-5, -2)]
                    mf[[1]] <- fr[[5]]
                    for(j in 3:length(mf)) mf[[j]] <- eval(mf[[j]], Zenv)
                    dcf <- function(p) {
                        mf[[2]] <- p
                        eval(mf, Zenv)
                    }
                    ncol <- 1
                    start <- 0
                    tname <- paste(tname, qlev*100, as.character(fr[[5]]), sep="_")
                } else {
                      ncol <- length(qlev)
                      tname <- paste(tname, qlev*100, sep="_")
                 }
            }    
        } else {
              if(fun %in% c("mls","dmls","fmls")) {
                  ll <- eval(fr[[3]],Zenv)
                  lagstruct <- switch(fun,
                                      fmls = 0:ll,
                                      dmls = 0:ll,
                                      mls = ll)                  
                  ncol <- length(lagstruct)
                  tname <- paste0(as.character(fr[[2]]), lagstruct)
              }                          
        }
        return(list(term_name = tname,
                    qlev = qlev,
                    lag_structure = lagstruct,
                    ncol = ncol,                    
                    discount_factor = dcf,
                    start <- rep(0, ncol)
                    ))        
    }

    rfd <- lapply(terms.lhs,dterm)
    if(attr(mt, "intercept") == 1) {
        intc <- dterm()
    }
    
    
    dfn <- !sapply(lapply(rfd, "[[", "discount_factor"),is.null)
    discountX <- any(dfn)
    
    if(!discountX) {
        ##We do not have to discount anything
        Ofunction <- "lm"                
    } else {
          terms_with_discount <- which(dfn)
          if(length(start) != length(terms_with_discount)) {
              stop("Please supply starting values for all terms which are discounted")
          }
          for(i in 1:length(terms_with_discount)) {
              rfd[[terms_with_discount[i]]][["start"]] <- c(rep(0,rfd[[terms_with_discount[i]]][["qlev"]]),start[[i]])
          }
      }

    
    indexes <- lapply(rfd, function(l) {
               if(is.null(l$discount_factor)) {
                   list(xdi = length(l$lag_structure),
                        tfi = length(l$lag_structure),
                        tmi = 1:length(l$lag_structure),
                        tdi = 0
                        )
               } else {
                     ##First coefficients in start are multiplication constants
                     st <- l$discount_factor(l$start[-(1:length(l$qlev))])
                     if(length(st)!=length(l$lag_structure)) stop()
                     list(xdi = length(l$qlev),
                          tfi = length(l$qlev) + length(l$start),
                          tmi = 1:length(l$qlev),
                          tdi = length(l$qlev) + 1:length(l$start)
                          )                     
                 }
           })
    
    
    build_indices <- function(ci) {
        inds <- cbind(c(1,ci[-length(ci)]+1),ci)
        inds <- apply(inds,1,function(x)list(x[1]:x[2]))
        inds <- lapply(inds,function(x)x[[1]])
        inds
    }

    xinds <- build_indices(cumsum(sapply(rfd, "[[", "ncol")))
    xdinds <- build_indices(cumsum(sapply(indexes, "[[", "xdi")))
    tfinds <- build_indices(cumsum(sapply(indexes, "[[", "tfi")))
    tdinds <- mapply(function(f, s)f[s], tfinds, lapply(indexes, "[[", "tdi"), SIMPLIFY = FALSE)
    tminds <- mapply(function(f, s)f[s], tfinds, lapply(indexes, "[[", "tdi"), SIMPLIFY = FALSE)
            
    rfd <- mapply(function(l, x, xd, tf, tm, td) c(l, list(index=list(x = x, xd = xd, tf = tf, tm = tm, td = td))),
                      rfd, xinds, xdinds, tf, tminds, tfinds, SIMPLIFY = FALSE)
    
    quantile_discounted<- function(X, p) {
        if(discountX) {
            do.call("cbind",
                lapply(rfd,
                       function(l) {
                           if(is.null(l$discount_factor)) return(X[, l$index$x])
                           dcf <- l$discount_factor(p[l$index$td])
                           mq(X[, l$index$x], l$qlev, l$lag_structure, dcf)                   
                       })
                    )
        }else X
    }

    tm <- sapply(lapply(rfd, "[[", "index"),"[[","tm")
    
    starto <- unlist(lapply(rfd, "[[", "start"))
    names(starto) <- unlist(lapply(rfd, "[[", "term_name"))
    
    if(Ofunction != "lm") {
        Xd <- quantile_discounted(X, p)
        starto[tm] <- coef(lsfit(Xd, y, intercept = FALSE))
    }
        
    rhs <- function(p, ...) {
        Xd <- quantile_discounted(X, p)
        coefs <- p[tm]
        X %*% coefs
    }
    
    fn0 <- function(p, ...) {
        r <- y - rhs(p)
        sum(r^2)
    }

    gr <- function(p)grad(fn0,p)
    hess <- function(x)numDeriv::hessian(fn0,x)
    
    control <- c(list(Ofunction=Ofunction),args)
    ##Override default method of optim. Use BFGS instead of Nelder-Mead
    if(!("method"%in% names(control)) & Ofunction=="optim") {        
        control$method <- "BFGS"
    }
    
    out <- list(coefficients = starto,
                model <- cbind(y,X),
                dX <- function(p) quantile_discounted(X,p),
                term_info = term_info,
                fn0 = fn0,
                rhs = rhs,
                opt = NULL,
                argmap_opt = control,
                start_opt = starto,
                start_list = start,
                terms = mt,
                call = cl,
                gradient = gr,
                hessiand= hess,
                Zenv = Zenv)                
}


mqr.fit <- function(c) {
    args <- x$argmap_opt
    function.opt <- args$Ofunction
    args$Ofunction <- NULL
   
    if(function.opt=="optim" | function.opt=="spg") {  
        args$par <- x$start_opt
        args$fn <- x$fn0
        opt <- try(do.call(function.opt,args),silent=TRUE)
        if(inherits(opt,"try-error")) {
            stop("The optimisation algorithm  failed with the following message:\n", opt,"\nPlease try other starting values or a different optimisation function")
        }
        par <- opt$par
        names(par) <- names(coef(x))
        x$convergence <- opt$convergence
    }
    
    if (function.opt == "lm") {
        fit <- lstfit(x$model[,-1], x$model[,1], intercept = FALSE)
        par <- coef(fit)
        names(par) <- names(coef(x))
        opt <- NULL
        x$convergence <- 0
    }

    if(function.opt=="nls") {
        rhs <- x$rhs
        y <- x$model[,1]
        args$formula <- formula(y~rhs(p))
        args$start <- list(p=x$start_opt)
        opt <- try(do.call("nls",args),silent=TRUE)
        if(inherits(opt,"try-error")) {
            stop("The optimisation algorithm of MIDAS regression failed with the following message:\n", opt,"\nPlease try other starting values or a different optimisation function")
        }
        par <- coef(opt)
        names(par) <- names(coef(x))
        x$convergence <- opt$convInfo$stopCode
    }
    
    x$opt <- opt
    x$coefficients <- par    
    ##Get Xd and get the coefficients corresponding to Xd
    tm <- sapply(lapply(x$term_info, "[[", "index"),"[[","tm")
    Xd <- x$Xd(par)
    x$fitted.values <- as.vector(Xd %*% par)
    x$residuals <- as.vector(x$model[, 1] - x$fitted.values)
    x
}
    
