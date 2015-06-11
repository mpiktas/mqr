#' @docType package
#' @name mqr
#' @useDynLib mqr
#' @importFrom Rcpp cppFunction 
NULL

##' Performs moving quantile regression
##'
##' Performs moving quantile regression
##' @title Moving quantile regression
##' @param formula a formula specifying moving quantile regression
##' @param data data data
##' @param start the starting values for discount function
##' @param Ofunction the name of optimisation function to use for NLS problem
##' @param ... additional parameters
##' @return an mqr object
##' @import numDeriv
##' @import optimx
##' @author Vaidotas Zemlys
##' @export
mqr <- function(formula, data, start=NULL, Ofunction="optim", ...) {

    Zenv <- new.env(parent=environment(formula))

    cl <- match.call()
    args <- list(...)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    ##We need to do na.action=na.pass here!
    mf[[1L]] <- as.name("model.frame") 
    mf[[4L]] <- as.name("na.pass")
    names(mf)[4] <- c("na.action")
    mf <- eval(mf, Zenv)
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    X <- model.matrix(mt, mf)
    out <- prepmqr(y, X, mt, Zenv, cl, args, start, Ofunction, mf)
    class(out) <- "mqr"
    mqr.fit(out)
}

prepmqr <- function(y, X, mt, Zenv, cl, args, start, Ofunction, mf) {
    
    terms.lhs <- as.list(attr(mt,"variables"))[-2:-1]
    term.labels <- attr(mt,"term.labels") 

    rfd <- vector("list",length(terms.lhs))

    dterm <- function(fr) {
        fun <- as.character(fr[[1]])
        tname <- deparse(fr)
        qlev <- dcf <- start <- NULL
        ncol <- 1
        lagstruct <- 0        

        if(fun == "mq") {
            qlev <- eval(fr[[3]], Zenv)
            lagstruct <- eval(fr[[4]], Zenv)
            tname <- paste0(as.character(fr[[2]]),"_q",paste(qlev*100,sep="_"),"_l",deparse(lagstruct))
            ncol <- length(qlev)
            if(length(fr)>=5) {
                discount <- eval(fr[[5]],Zenv)
                tname <- paste0(tname,"_d",deparse(fr[[5]]))
                if(is.function(discount)) {
                    mf <- fr[c(-5, -2)]
                    mf[[1]] <- fr[[5]]
                    for(j in 3:length(mf)) mf[[j]] <- eval(mf[[j]], Zenv)
                    mf[[3]] <- length(mf[[3]])
                    dcf <- function(p) {
                        mf[[2]] <- eval(p)                        
                        eval(mf, Zenv)
                    }
                    ncol <- 1
                    start <- 0
                    tname <- paste(tname, qlev*100, as.character(fr[[5]]), sep="_")
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
                    start = rep(0, ncol)
                    ))        
    }

    rfd <- lapply(terms.lhs,dterm)
    if(attr(mt, "intercept") == 1) {
        intc <- dterm(expression(1))
        intc$term_name <- "(Intercept)"
        rfd <- c(list(intc),rfd)
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
              twd <- terms_with_discount[[i]]
              rfd[[twd]][["start"]] <- c(rep(0,length(rfd[[twd]][["qlev"]])),start[[i]])
              rfd[[twd]][["term_name"]] <- c(rfd[[twd]][["term_name"]],paste0("df_",i,"_",1:length(start[[i]])))
          }
      }

    indexes <- lapply(rfd, function(l) {
                          if(is.null(l$discount_factor)) {                              
                   list(xdi = l$ncol,
                        tfi = l$ncol,
                        tmi = 1:l$ncol,
                        tdi = 0
                        )
               } else {
                     ##First coefficients in start are multiplication constants
                     st <- l$discount_factor(l$start[-(1:length(l$qlev))])
                     if(length(st)!=length(l$lag_structure)) stop()
                     list(xdi = length(l$qlev),
                          tfi = length(l$start),
                          tmi = 1:length(l$qlev),
                          tdi = (length(l$qlev)+1):length(l$start)
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
    tminds <- mapply(function(f, s)f[s], tfinds, lapply(indexes, "[[", "tmi"), SIMPLIFY = FALSE)
    tdinds <- mapply(function(f, s)f[s], tfinds, lapply(indexes, "[[", "tdi"), SIMPLIFY = FALSE)
            
    rfd <- mapply(function(l, x, xd, tf, tm, td) c(l, list(index=list(x = x, xd = xd, tf = tf, tm = tm, td = td))),
                      rfd, xinds, xdinds, tfinds, tminds, tdinds, SIMPLIFY = FALSE)

    discount_term <- function(l,p) {
        if(is.null(l$discount_factor)) return(X[, l$index[["x"]]])                          
        dcf <- l$discount_factor(p[l$index$td])
        mqC(X[, l$index$x], l$qlev, l$lag_structure, dcf)
    }
    
    quantile_discounted<- function(X, par) {
        if(discountX) {
            xcols <- vector("list",length(rfd))
            for(i in 1:length(rfd)) xcols[[i]] <- discount_term(rfd[[i]],par)
            do.call("cbind", xcols)                    
        } else X
    }
    
    tm <- unlist(lapply(lapply(rfd, "[[", "index"),"[[","tm"))

    starto <- unlist(lapply(rfd, "[[", "start"))
    names(starto) <- unlist(lapply(rfd, "[[", "term_name"))
    Xd <- quantile_discounted(X,starto)

    md <- na.omit(cbind(y,Xd))
    if(is.null(attributes(md)[["na.action"]])) {
        ind <- 1:length(y)
    } else {
          ind <- (1:length(y))[-attr(md,"na.action")]
      }    

    if(Ofunction != "lm") {
        starto[tm] <- coef(lsfit(Xd[ind,], y[ind], intercept = FALSE))
    } 

    rhs <- function(p, ...) {
        Xd <- quantile_discounted(X, p)[ind, ]
        coefs <- p[tm]
        Xd %*% coefs
    }
    
    fn0 <- function(p, ...) {
        r <- y[ind] - rhs(p)
        sum(r^2, na.rm = TRUE)
    }

    gr <- function(p)grad(fn0,p)
    hess <- function(x)numDeriv::hessian(fn0,x)
    
    control <- c(list(Ofunction=Ofunction),args)
    ##Override default method of optim. Use BFGS instead of Nelder-Mead
    if(!("method"%in% names(control)) & Ofunction=="optim") {        
        control$method <- "BFGS"
    }
    
    out <- list(coefficients = starto,
                model = mf,
                dX = function(p) quantile_discounted(X, p)[ind, ],
                term_info = rfd,
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
                estimation_sample = ind,
                Zenv = Zenv)                
}


mqr.fit <- function(x) {
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
        y <- model.response(x$model, "numeric")[x$estimation_sample]
        fit <- lsfit(x$dX(coef(x)), y, intercept = FALSE)
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
            stop("The optimisation algorithm of Moving quantile regression failed with the following message:\n", opt,"\nPlease try other starting values or a different optimisation function")
        }
        par <- coef(opt)
        names(par) <- names(coef(x))
        x$convergence <- opt$convInfo$stopCode
    }
    
    x$opt <- opt
    x$coefficients <- par    
    
    x$fitted.values <- as.vector(x$rhs(coef(x)))

    y <- model.response(x$model, "numeric")[x$estimation_sample]
    x$residuals <- as.vector(y - x$fitted.values)
    x
}
    
