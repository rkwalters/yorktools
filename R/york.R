#' York Regression
#'
#' @description
#' Regression with correlated sampling error in X and Y, following York (1968) and extensions
#'
#' @details
#' ...
#'
#' @param data
#' Data matrix. Assumed to be a matrix with columns containing (in order):
#' X, standard error of X, Y, standard error of Y, and residual correlation.
#' The correlation column may be omitted, in which case the correlation is
#' assumed to be 0.
#' @param intercept
#' Regression intercept. If 'NA' intercept is freely estimated, otherwise fixed to the given value.
#' @param method
#' Method for computing standard errors. Can be one
#' of "LS" (original least squares SE from York) or "ML" (Titterington).
#' @param tol
#' Tolerance for evaluating convergence, based on relative change in the slope.
#' @param maxit
#' Maximum number of iterations
#' @param verbose
#' Boolean, whether to print parameter estimates and loss function at each iteration step
#' @param gridstartvals
#' Boolean, whether to pick initial slope estimate from a grid around a set of possible WLS regression fits.
#' Default (false) means start value are chosen from WLS with weights proportional to variance of X-Y.
#'
#' @return
#' List with elements:
#' \describe{
#'
#'   \item{\code{converge}}{list indicating whether the regression fit converged and the tolerance criteria used}
#'
#'   \item{\code{a}}{fitted intercept and standard error}
#'
#'   \item{\code{b}}{fitted slope and standard error}
#'
#'   \item{\code{cov.ab}}{covariance of the slope and intercept estimates}
#'
#'   \item{\code{chi2}}{goodness of fit statistic}
#'
#'   \item{\code{df}}{degrees of freedom}
#'
#'   \item{\code{mswd}}{weighted mean squared residual}
#'
#'   \item{\code{p.value}}{test of goodness of fit}
#'
#' }
#'
#' @export
#'
#' @references
#' York, D., 1968. Least squares fitting of a straight line with
#' correlated errors. Earth and Planetary Science Letters, 5, pp.320-324.
#'
#' Titterington, D.M. and Halliday, A.N., 1979. On the fitting of
#' parallel isochrons and the method of maximum likelihood. Chemical
#' Geology, 26(3), pp.183-195.
#'
#' Minster, J.F., et al., 1979. 87Rb-87Sr chronology of enstatite
#' meteorites. Earth and Planetary Science Letters, 44(3), pp.420-440.
#'
#' York, D., et al., 2004. Unified equations for the slope,
#' intercept, and standard errors of the best straight line. American
#' Journal of Physics 72.3, pp.367-375.
#'
#'
#' @examples
#'
#' n <- 50
#' a0 <- 1
#' b <- .75
#' xtrue <- rnorm(n,4,1.5)
#' ytrue <- a0 + b*xtrue
#' xsd <- sqrt(runif(n,.05,1))
#' ysd <- sqrt(runif(n,.05,0.5))
#' re <- rep(0.3,n)
#' obs <- t(sapply(1:n,
#'     function(i){
#'         MASS::mvrnorm(1,
#'             mu=c(xtrue[i],ytrue[i]),
#'             Sigma = matrix(c(xsd[i]^2, re[i]*xsd[i]*ysd[i], re[i]*xsd[i]*ysd[i], ysd[i]^2),2,2)
#'         )
#'     }))
#' dat <-data.frame(x=obs[,1],sx=xsd,y=obs[,2],sy=ysd,re=re)
#'
#' # IsoplotR::york(dat)
#'
#' f0 <- lm(y~x, data=dat)
#' f1 <- lm(y~x, data=dat, weight=1/dat$sy^2)
#' f2 <- york(dat)
#' f3 <- york(dat, intercept=0)
#'
#' \dontrun{
#' plot(dat$x, dat$y, xlim=c(0,7.5), ylim=c(0,7.5), pch=20)
#' arrows(x0=dat$x-dat$sx,x1=dat$x+dat$sx,y0=dat$y,y1=dat$y,code=3,angle=90,length=0.02)
#' arrows(x0=dat$x,x1=dat$x,y0=dat$y-dat$sy,y1=dat$y+dat$sy,code=3,angle=90,length=0.02)
#' abline(a0, b, col="black", lwd=2)
#' abline(f0$coefficients[1], f0$coefficients[2], col="blue", lty=2)
#' abline(f1$coefficients[1], f1$coefficients[2], col="red", lty=2)
#' abline(f2$a[1], f2$b[1], col="green")
#' abline(f3$a[1], f3$b[1], col="purple")
#' dev.off()
#' }
#'
#'


# adapted from isoplotR function york()
york <- function(data,intercept=NA,method="LS",tol=1e-12,maxit=5000,verbose=FALSE,gridstartvals=TRUE){

  # input check
  dat <- as.data.frame(data)

  if(!(dim(dat)[2] %in% c(4,5))){
    stop("data must have 4 or 5 columns.")
  }

  if(!(method=="LS" || method=="ML")){
    stop("method must be one of 'LS' (least squares) or 'ML' (maximum likelihood)")
  }

  if(method=="ML" && !is.na(intercept)){
    warning("ML standard errors for fixed intercept regression are unverified")
  }

  # get data from data frame
  if (ncol(dat)==4){
    dat <- cbind(dat,0)
  }
  colnames(dat) <- c('X','sX','Y','sY','rXY')

  # drop NAs
  na_rows <- is.na(rowSums(dat))
  if(any(na_rows)){
    n_drop <- sum(na_rows)
    warning("Omitting ",n_drop," entries with missing data.")
    dat <- dat[!na_rows,]
  }

  # extract values
  X <- dat[,'X']
  sX <- dat[,'sX']
  wX <- 1/sX^2
  Y <- dat[,'Y']
  sY <- dat[,'sY']
  wY <- 1/sY^2
  rXY <- dat[,'rXY']
  covXY <- rXY*sX*sY

  # initial guesses
  if(!gridstartvals){
    if(is.na(intercept)){
      mod <- stats::lm(Y ~ X, weights=1/((sX^2) + (sY^2) - 2*covXY))$coefficients
      if(any(is.na(mod))){
        stop('Cannot fit a straight line through these data')
      }
      b <- mod[2]
      a <- mod[1]
    }else{
      modb <- stats::lm(Y ~ 0 + X, offset = rep(intercept,nrow(dat)), weights=1/((sX^2) + (sY^2) - 2*covXY))$coefficients
      if(any(is.na(modb))){
        stop('Cannot fit a straight line through these data')
      }
      b <- modb[1]
      a <- intercept
    }

    Winit <- 1/(sY^2 + (b^2)*sX^2 - 2*b*covXY)
    loss <- sum(Winit*(Y-X*b-a)^2)

  }else{
    # test multiple lm fits + surrounding grid to help avoid local minima
    if(is.na(intercept)){
      mod1b <- stats::lm(Y ~ X)$coefficients[2]
      if(any(is.na(mod1b))){
        stop('Cannot fit a straight line through these data')
      }
      mod2b <- stats::lm(Y ~ X, weights=wY)$coefficients[2]
      mod3b <- stats::lm(Y ~ X, weights=wX)$coefficients[2]
      mod4b <- stats::lm(Y ~ X, weights=1/((sX^2) + (sY^2)))$coefficients[2]
      mod5b <- stats::lm(Y ~ X, weights=1/((sX^2) + (sY^2) - 2*covXY))$coefficients[2]

    }else{
      mod1b <- stats::lm(Y ~ 0 + X, offset = rep(intercept,nrow(dat)))$coefficients[1]
      if(any(is.na(mod1b))){
        stop('Cannot fit a straight line through these data')
      }
      mod2b <- stats::lm(Y ~ 0 + X, offset = rep(intercept,nrow(dat)), weights=wY)$coefficients[1]
      mod3b <- stats::lm(Y ~ 0 + X, offset = rep(intercept,nrow(dat)), weights=wX)$coefficients[1]
      mod4b <- stats::lm(Y ~ 0 + X, offset = rep(intercept,nrow(dat)), weights=1/((sX^2) + (sY^2)))$coefficients[1]
      mod5b <- stats::lm(Y ~ 0 + X, offset = rep(intercept,nrow(dat)), weights=1/((sX^2) + (sY^2) - 2*covXY))$coefficients[1]
    }

    maxb <- max(abs(c(mod1b, mod2b, mod3b, mod4b, mod5b)))
    qq <- seq(-5*maxb, 5*maxb, length=10)
    qq <- c(qq,0,mod1b, mod2b, mod3b, mod4b, mod5b)

    init_loss <- rep(NA, length(qq))
    ainit <- rep(NA, length(qq))

    for(i in 1:length(qq)){
      Winit <- wX*wY/(wX+(qq[i]^2)*wY-2*qq[i]*dat[,'rXY']*sqrt(wX*wY))
      if(is.na(intercept)){
        Xbar <- sum(Winit*X)/sum(Winit)
        Ybar <- sum(Winit*Y)/sum(Winit)
        ainit[i] <- Ybar-qq[i]*Xbar
      }else{
        ainit[i] <- intercept
      }
      init_loss[i] <- sum(Winit*(Y-X*qq[i]-ainit[i])^2)
    }

    best_start <- which.min(init_loss)
    b <- qq[best_start]
    a <- ainit[best_start]
    loss <- init_loss[best_start]
  }

  # iterative fit
  converge=FALSE
  if(verbose){
    print(paste0("starting values: a: ",a,"; b: ",b))
  }
  for (i in 1:maxit){
    bold <- b
    lossold <- loss
    A <- sqrt(wX*wY)
    W <- 1/(sY^2 + (b^2)*sX^2 - 2*b*covXY)
    # possible safeguard, but shouldn't be necessary (and probably right to error out if it happens)
    # Xbar <- sum(W*X,na.rm=TRUE)/sum(W,na.rm=TRUE)
    # Ybar <- sum(W*Y,na.rm=TRUE)/sum(W,na.rm=TRUE)
    Xbar <- sum(W*X)/sum(W)
    Ybar <- sum(W*Y)/sum(W)

    if(is.na(intercept)){
      U <- X - Xbar
      V <- Y - Ybar
      Beta <- W*(U*sY^2 + b*V*sX^2 - (b*U+V)*covXY)

      b <- sum(W*Beta*V)/sum(W*Beta*U)
      a <- Ybar-b*Xbar
    }else{
      a <- intercept # redundant, but for clarity
      V <- Y-a

      b <- sum(W^2*V*(X*sY^2 + b*V*sX^2 - V*covXY))/sum(W^2*X*(X*sY^2 + b*V*sX^2 - b*X*covXY))
    }

    loss <- sum(W*(Y-X*b-a)^2)
    # if(is.na(intercept)){
    #   xpred <- Xbar + Beta
    #   ypred <- Ybar + b*Beta
    # }else{
    #   xpred <- W*(X*(sY^2-b*covXY)+(Y-intercept)*(b*sX^2-covXY))
    #   ypred <- b*xpred+a
    # }
    # loss2 <- sum(((X-xpred)^2/sX^2 + (Y-ypred)^2/sY^2 - 2*dat[,'rXY']*(X-xpred)*(Y-ypred)/(sX*sY))/(1-dat[,'rXY'])^2)

    if(verbose){
      print(paste0("iteration: ",i,"; a: ",a,"; b: ",b,"; loss: ",loss)) #, "; loss2: ",loss2))
    }
    if ((bold/b-1)^2 < tol){ # convergence reached
      converge <- TRUE
      break
    }
  }

  # standard error
  if(converge){
    if(is.na(intercept)){

      df <- nrow(dat)-2

      if(method=="LS"){
        BBar <- sum(W*Beta)/sum(W)
        num <- sum(W^2*(U^2*sY^2 + V^2*sX^2 - 2*U*V*covXY))
        d1 <- (1/b)*sum(W*U*V)
        d2 <- 4*sum(W*(Beta-U)*(Beta-BBar))
        d3 <- (1/b)*sum(W^2*(b*U-V)^2*covXY)
        D <- d1+d2-d3
        sb <- sqrt(num/(D^2))

        sa <- sqrt((1/sum(W)) + (Xbar+2*BBar)^2*sb^2 + (2*(Xbar-2*BBar)*BBar)/D)

        cov.ab <- -(Xbar+2*BBar)*sb^2 - BBar/D # Minster, plus York 2004 comment of error in Q

      }else if(method=="ML"){
        xpred <- Xbar + Beta
        ypred <- Ybar + b*Beta

        xpredbar <- sum(W*xpred)/sum(W)
        upred <- xpred-xpredbar

        sb <- sqrt(1/(sum(W*upred^2)))
        sa <- sqrt((1/sum(W))+xpredbar^2*sb^2)

        cov.ab <- -xpredbar*sb^2

      }
    }else{
      if(method=="LS"){
        Beta <- W*(X*sY^2 + b*V*sX^2 - (b*X+V)*covXY)
        num <- sum(W^2*(X^2*sY^2 + V^2*sX^2 - 2*X*V*covXY))
        d1 <- (1/b)*sum(W*X*V)
        d2 <- 4*sum(W*(Beta-X)*Beta)
        d3 <- (1/b)*sum(W*(b*X-V)*Beta)
        d4 <- (1/b)*sum(W^2*(b*X-V)^2*covXY)
        D <- d1+d2+d3-d4
        sb <- sqrt(num/(D^2))

      }else if(method=="ML"){
        xpred <- W*(X*sY^2 + b*V*sX^2 - (b*X+V)*covXY)
        sb <- 1/sum(W*xpred^2)
      }

      sa <- NA
      cov.ab <- NA
      df <- nrow(dat)-1

    }

    chi2 <- sum(W*(Y-X*b-a)^2)


  }else{
    sb <- NA
    sa <- NA
    cov.ab <- NA
    chi2 <- NA
    df <- NA
    warning("York with fixed intercept failed to converge")
  }

  # build output
  out <- list()
  out$a <- c(a,sa)
  out$b <- c(b,sb)
  out$cov.ab <- cov.ab
  names(out$a) <- c('a','s[a]')
  names(out$b) <- c('b','s[b]')

  out$df <- df
  out$chi2 <- chi2
  out$mswd <- chi2/df
  out$p.value <- stats::pchisq(chi2,out$df,lower=F)

  out$converge <- list(converged=converge, tol=tol)
  out$type <- 'york'
  class(out) <- 'york'

  return(out)
}


#' York Regression predicted values
#'
#' @description
#' Fitted values method for use with york regression output
#'
#' @details
#' ...
#'
#' @param data
#' Data matrix. Assumed to be a matrix with columns containing (in order):
#' X, standard error of X, Y, standard error of Y, and residual correlation.
#' The correlation column may be omitted, in which case the correlation is
#' assumed to be 0.
#'
#' @param model
#' York regression model, as output by `york()`. Alternatively, can specify model by providing `slope` and `intercept` values.
#'
#' @param intercept
#' Scalar, York regression model intercept.
#'
#' @param slope
#' Scalar, York regression model slope.
#'
#' @return
#' Data frame of predicted values for each observation with columns:
#'
#' @export
pred.york <- function(data, model=NULL, intercept=NA, slope=NA){

  if(is.null(model) && is.na(slope)){
    stop("must specify either model or slope")
  }

  if(!is.null(model) && !is.na(slope)){
    stop("must specify only one of 'model' or 'slope'")
  }

  if(!is.null(model) && (class(model)!="york") ){
    stop("model must be from york regression")
  }

  if(!is.null(model)){
    slope <- model$b[1]
    intercept <- model$a[1]
  }

  dat <- as.data.frame(data)

  if(!(dim(dat)[2] %in% c(4,5))){
    stop("data must have 4 or 5 columns.")
  }

  # get data from data frame
  if (ncol(dat)==4){
    dat <- cbind(dat,0)
  }
  colnames(dat) <- c('X','sX','Y','sY','rXY')



  W <- 1/(dat$sY^2 + (slope^2)*dat$sX^2 - 2*slope*dat$rXY*dat$sY*dat$sX)

  xpred <- W*(dat$X*(dat$sY^2-slope*dat$rXY*dat$sY*dat$sX)+(dat$Y-intercept)*(slope*dat$sX^2-dat$rXY*dat$sY*dat$sX))
  ypred <- intercept + slope*xpred

  out <- data.frame(xpred=xpred, ypred=ypred)
  row.names(out) <- row.names(data)
  return(out)
}

#' York Regression predicted values
#'
#' @description
#' Fitted values method for use with york regression output
#'
#' @details
#' Wrapper for more standard interface to `pred.york()`
#'
#' @param object
#' York regression model, as output by `york()`
#'
#' @param ...
#' Must include data matrix named `data`. Assumed to be a matrix with columns containing (in order):
#' X, standard error of X, Y, standard error of Y, and residual correlation.
#' The correlation column may be omitted, in which case the correlation is
#' assumed to be 0. All other arguments ignored.
#'
#' @return
#' Data frame of predicted values for each observation with columns:
#'
#' @export
fitted.york <- function(object, ...){

  args <- list(...)

  if(exists(args$data)){
    pred.york(args$data, model=object)
  }else{
    stop("York regression fitted values requires specifying separate `data`")
  }
}


#' York Regression coefficients
#'
#' @description
#' Coefficient methods to extract York regression model parameters
#'
#' @details
#' ...
#'
#' @param object
#' York regression model, as output by `york()`
#'
#' @param ...
#' Ignored
#'
#' @return
#' Estimates and standard errors from the fitted York regression model
#'
#'
#' @export
coef.york <- function(object, ...){
  if(class(object)!="york" ){
    stop("model must be from york regression")
  }

  out <- data.frame(Estimate=c(object$a[1], object$b[1]), `Std. Error`=c(object$a[2], object$b[2]))
  rownames(out) <- c("Intercept","x")

  return(out)
}


#' York Regression loss
#'
#' @description
#' Residuals and loss function values for each observation in a fitted York regression
#'
#' @details
#' ...
#'
#' @param data
#' Data matrix. Assumed to be a matrix with columns containing (in order):
#' X, standard error of X, Y, standard error of Y, and residual correlation.
#' The correlation column may be omitted, in which case the correlation is
#' assumed to be 0.
#'
#' @param model
#' York regression model, as output by `york()`. Alternatively, can specify model by providing `slope` and `intercept` values.
#'
#' @param intercept
#' Scalar, York regression model intercept.
#'
#' @param slope
#' Scalar, York regression model slope.
#'
#' @param func
#' Which loss function to report. Either "line" or "predicted_values". See Details.
#'
#' @return
#' Data frame of residuals for each observation with columns:
#'
#'
#' @export
loss_per_obs.york <- function(data, model=NULL, intercept=NA, slope=NA, func="line"){

  if(is.null(model) && is.na(slope)){
    stop("must specify either model or slope")
  }

  if(!is.null(model) && !is.na(slope)){
    stop("must specify only one of 'model' or 'slope'")
  }

  if(!is.null(model) && (class(model)!="york") ){
    stop("model must be from york regression")
  }

  if(!(func %in% c("line","predicted_values"))){
    stop("func must be either 'line' or 'predicted_values'")
  }

  if(!is.null(model)){
    slope <- model$b[1]
    intercept <- model$a[1]
  }

  dat <- as.data.frame(data)

  if(!(dim(dat)[2] %in% c(4,5))){
    stop("data must have 4 or 5 columns.")
  }

  # get data from data frame
  if (ncol(dat)==4){
    dat <- cbind(dat,0)
  }
  colnames(dat) <- c('X','sX','Y','sY','rXY')


  if(func=="line"){

    W <- 1/(dat$sY^2 + (slope^2)*dat$sX^2 - 2*slope*dat$rXY*dat$sY*dat$sX)

    resid <- dat$Y-intercept-slope*dat$X
    loss <- W*resid^2

    out <- data.frame(resid=resid, loss=loss)
    rownames(out) <- rownames(data)


  }else if(func=="predicted_values"){

    pred <- pred.york(dat, intercept=intercept, slope=slope)

    xresid <- dat$X-pred$xpred
    yresid <- dat$Y-pred$ypred

    loss <- (xresid^2/dat$sX^2 + yresid^2/dat$sY^2 - 2*dat$rXY*xresid*yresid/(dat$sX*dat$sY))/(1-dat$rXY)^2

    out <- data.frame(xresid=xresid, yresid=yresid, loss=loss)
    rownames(out) <- rownames(data)

  }

  return(out)

}



#' York Regression residuals
#'
#' @description
#' Residuals method for use with york regression output
#'
#' @details
#' ...
#'
#' @param object
#' York regression model, as output by `york()`
#'
#' @param ...
#' Must include data matrix named `data`. Assumed to be a matrix with columns containing (in order):
#' X, standard error of X, Y, standard error of Y, and residual correlation.
#' The correlation column may be omitted, in which case the correlation is
#' assumed to be 0. All other arguments ignored.
#'
#' @return
#' Data frame of residuals for each observation with columns:
#' \describe{
#'
#'   \item{\code{xresid}}{difference between X and predicted X value}
#'
#'   \item{\code{yresid}}{difference between Y and predicted Y value}
#'
#'   \item{\code{loss}}{weighted misfit}
#'
#' }
#'
#' @export

residuals.york <- function(object, ...){

  args <- list(...)

  if(exists(args$data)){
    loss_per_obs.york(args$data, model=object, func="predicted_values")
  }else{
    stop("York regression residuals requires specifying separate `data`")
  }
}
