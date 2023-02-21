#' York Regression
#'
#' @param dat
#' Data matrix. Assumed to be a matrix with columns containing (in order):
#' X, standard error of X, Y, standard error of Y, and residual correlation.
#' The correlation column may be omitted, in which case the correlation is
#' assumed to be 0.
#' @param intercept
#' Regression intercept. If 'NA' intercept is freely estimated, otherwise fixed to the given value.
#' @param method
#' Method for computing standard errors. Can be one of
#' "LS" (original least squares SE from York) or "ML" (Titterington).
#' @param tol
#' Tolerance for evaluating convergence, based on relative change in the slope.
#' @param maxit
#' Maximum number of iterations
#'
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
#' York, D., et al., 2004. Unified equations for the slope,
#' intercept, and standard errors of the best straight line. American
#' Journal of Physics 72.3, pp.367-375.
#'
#' Mahon, K.I., 1996. The new “York” regression: application of an
#' improved statistical method to geochemistry. International Geology
#' Review, 38(4), pp. 293-303.
#'
#'
#' @examples
#'
#' n <- 50
#' xt <- rnorm(n)
#' xs <- runif(n,.005,.05)
#'
#' yt <- sqrt(.9)*xt + sqrt(.1)*rnorm(n)
#' ys <- runif(n,.005,.05)
#'
#' rr <- .5
#'
#' rand_biv <- function(i){
#'   ss <- matrix(c(xs[i]^2, rr*xs[i]*ys[i], rr*xs[i]*ys[i], ys[i]), 2, 2)
#'   ee <- MASS::mvrnorm(1, mu=c(xt[i],yt[i]), Sigma=ss)
#' }
#'
#' eo <- t(sapply(1:n, rand_biv))
#'
#' xo <- xt + eo[,1]
#' yo <- yt + eo[,2]
#'
#' yo <- sign(xo)*yo
#' xo <- sign(xo)*xo
#'
#' dat <- data.frame(xo,xs,yo,ys,rep(rr,n))
#'
#' IsoplotR::york(dat)
#'
#' york(dat)
#'
#' york(dat, intercept=0)
#'
#' york(dat, method="ML")
#'


# adapted from isoplotR function york()
york <- function(dat,intercept=0,method="LS",tol=1e-15,maxit=50){

  # get data from data frame
  if (ncol(dat)==4) dat <- cbind(dat,0)
  colnames(dat) <- c('X','sX','Y','sY','rXY')
  dat <- dat[!is.na(rowSums(dat)),]
  X <- dat[,'X']
  sX <- dat[,'sX']
  wX <- 1/sX^2
  Y <- dat[,'Y']
  sY <- dat[,'sY']
  wY <- 1/sY^2
  covXY <- dat[,'rXY']*sX*sY

  # initial guesses
  ab <- stats::lm(Y ~ 0 + X, offset = rep(intercept,nrow(dat)))$coefficients
  b <- ab[1]
  if (any(is.na(ab)))
    stop('Cannot fit a straight line through these data')

  # iterative fit
  converge=FALSE
  for (i in 1:maxit){
    bold <- b
    A <- sqrt(wX*wY)
    W <- wX*wY/(wX+(b^2)*wY-2*b*dat[,'rXY']*A)
    Xbar <- sum(W*X,na.rm=TRUE)/sum(W,na.rm=TRUE)
    Ybar <- sum(W*Y,na.rm=TRUE)/sum(W,na.rm=TRUE)
    b <- (Ybar - intercept)/Xbar
    if ((bold/b-1)^2 < tol){ # convergence reached
      converge <- TRUE
      break
    }
  }

  # standard error
  if(converge){
    dx <- b*W / sum(W)
    dy <- -1*W / sum(W)
    dw <- 2*b*sX^2 - 2*covXY
    sum_dyw2 <- sum(dw * Y * (W^2) )
    sum_bdxw2 <- sum(b * dw * X * (W^2) )
    sum_bxw <- sum(b * X * W )
    sum_yw <- sum(Y * W)
    sum_w <- sum(W)
    sum_dw2 <- sum(dw * (W^2) )
    dtheta <- Xbar + ((sum_w*(sum_dyw2-sum_bdxw2) + sum_dw2*(sum_bxw-sum_yw))/(sum_w^2))

    sb <- sqrt( sum((dx^2)*(sX^2) + (dy^2)*(sY^2) + 2*dx*dy*covXY ) / dtheta^2)

    chi2 <- sum(W*(Y-X*b)^2)
    df <- nrow(dat)-1

  }else{
    sb <- NA
    stop("York with fixed intercept failed to converge")
  }

  # build output
  out <- list()
  out$converge <- list(converged=converge, tol=tol)
  out$a <- c(intercept,NA)
  out$b <- c(b,sb)
  out$cov.ab <- NA
  names(out$a) <- c('a','s[a]')
  names(out$b) <- c('b','s[b]')
  out$type <- 'york'

  out$chi2 <- chi2
  out$df <- nrow(dat)-1
  out$mswd <- chi2/df
  out$p.value <- as.numeric(1-stats::pchisq(chi2,out$df))

  return(out)
}
