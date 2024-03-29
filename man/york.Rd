% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/york.R
\name{york}
\alias{york}
\title{York Regression}
\usage{
york(data, intercept = NA, method = "LS", tol = 1e-15, maxit = 50,
  verbose = FALSE, gridstartvals = TRUE)
}
\arguments{
\item{data}{Data matrix. Assumed to be a matrix with columns containing (in order):
X, standard error of X, Y, standard error of Y, and residual correlation.
The correlation column may be omitted, in which case the correlation is
assumed to be 0.}

\item{intercept}{Regression intercept. If 'NA' intercept is freely estimated, otherwise fixed to the given value.}

\item{method}{Method for computing standard errors. Can be one
of "LS" (original least squares SE from York) or "ML" (Titterington).}

\item{tol}{Tolerance for evaluating convergence, based on relative change in the slope.}

\item{maxit}{Maximum number of iterations}

\item{verbose}{Boolean, whether to print parameter estimates and loss function at each iteration step}

\item{gridstartvals}{Boolean, whether to pick initial slope estimate from a grid around a set of possible WLS regression fits.
Default (false) means start value are chosen from WLS with weights proportional to variance of X-Y.}
}
\value{
List with elements:
\describe{

  \item{\code{converge}}{list indicating whether the regression fit converged and the tolerance criteria used}

  \item{\code{a}}{fitted intercept and standard error}

  \item{\code{b}}{fitted slope and standard error}

  \item{\code{cov.ab}}{covariance of the slope and intercept estimates}

  \item{\code{chi2}}{goodness of fit statistic}

  \item{\code{df}}{degrees of freedom}

  \item{\code{mswd}}{weighted mean squared residual}

  \item{\code{p.value}}{test of goodness of fit}

}
}
\description{
Regression with correlated sampling error in X and Y, following York (1968) and extensions
}
\details{
...
}
\examples{

n <- 50
a0 <- 1
b <- .75
xtrue <- rnorm(n,4,1.5)
ytrue <- a0 + b*xtrue
xsd <- sqrt(runif(n,.05,1))
ysd <- sqrt(runif(n,.05,0.5))
re <- rep(0.3,n)
obs <- t(sapply(1:n,
    function(i){
        MASS::mvrnorm(1,
            mu=c(xtrue[i],ytrue[i]),
            Sigma = matrix(c(xsd[i]^2, re[i]*xsd[i]*ysd[i], re[i]*xsd[i]*ysd[i], ysd[i]^2),2,2)
        )
    }))
dat <-data.frame(x=obs[,1],sx=xsd,y=obs[,2],sy=ysd,re=re)

# IsoplotR::york(dat)

f0 <- lm(y~x, data=dat)
f1 <- lm(y~x, data=dat, weight=1/dat$sy^2)
f2 <- york(dat)
f3 <- york(dat, intercept=0)

\dontrun{
plot(dat$x, dat$y, xlim=c(0,7.5), ylim=c(0,7.5), pch=20)
arrows(x0=dat$x-dat$sx,x1=dat$x+dat$sx,y0=dat$y,y1=dat$y,code=3,angle=90,length=0.02)
arrows(x0=dat$x,x1=dat$x,y0=dat$y-dat$sy,y1=dat$y+dat$sy,code=3,angle=90,length=0.02)
abline(a0, b, col="black", lwd=2)
abline(f0$coefficients[1], f0$coefficients[2], col="blue", lty=2)
abline(f1$coefficients[1], f1$coefficients[2], col="red", lty=2)
abline(f2$a[1], f2$b[1], col="green")
abline(f3$a[1], f3$b[1], col="purple")
dev.off()
}


}
\references{
York, D., 1968. Least squares fitting of a straight line with
correlated errors. Earth and Planetary Science Letters, 5, pp.320-324.

Titterington, D.M. and Halliday, A.N., 1979. On the fitting of
parallel isochrons and the method of maximum likelihood. Chemical
Geology, 26(3), pp.183-195.

Minster, J.F., et al., 1979. 87Rb-87Sr chronology of enstatite
meteorites. Earth and Planetary Science Letters, 44(3), pp.420-440.

York, D., et al., 2004. Unified equations for the slope,
intercept, and standard errors of the best straight line. American
Journal of Physics 72.3, pp.367-375.
}
