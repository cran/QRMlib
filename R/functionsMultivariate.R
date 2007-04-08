# QRMlib: version 1.4
# this file is a component of QRMlib 

# Copyright (C) 2005-06 Alexander McNeil 
# R-language additions Copyright (C) 2006-2007 Scott Ulman

# This program is free software; you can redistribute it and/or 
# modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation; either version 2 
# of the License, or (at your option) any later version. 

# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

# You should have received a copy of the GNU General Public License 
# along with this program; if not, write to the Free Software 
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA

# Contact: Alexander J. McNeil:  mcneil@math.ethz.ch */
# R-language contact: Scott Ulman : scottulman@hotmail.com
# Note in the R-translations that TRUE has been substituted throughout the 
# code for T (and FALSE for F) when setting default parameter values as 
# specified in section 2.4 of an "Introduction to R" (R-intro.pdf).
# Otherwise "R CMD check -QRMlib" returned the following type of error when
# trying to run an example:
# > BiDensPlot(func=dmnorm,mu=c(0,0),Sigma=equicorr(2,-0.7))
#   Error in func(data, ...) : F used instead of FALSE
#   Execution halted
############################################################
#R-language already has a function in the stats library called
# cov2cor() which performs this functionality.  Hence this function 
# is redundant in the R-language.
CovToCor <- function(mat)
{
  sds <- diag(1./sqrt(diag(mat)))
  out <- sds %*% mat %*% sds
  dimnames(out) <- dimnames(mat)
  out
}
##################################################

eigenmeth <- function(mat, delta = 0.001)
{
  decomp <- eigen(mat)
  Lambda <- decomp$values
  Lambda[Lambda < 0] <- delta
  Gamma <- decomp$vectors
  newmat <- Gamma %*% diag(Lambda) %*% t(Gamma)
  D <- 1/sqrt(diag(newmat))
  diag(D) %*% newmat %*% diag(D)
}

########################################################################

Spearman <- function(data){

  cor(apply(data, 2, rank))
}

#######################################################################

Kendall <- function(data, noforce = TRUE)
{
  n <- dim(data)[1]
  if((n > 5000) & (noforce))
    stop("Too many data - would be slow. Set noforce to FALSE for override"
         )
  d <- dim(data)[2]
  out <- matrix(1, nrow = d, ncol = d)
  for(i in (1:(d - 1))) {
    for(j in ((i + 1):d)) {
      out[i, j] <- cor.test(data[, i], data[, j], method = 
                            "kendall")$estimate
      out[j, i] <- out[i, j]
    }
  }
  nms <- dimnames(data)[[2]]
  dimnames(out) <- list(nms, nms)
  out
}

################################################################

equicorr <- function(d, rho)
{
  if(rho < ( - (d - 1.)^(-1.)))
    stop(paste("rho must be at least",  - (d - 1.)^(-1.)))
  J <- matrix(rho, nrow = d, ncol = d)
  D <- diag(rep(1. - rho, d))
  J + D
}

########################################################################

rmnorm <- function(n, Sigma = equicorr(d, rho), mu = rep(0, d), d=2, rho=0.7)
{
  d <- dim(Sigma)[1]
  A <- t(chol(Sigma))
  X <- matrix(rnorm(n * d), nrow = n, ncol = d)
  mu.matrix <- matrix(mu, nrow = n, ncol = d, byrow = TRUE)
  return(t(A %*% t(X)) + mu.matrix)
}

#################################################

rmt <- function(n, df = 4, Sigma = equicorr(d, rho), mu = rep(0, d), d=2, rho=0.7)
{
  d <- dim(Sigma)[1.]
  chi <- 2. * rgamma(n, shape = df/2.)
  m1 <- rmnorm(n, Sigma = Sigma)
  m2 <- matrix(rep(sqrt(df)/sqrt(chi), d), ncol = d)
  mu.matrix <- matrix(mu, nrow = n, ncol = d, byrow = TRUE)
  return(m1 * m2 + mu.matrix)
}
  
############################################################################################

dmnorm <- function(x,mu,Sigma,logvalue=FALSE)
{
  d <- dim(x)[2]
  Q <- mahalanobis(x,mu,Sigma)
  log.const.bottom <- d*log(2*pi)/2 +0.5*log(det(Sigma))
  log.top <- -Q/2
  out <- log.top-log.const.bottom
  if(!(logvalue))
    out <- exp(out)
  out
}

##########################################################################################

dmt <- function(x, nu, mu, Sigma, logvalue = FALSE)
{
  d <- dim(x)[2]
  Q <- mahalanobis(x, mu, Sigma)
  log.const.top <- log(gamma((nu + d)/2))
  log.const.bottom <- log(gamma(nu/2)) + (d * log(pi * nu))/2 + 0.5 *
    log(det(Sigma))
  log.top <- ( - (nu + d) * log(1 + Q/nu))/2
  out <- log.const.top + log.top - log.const.bottom
  if(!(logvalue))
    out <- exp(out)
  out
}

#############################################################

fit.norm <- function(data)
{
  #SU: 06/15/2006: the following works in R only if you have the fSeries package loaded.
  # It loads the RMetrics fCalendar timeSeries class.  Since the other fixes for timeSeries
  # assume we are using the RMetrics package, this is OK.
  if(class(data) == "timeSeries")
    data <- seriesData(data)
  if (is.matrix(data)==FALSE)
    data <- as.matrix(data)
  mu <- apply(data, 2., mean)

  #SU 06/27/2006: The R-function var() does not have an unbiased=F alternative. 
  #In S-Plus, unbiased=F implies the sum will be divided by N (the series length) rather than 
  # by N-1 which is always used by R. The value will be biased and hence smaller.
  #Hence we need to multiply the R-var() by (N-1)/N to get the unbiased=F (biased)S-Plus-var()to match
  # McNeil's code.  Hence we comment out the following S-Plus code from McNeil.
  #Sigma <- var(data, unbiased = F)
  #and replace it with this code:
  N = length(data)
  if(N==1) stop ("only one observation in data sent to fit.norm")
  Sigma <- (N-1)*var(data)/N  # creates the 'biased' (smaller) version which is divided by N
  
  d <- dim(data)[2]
  cor <- NA
  if(d > 1.) cor <- CovToCor(Sigma)
  maxloglik <- sum(dmnorm(data,mu,Sigma,logvalue=TRUE))
  out <- list(mu = mu, Sigma = Sigma, cor = cor, ll.max
              = maxloglik)
  out
}


##########################################

jointnormalTest <- function(data,dist="chisquare")
{
  #SU: 06/15/2006: the following works in R only if you have the fSeries package loaded.
  # It loads the RMetrics fCalendar timeSeries class.  Since the other fixes for timeSeries
  # assume we are using the RMetrics package, this is OK.
  if(class(data) == "timeSeries")
    data <- seriesData(data)
  d <- dim(data)[2]
  n <- dim(data)[1]
  mu <- apply(data, 2, mean)
  Sigma <- var(data)
  if (dist=="chisquare"){
    D <- mahalanobis(data, mu, Sigma)
    plot(qchisq(ppoints(D), d), sort(D), xlab = paste("Chi-sq(", d, ") quantiles", sep = ""), ylab = paste("Ordered D data", sep = ""))
    #SU: 06/15/2006: this is a single-sample Kolmogorov-Smirnoff goodness-of-fit test.  The S-Plus
    #function ks.gof() must be replaced by the R-function ks.test() contained in the stats package which will
    #normally be loaded:
    # pval <- ks.gof(D, distribution = "chisquare", df = d)$p.value  #S-Plus version
    #Note that in stats, the distribution is 'pchisq' with a preceding p indicating the probability 
    #distribution; you must use the preceding p. Also help("pchisq") shows the parameter name is df which you
    #must also set.
    #The '$p.value' will print only the p-value and ignore the other list elements which ks.test() returns
    pval <- ks.test(D, "pchisq", df = d)$p.value
  } 
  else if (dist=="beta"){
    D <- mahalanobis(data, mu, Sigma)*n*(n-1)^(-2)
    plot(qbeta(ppoints(D), d/2,(n-d-1)/2), sort(D), xlab = paste("Beta(", d/2,",",(n-d-1)/2,") quantiles", sep = ""),
         ylab = paste("Ordered D^2 data (scaled)", sep =""))
    #SU: 06/15/2006: this is a single-sample Kolmogorov-Smirnoff goodness-of-fit test.  The S-Plus
    #function ks.gof() must be replaced by the R-function ks.test() contained in the stats package which will
    #normally be loaded:
    #pval <- ks.gof(D, distribution = "beta", shape1=d/2,shape2=(n-d-1)/2)$p.value
    #Note that in stats, the distribution is 'pbeta' with a preceding 'p' for probability or 
    # distribution function; you must match the spelling including the preceding p. Also according to
    #help("pbeta"), the parameters are named shape1 and shape2.
    #The '$p.value' will print only the p-value and ignore the other list elements which ks.test() returns
    pval <- ks.test(D, "pbeta", shape1=d/2, shape2=(n-d-1)/2)$p.value

  }
  else stop("Unknown reference distribution")
  abline(0, 1)
  return(paste("KS pvalue",round(pval,3)))
}

######################################################

MardiaTest <- function(data)
  {
    #SU: 06/15/2006: the following works in R only if you have the fSeries package loaded.
    # It loads the RMetrics fCalendar timeSeries class.  Since the other fixes for timeSeries
    # assume we are using the RMetrics package, this is OK.
    if(class(data) == "timeSeries")
		data <- seriesData(data)
    d <- dim(data)[2]
    n <- dim(data)[1]
    Xbar <- apply(data,2,mean)
    Xbar.matrix <- matrix(Xbar,nrow=n,ncol=d,byrow=TRUE)
    standardised <- data -Xbar.matrix
    S <- var(data)
    A <- t(chol(S))
    Ainv <- solve(A)
    Zdata <- Ainv %*% t(standardised)
    Zdata <- t(Zdata)
    D2.check <- mahalanobis(data,Xbar,S)
    Dij <- Zdata %*% t(Zdata)
    D2 <- diag(Dij)
    K3 <- mean(as.vector(Dij)^3)
    K4 <- mean(D2^2)
    statK3 <- n*K3/6
    df <- d*(d+1)*(d+2)/6
    K3.pval <- 1-pchisq(statK3,df)
    mn <- d*(d+2)
    vr <- 8*d*(d+2)/n
    K4.stat <- (K4-mn)/sqrt(vr)
    K4.pval <- 1-pnorm(abs(K4.stat))
    c(K3,K3.pval,K4,K4.pval)
  }

#########################################################

#SU: 06/15/2006: R-language reports the following error from this S-Plus function:
#Error in matrix(z, nrow = length(x), nrow = length(y), byrow = TRUE) : 
#        formal argument "nrow" matched by multiple actual arguments
BiDensPlot <- function(func,xpts =  c(-2, 2), ypts = c(-2, 2), npts=50, type="persp", ...)
{
  #SU: 06/15/2006: the function is seq(from, to, length.out= ) rather than 
  # the S-Plus version of seq(from, to , length =)
  # x <- seq(from = xpts[1], to = xpts[2], length = npts) #S-Plus version
  #y <- seq(from = ypts[1], to = ypts[2], length = npts)  #S-Plus version
  x <- seq(from = xpts[1], to = xpts[2], length.out = npts)
  y <- seq(from = ypts[1], to = ypts[2], length.out = npts)

  xval <- rep(x, length(y))
  yval <- rep(y, rep(length(x), length(y)))
  data <- cbind(xval, yval)
  z <- eval(func(data, ...))
  #SU: 06/15/2006: major error here: 2nd nrow= must be replaced by ncol=
  #z <- matrix(z, nrow = length(x), nrow = length(y), byrow = T)
  z <- matrix(z, nrow = length(x), ncol = length(y), byrow = TRUE)
  switch(type,persp = persp(x, y, z), contour =contour(x, y, z))
  return(invisible())
}

#################################################################

