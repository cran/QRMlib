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

# Contact: Alexander J. McNeil:  mcneil@math.ethz.ch 
# R-language contact: Scott Ulman : scottulman@hotmail.com 
# Note in the R-translations that TRUE has been substituted throughout the 
# code for T (and FALSE for F) when setting default parameter values as 
# specified in section 2.4 of an "Introduction to R" (R-intro.pdf).
# Otherwise "R CMD check -QRMlib" returned the following type of error when
# trying to run an example:
# > BiDensPlot(func=dmnorm,mu=c(0,0),Sigma=equicorr(2,-0.7))
#   Error in func(data, ...) : F used instead of FALSE
#   Execution halted
#############################################################################

rcopula.gauss <- function(n,  Sigma = equicorr(d, rho), d=2, rho=0.7)
{
  d <- dim(Sigma)[1]
  if(sum(diag(Sigma)) != d)
    stop("Sigma should be correlation matrix")
  mnorm <- rmnorm(n, Sigma = Sigma)
  matrix(pnorm(mnorm), ncol = d)
}

#################################################################################

rcopula.t <- function(n, df, Sigma = equicorr(d, rho), d=2, rho=0.7)
{
  d <- dim(Sigma)[1]
  if(sum(diag(Sigma)) != d)
    stop("Sigma should be correlation matrix")
  tmp <- rmt(n, df, Sigma)
  matrix(pt(tmp, df), ncol = d)
}

###################################################################################

rAC <- function(name, n,d,theta)
  #generic function for Archimedean copula simulation
  {
    illegalpar <- switch(name,
                         clayton=(theta<0),
                         gumbel=(theta<1),
                         frank=(theta<0),
                         BB9=((theta[1]<1) | (theta[2] <0)))
    if (illegalpar) stop("Illegal parameter value")
    independence <- switch(name,
                           clayton=(theta==0),
                           gumbel=(theta==1),
                           frank=(theta==0),
                           BB9=(theta[1]==1))
    U <- runif(n*d)
    U <- matrix(U,nrow=n,ncol=d)
    if (independence) return(U)
    Y <- switch(name,
                clayton=rgamma(n,1/theta),
                gumbel=rstable(n,1/theta)*(cos(pi/(2*theta)))^theta,
                frank=rFrankMix(n,theta),
                BB9=rBB9Mix(n,theta))
    Y <- matrix(Y, nrow=n,ncol=d)
    phi.inverse <- switch(name,
                          clayton=function(t,theta){(1+t)^(-1/theta)},
                          gumbel=function(t,theta){exp(-t^(1/theta))},
                          frank=function(t,theta){(-1/theta)*log(1-(1-exp(-theta))*exp(-t))},
                          BB9=function(t,theta){exp(-(theta[2]^theta[1]+t)^(1/theta[1])+theta[2])})
    phi.inverse(-log(U)/Y,theta)
  }

#####################################################################################

rcopula.gumbel <- function(n, theta, d)
  #special call to rAC for backwards compatibility
{
  rAC("gumbel",n,d,theta)
}

####################################################################################

rcopula.clayton <- function(n, theta, d)
  #special call to rAC for backwards compatibility
{
  rAC("clayton",n,d,theta)
}

####################################################################################

rcopula.frank <- function(n, theta, d)
  #special call to rAC for backwards compatibility
{
  rAC("frank",n,d,theta)
}

#######################################################################################

rstable <- function(n,alpha,beta=1){
  t0 <- atan(beta*tan((pi * alpha)/2))/alpha
  Theta <- pi * (runif(n)-0.5)
  W <-  - log(runif(n))
  term1 <- sin(alpha*(t0+Theta))/(cos(alpha*t0)*cos(Theta))^(1/alpha)
  term2 <- ((cos(alpha*t0+(alpha-1)*Theta))/W)^((1-alpha)/alpha)
  term1*term2
}

#######################################################################

rFrankMix <- function(n,theta)
{
# SU 05/25/2006: made changes for consistency with R rather than S-Plus.
# The C() method call is void frank(long* n, double *theta, double *output)
# We have passed in the length of the vector to use and the value of theta.
# We will return a vector of length n.
# create an internal vector of length n with all 0s which we will modify and return
    vector = rep(0,n) 
                       
# Note that the function arguments in a .C() call are a LIST.  That means the 
# $ operator can be used to specify which arguments to return.  By default, all the 
# arguments will be returned in a list.  Hence we could call
#   tmp <- .C(name = as.character("frank"),
#             length = as.integer(n),
#             theta = as.double(theta),
#             output = as.double(vector),
#             PACKAGE = "QRMlib")
# Note the parameters are in an 'input-output' type of framework.  E.g. we have
# 'output = as.double(vector)' where the LHS 'output' is the output and the RHS
# 'as.double(vector)' is the input.
# The return value from this function is tmp which is a LIST.  It would contain
# tmp$name, tmp$length, tmp$theta, tmp$output, tmp$PACKAGE which are the named
# outputs from the function (the LHS of all the argument 'equations').
# On some occasions, you may want to return only one of the parameters as an argument.
# Suppose we want only to return vector. Then we write 'res = as.double(vector)' implyin
# 'res' is the output from the function which inputs 'as.double(vector)".  We then indicate
# we want 'tmp' list to contain only the 'res' element by placing a $res AFTER the end of
# the function.                
    tmp <- .C("frank",
            as.integer(n),
            as.double(theta),
            res= as.double(vector),  
#           Also we are concerned only with the 'vector' parameter which we want to denote as 
#           the 'res' element of the tmp list. We indicate we care only about the single output
#           parameter (renamed to 'res') by using the $res notation following the function call.
            PACKAGE="QRMlib")$res
    return(tmp)

#   The old S+ way with its arguments follow.  These do not work in R. 
#   tmp <- .C("frank",
#            n = n,
#            theta = theta,
#            output = rep(0,n),
#            CLASSES = c("long", "double", "double"))
#    tmp$output

}

########################################################################

rBB9Mix <- function(n,theta){
  out <- rep(NA,n)
  for (i in 1:n){
    X <- 0
    U <- 2
    while (U > exp(-X*theta[2]^theta[1])){
      X <- rstable(1,1/theta[1])
      U <- runif(1)
    }
    out[i] <- X
  }
  #If no return(), the last item evaluated is returned. We could use return(out).
  out
}

###############################################################################

rcopula.AGumbel <- function(n,theta,alpha=rep(1,d),d=2)
{
  d <- length(alpha)
  U <- rcopula.gumbel(n=n,theta=theta,d=d)
  U2 <- rcopula.gumbel(n=n,theta=1,d=d)
  alphamat <- matrix(alpha,ncol=d,nrow=n,byrow=TRUE)
  U <- U^(1/alphamat)
  U2 <- U2^(1/(1-alphamat))
  pmax(U,U2)
}

######################################################################################

rcopula.Gumbel2Gp <- function(n = 1000, gpsizes =c(2,2), theta =c(2,3,5))
{
  Y <- rstable(n,1/theta[1])*(cos(pi/(2*theta[1])))^theta[1]
  innerU1 <- rcopula.gumbel(n,theta[2]/theta[1],gpsizes[1])
  innerU2 <- rcopula.gumbel(n,theta[3]/theta[1],gpsizes[2])
  U <- cbind(innerU1,innerU2)
  Y <- matrix(Y, nrow = n, ncol = sum(gpsizes))                               
  out <- exp( - ( - log(U)/Y)^(1/theta[1]))
  #If no return(), the last item evaluated is returned. We could use return(out).
  out
}

####################################################################################

rcopula.GumbelNested <- function(n, theta)
{
  d <- length(theta)+1
  if (d==2)
    out <- rcopula.gumbel(n,theta,d)
  else if (d>2){
    Y <- rstable(n,1/theta[1])*(cos(pi/(2*theta[1])))^theta[1]
    U <- runif(n)
    innerU <- rcopula.GumbelNested(n,theta[-1]/theta[1])
    U <- cbind(U,innerU)
    Y <- matrix(Y, nrow = n, ncol = d)                               
    out <- exp( - ( - log(U)/Y)^(1/theta[1]))
  }
  #If no return(), the last item evaluated is returned. We could use return(out).
  out
}

########################################################################################

########Copulas Densities
  
############################################################################################

dcopula.gauss <- function(u,P,logvalue=FALSE)
{
  d <- dim(u)[2]
  Qdata <- apply(u,2,qnorm)
  out <- dmnorm(Qdata,rep(0,d),P,logvalue=TRUE) - apply(log(apply(Qdata,2,dnorm)),1,sum)
  if(!(logvalue))
    out <- exp(out)
  #If no return(), the last item evaluated is returned. We could use return(out).
  out
}

#####################################################################################

dcopula.t <- function(u,nu,P,logvalue=FALSE)
{
  d <- dim(u)[2]
  Qdata <- apply(u,2,qt,df=nu)
  out <- dmt(Qdata,nu,rep(0,d),P,logvalue=TRUE) - apply(log(apply(Qdata,2,dt,df=nu)),1,sum)
  if(!(logvalue))
    out <- exp(out)
  #If no return(), the last item evaluated is returned. We could use return(out).
  out
}

###############################################################################

#### Copula Fitting

#########################################################################################
fit.gausscopula <- function(Udata)
{
  gausscopula.negloglik <- function(theta)
  {
    P <- Pconstruct(theta)
    -sum(dcopula.gauss(data.nlmin,P,logvalue=TRUE))
  }
  #SU: 06/06/2006. All assign(...,frame=1) calls from S-Plus must be replaced by either
  #assign(...,env=parent.frame()) or by assign(...,env= .GlobalEnv)
  #The original S-PLUS line is commented out:
  #assign("data.nlmin",Udata,frame=1)
  #replaced by the R-code line
  assign("data.nlmin",Udata, env = .GlobalEnv)  #env = parent.frame())

  d <- dim(Udata)[2]
  theta <- Pdeconstruct(Spearman(data.nlmin))  #Udata))
  # SU: 06/07/2006: The S-PLUS calls:
  #out <- nlmin(gausscopula.negloglik,theta)
  #theta <- out$x
  #must be replaced by either the R-language calls nlm() or optim().  nlm() has failed in
  #several cases in the code so I am using optim():
  #out <- nlm(gausscopula.negloglik,theta)
  #theta <- out$estimate 
  optimout <- optim(theta, gausscopula.negloglik, method="BFGS")
  theta <- optimout$par 
  P <- Pconstruct(theta)
  #The S-Plus parameter in the following must be replaced: 
  #list(P=P,converged=out$converged,ll.max=-gausscopula.negloglik(theta))
  #by the following:
  if(optimout$convergence == 0) #if code is 0 we are OK; 
    converged <- TRUE
  else
    converged <- FALSE
  ll.max=-gausscopula.negloglik(theta)
  #7/18/2006: removed all 'global' variables 'assign()'ed in this routine:
  if(!is.null(data.nlmin))
     rm(data.nlmin, envir = .GlobalEnv)
  #Replace list commented out above:
  list(P=P,converged=converged,ll.max=ll.max)
}

###############################################################################################

fit.tcopula <- function(Udata)
{
  tcopula.negloglik <- function(theta)
  {
    nu <- theta[1]  
    P <- Pconstruct(theta[-1])
    -sum(dcopula.t(data.nlmin,nu,P,logvalue=TRUE))
  }
  #SU: 06/06/2006. All assign(...,frame=1) calls from S-Plus must be replaced by either
  #assign(...,env=parent.frame()) or by assign(...,env= .GlobalEnv)
  #The original S-PLUS line is commented out:
  #assign("data.nlmin",Udata,frame=1)
  #replaced by the R-code line
  assign("data.nlmin",Udata,env = parent.frame())

  d <- dim(Udata)[2]
  theta <- c(5,Pdeconstruct(Spearman(data.nlmin)))  #Udata)))

  #SU: 06/07/2006. All S-Plus nlmin(f,x) calls must be replaced by calls to either
  # R-function nlm(f,p) or R-function optim(p,f,method="BFGS")
  # Additionally, the return variables from the R-call out <- nlm() come as a list "out" 
  #and are out$estimate, out$code,and out$minimum (the function value at the estimate).  
  #These must replace the S-Plus outputs out$x, out$converged (there is no S-Plus analog 
  #for the minimum function value). Since nlm() has failed to converge in
  #several cases in the code, I am using optim().

  # The S-PLUS calls:
  #out <- nlmin(tcopula.negloglik,theta,max.it=500,max.fcal=500,verbose=T)
  #theta <- out$x
  #must be replaced by the R-language code:
  optimout <- optim(theta,tcopula.negloglik,method="BFGS")
  theta <- optimout$par

  nu <- theta[1]
  P <- Pconstruct(theta[-1])
  #S-Plus output out$converged must be replaced in list() command:
  #list(P=P,nu=nu,converged=out$converged,ll.max=-tcopula.negloglik(theta))
  if(optimout$convergence == 0) #if code is 0 we are OK; 
    converged <- TRUE
  else
    converged <- FALSE

  ll.max=-tcopula.negloglik(theta)
  #7/18/2006: removed all 'global' variables 'assign()'ed in this routine:
  if(!is.null(data.nlmin))
     rm(data.nlmin, envir = .GlobalEnv)
  #Replace list commented out above:
  list(P=P,nu=nu,converged=converged,ll.max=ll.max)
}

################################################################################################
# The default method is assigned as Kendall if the argument method is not passed
fit.tcopula.rank <- function(Udata,method="Kendall")
{
  if (method=="Kendall")
  {
    Rtau <- Kendall(Udata)
    P <- sin(pi*Rtau/2)
  }
  else
    P <- Spearman(Udata)
  if (min(eigen(P)$values)<0) stop ("Non p.s.d. covariance matrix")
  tcopula.negloglik <- function(nu){
    -sum(dcopula.t(data.nlmin,nu,P.fix,logvalue=TRUE))
  }
  #SU: 06/06/2006. All assign(...,frame=1) calls from S-Plus must be replaced by either
  #assign(...,env=parent.frame()) or by assign(...,env= .GlobalEnv)
  #The original S-PLUS lines are commented out:
  #assign("P.fix",P,frame=1)
  #assign("data.nlmin",Udata,frame=1)
  assign("P.fix", P, env = parent.frame())
  assign("data.nlmin", Udata, env = parent.frame())
  nu <- 5
  ##SU: 06/07/2006. All S-Plus nlmin(f,x) calls must be replaced by calls to either
  # R-function nlm(f,p) or R-function optim(p,f,method="BFGS")
  # Additionally, the return variables from the R-call out <- nlm() come as a list "out" 
  #and are out$estimate, out$code,and out$minimum (the function value at the estimate).  
  #These must replace the S-Plus outputs out$x, out$converged (there is no S-Plus analog 
  #for the minimum function value. Since nlm() has failed to converge in
  #several cases in the code, I am using optim().

  # 06/08/2006: The S-PLUS calls:
  #out <- nlmin(tcopula.negloglik,nu,max.it=500,max.fcal=500,verbose=T)
  #nu <- out$x
  #list(P=P,nu=nu,converged=out$converged,ll.max=-tcopula.negloglik(nu))
  # are replaced by the following R-language code:
  optimout <- optim(nu, tcopula.negloglik,method="BFGS")
  nu <- optimout$par
  if(optimout$convergence == 0) #if code is 0 we are OK; 
    converged <- TRUE
  else
    converged <- FALSE
  ll.max=-tcopula.negloglik(nu)
  #Replace list commented out above:
  list(P=P,nu=nu,converged=converged,ll.max=ll.max)

}

###############################################

dcopula.clayton <- function(u,theta,logvalue=FALSE)
{
  d <- dim(u)[2]
  if (d>2) stop("Clayton copula density only implemented for d=2")
  u1 <- u[,1]
  u2 <- u[,2]
  out <- log(1 + theta) + (-1 - theta)*log(u1*u2) + (-2 - 1/theta)*log(u1^( - theta) + u2^( - theta) - 1)
  if (!(logvalue))
    out <- exp(out)
  #If no return(), the last item evaluated is returned. We could use return(out).
  out
}

##############################################

dcopula.gumbel <- function(u,theta,logvalue=FALSE)
{
  d <- dim(u)[2]
  if (d>2) stop("Gumbel copula density only implemented for d=2")
  u1 <- u[,1]
  u2 <- u[,2]
  innerfunc <- function(x,y,theta){((-log(x))^theta + (-log(y))^theta)^(1/theta)}
  out <- -innerfunc(u1,u2,theta) -log(u1*u2) + (theta-1)*log(log(u1)*log(u2))+ log(theta-1+innerfunc(u1,u2,theta)) +(1-2*theta)*log(innerfunc(u1,u2,theta))
  if (!(logvalue))
    out <- exp(out)
  #If no return(), the last item evaluated is returned. We could use return(out).
  out
}

############################################

fit.Archcopula2d <- function(data,name)
{
  if (name=="clayton")
    negloglik <- function(x)
      {
       -sum(dcopula.clayton(data.tmp,x,logvalue=TRUE))
      }
  else if (name=="gumbel")
    negloglik <- function(x)
      {
       -sum(dcopula.gumbel(data.tmp,x,logvalue=TRUE))
      }
  else stop ("copula not implemented")
  #SU: 06/06/2006. All assign(...,frame=1) calls from S-Plus must be replaced by either
  #assign(...,env=parent.frame()) or by assign(...,env= .GlobalEnv)
  #assign("data.tmp", data, frame = 1)
  assign("data.tmp", data, env = parent.frame())

  #SU: 06/07/2006. All S-Plus nlmin(f,x) calls must be replaced by calls to R-function nlm(f,p)
  # Additionally, the return variables from the R-call out <- nlm() come as a list "out" 
  #and are out$estimate, out$code,and out$minimum (the function value at the estimate).  
  #These must replace the S-Plus outputs out$x, out$converged (there is no S-Plus analog 
  #for the minimum function value). Because nlm() failed many times, I am using optim().
  # The S-PLUS calls:
  #fit <- nlmin(negloglik, x = 2)
  #theta <- fit$x
  #converged <- fit$converged
  # are replaced by the following for R:
  #note we do not use p=2 as argument in function invocation due to R lazy argument 
  #implementation.  See pp. 22-3 in R-lang.pdf (section 4.3.3 on argument evaluation)
  fit <- optim(2, negloglik, method="BFGS")
  theta <- fit$par
  
  if(fit$convergence == 0) #if convergence is 0 we are OK
     converged <- TRUE
  else
    converged <- FALSE


  # Function hessb() is called from funtionsUtility.R
  hessianmatrix <- hessb(negloglik,theta)
  varcov <- solve(hessianmatrix)
  ses <- sqrt(diag(varcov))
  ll.max <- -negloglik(theta)
  #See pp. 22-3 in R-lang.pdf (section 4.3.3 on argument evaluation. Due to lazy argument
  #evaluation in R, we don't use z <- foo(x=y)in implementation but merely pass parameters
  out <- list(ll.max=ll.max, theta=theta, ses=ses, converged=converged)
  #If no return(), the last item evaluated is returned. We could use return(out).
  out
}

##########################################################

Pconstruct <- function(theta){
  n <- length(theta)
  d <- (1 + sqrt(1+8*n))/2
  A <- matrix(0,nrow=d,ncol=d)
  A[lower.tri(A)] <- theta
  diag(A) <- rep(1,d)
  Q <- A %*% t(A)
  P <- CovToCor(Q)
  P
}

###########################################################################################
Pdeconstruct <- function(P){
  A <- t(chol(P))
  Adiag <- diag(diag(A))
  Astar <- solve(Adiag) %*% A
  Astar[lower.tri(Astar)]
}
