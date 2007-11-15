# QRMlib: version 1.4.2
# this file is a component of QRMlib 

# Copyright (C) 2005-07 Alexander McNeil
# R-language translation/additions Copyright (C) 2006-2007 Scott Ulman

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

# Contact: Alexander J. McNeil:  A.J.McNeil@hw.ac.uk
# R-language contact: Scott Ulman : scottulman@hotmail.com 
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
#Simulate Archimedean Copulas.   It may be called instead of functions rcopula.clayton()
#rcopula.gumbel(), rcopula.frank().  You may also simulate BB9 and GIG copulas which
#have no direct calls of their own. Added type GIG since version 1.4. 
rAC <- function(name, n, d, theta){
#generic function for Archimedean copula simulation
  illegalpar <- switch(name,
                       clayton = (theta < 0),
                       gumbel = (theta < 1),
                       frank = (theta < 0),
                       BB9 = ((theta[1] < 1) | (theta[2] < 0)),
                       GIG = ((theta[2] < 0) | (theta[3] < 0) | ((theta[1]>0) & (theta[3]==0)) | ((theta[1]<0) & (theta[2]==0))))
  if(illegalpar)
    stop("Illegal parameter value")
  independence <- switch(name,
                         clayton = (theta == 0),
                         gumbel = (theta == 1),
                         frank = (theta == 0),
                         BB9 = (theta[1] == 1),
                         GIG=FALSE)
  U <- runif(n * d)
  U <- matrix(U, nrow = n, ncol = d)
  if(independence)
    return(U)
  Y <- switch(name,
              clayton = rgamma(n, 1/theta),
              gumbel = rstable(n, 1/theta) * (cos(pi/(2 * theta)))^theta,
              frank = rFrankMix(n, theta),
              BB9 = rBB9Mix(n, theta),
              GIG = rGIG(n,theta[1],theta[2],theta[3]))
  Y <- matrix(Y, nrow = n, ncol = d)
  phi.inverse <- switch(name,
                        clayton = function(t, theta)
                        {
			      (1 + t)^(-1/theta)
                        }
                        ,
                        gumbel = function(t, theta)
                        {
			      exp( - t^(1/theta))
                        }
                        ,
                        frank = function(t, theta)
                        {
                          (-1/theta) * log(1 - (1 - exp( - theta)) * exp( - t))
                        }
                        ,
                        BB9 = function(t, theta)
                        {
                          exp( - (theta[2]^theta[1] + t)^(1/theta[1]) + theta[2])
                        }
                        ,
                        GIG = function(t, theta)
                        {
                          lambda <- theta[1]
                          chi <- theta[2]
                          psi <- theta[3]
                          if (chi==0)
                            out <- (1+2*t/psi)^(-lambda)
                          else if (psi==0)
                            out <- 2^(lambda+1)*exp(besselM3(lambda,sqrt(2*chi*t),logvalue=TRUE)-lambda*log(2*chi*t)/2)/gamma(-lambda)
                          else
                            out <- exp(besselM3(lambda,sqrt(chi*(psi+2*t)),logvalue=TRUE)+lambda*log(chi*psi)/2-besselM3(lambda,sqrt(chi*psi),logvalue=TRUE)-lambda*log(chi*(psi+2*t))/2)
                          out
                        }
                        )
  phi.inverse( - log(U)/Y, theta)
}

###################################################################################
#Simulate weighted Archimedean copulas:
rACp <- function(n,d,theta,A,name="gumbel")
{
  p <- length(theta)
  if ((dim(A)[1] != d) | (dim(A)[2] !=p)) stop("Weight matrix A has incorrect dimensions")
  sumcheck <- apply(A,1,sum) - rep(1,d)
  if (sum(sumcheck^2) != 0) stop("Weights do not sum to one")
  for (j in 1:p)
    {
      tmp <- rAC(name=name,n=n,d=d,theta=theta[j])
      Amat <- matrix(A[,j],ncol=d,nrow=n,byrow=TRUE)
      tmp <- tmp^(1/Amat)
      eval(parse(text=paste("U",j," <- tmp",sep="")))
    }
  args <- paste("U",1:p,sep="",collapse=",")
  result <- parse(text=paste("pmax(",args,")"))
  eval(result)
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
#  #10/12/2007: separate methods for SPlus and R:
#  if (exists("is.R") && is.function(is.R) && is.R()) 
#  {
    #In R-2.6.0, we must generate a vector of output values to pass to the 
    #C-frank function or we will receive the following error message:
    #"Error in as.double(vector) : cannot coerce to vector"
    output <- rep(0,n); 
    tmp <- .C("frank",
            as.integer(n),
            as.double(theta),
            res= as.double(output), #this variable will be changed in the C-call 
#           we are concerned only with returning the 'vector' parameter 'res'
#           which we denote as the 'res' element. This tells R to return only
#           tmp$res as the value of tmp. 
            PACKAGE="QRMlib")$res;
    return(tmp);
#  } else
#  {
#    tmp <- .C("frank",
#            n = n,
#            theta = theta,
#            output = rep(0,n),
#            CLASSES = c("long", "double", "double"))
#    tmp$output
#  }
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
  out
}

###############################################################################
#This function has been deprecated and replaced by rAC()
# rcopula.AGumbel.old <- function(n,theta,alpha=rep(1,d),d=2)
# {
#  d <- length(alpha)
#  U <- rcopula.gumbel(n=n,theta=theta,d=d)
#  U2 <- rcopula.gumbel(n=n,theta=1,d=d)
#  alphamat <- matrix(alpha,ncol=d,nrow=n,byrow=TRUE)
#  U <- U^(1/alphamat)
#  U2 <- U2^(1/(1-alphamat))
#  pmax(U,U2)
#}


######################################################################################

rcopula.Gumbel2Gp <- function(n = 1000, gpsizes =c(2,2), theta =c(2,3,5))
{
  Y <- rstable(n,1/theta[1])*(cos(pi/(2*theta[1])))^theta[1]
  innerU1 <- rcopula.gumbel(n,theta[2]/theta[1],gpsizes[1])
  innerU2 <- rcopula.gumbel(n,theta[3]/theta[1],gpsizes[2])
  U <- cbind(innerU1,innerU2)
  Y <- matrix(Y, nrow = n, ncol = sum(gpsizes))                               
  out <- exp( - ( - log(U)/Y)^(1/theta[1]))
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
  out
}

########################################################################################

########Copulas Densities
  
############################################################################################

dcopula.gauss <- function(u,P,logvalue=FALSE){
  d <- dim(u)[2]
  Qdata <- apply(u,2,qnorm)
  out <- dmnorm(Qdata,rep(0,d),P,logvalue=TRUE) - apply(log(apply(Qdata,2,dnorm)),1,sum)
  if(!(logvalue))
    out <- exp(out)
  out
}

#####################################################################################

dcopula.t <- function(u,nu,P,logvalue=FALSE){
  d <- dim(u)[2]
  Qdata <- apply(u,2,qt,df=nu)
  out <- dmt(Qdata,nu,rep(0,d),P,logvalue=TRUE) - apply(log(apply(Qdata,2,dt,df=nu)),1,sum)
  if(!(logvalue))
    out <- exp(out)
  out
}

###############################################################################

#### Copula Fitting

#########################################################################################

fit.gausscopula <- function(Udata)
{
  gausscopula.negloglik <- function(theta,data)
  {
    P <- Pconstruct(theta)
    -sum(dcopula.gauss(data,P,logvalue=TRUE))
  }
  theta <- Pdeconstruct(Spearman(Udata))
  out <- nlminb(theta,gausscopula.negloglik,data=Udata)
  #10/12/2007: different parameter returns for SPlus and R:
  if (exists("is.R") && is.function(is.R) && is.R()) 
  {
     theta <- out$par;
     P <- Pconstruct(theta);
     if(out$convergence == 0) conv=TRUE
     else                     conv=FALSE;
     list(P=P,converged=conv,ll.max=-out$objective);
  } else
  {
    theta <- out$parameters
    P <- Pconstruct(theta)
    list(P=P,converged=out$message,ll.max=-out$objective)
  }
}

###############################################################################################

fit.tcopula <- function(Udata)
{
  tcopula.negloglik <- function(theta,data)
  {
    nu <- theta[1]  
    P <- Pconstruct(theta[-1])
    -sum(dcopula.t(data,nu,P,logvalue=TRUE))
  }
  theta <- c(5,Pdeconstruct(Spearman(Udata)))
  out <- nlminb(theta,tcopula.negloglik,data=Udata)
  #10/12/2007: different parameter returns for SPlus and R:
  if (exists("is.R") && is.function(is.R) && is.R()) 
  {
     theta <- out$par;
     nu <- theta[1];
     P <- Pconstruct(theta[-1]);
     if(out$convergence == 0) conv=TRUE
     else                     conv=FALSE;
     list(P=P,nu=nu,converged=conv,ll.max=-out$objective);
  } else
  {
     theta <- out$parameters
     nu <- theta[1]
     P <- Pconstruct(theta[-1])
     list(P=P,nu=nu,converged=out$message,ll.max=-out$objective)
  }
}

################################################################################################

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
  tcopula.negloglik <- function(nu,data,Pest){
    -sum(dcopula.t(data,nu,Pest,logvalue=TRUE))
  }
  nu <- 5
  out <- nlminb(nu,tcopula.negloglik,data=Udata,Pest=P)
  #10/12/2007: different parameter returns for SPlus and R:
  if (exists("is.R") && is.function(is.R) && is.R()) 
  {
     nu <- out$par;
     if(out$convergence == 0) conv=TRUE
     else                     conv=FALSE;
     list(P=P,nu=nu,converged=conv,ll.max=-out$objective)
  } else
  {
     nu <- out$parameters
     list(P=P,nu=nu,converged=out$message,ll.max=-out$objective)
  }
}

###############################################
#This method supplants the old dcopula.gumbel and dcopula.clayton functions
dcopula.AC <- function(u,theta,name,logvalue=TRUE){
  d <- dim(u)[2]
  if ((name=="gumbel") & (d>2)) stop("Only bivariate Gumbel implemented")
  illegalpar <- switch(name,
		clayton = (theta <= 0),
		gumbel = (theta <= 1))
  if (illegalpar) out <- NA
  else{
    phi <- switch(name,
                  clayton = function(u,theta){(u^(-theta)-1)/theta},
                  gumbel = function(u,theta){(-log(u))^theta})
    lnegphidash <-switch(name,
                         clayton = function(u,theta){(-theta-1)*log(u)},
                         gumbel = function(u,theta){log(theta)+(theta-1)*log(-log(u))-log(u)})
    loggfunc <- switch(name,
                       clayton = function(t,d,theta){d*log(theta)+sum(log((1:d)+1/theta-1))-(d+1/theta)*log(t*theta+1)},
                       gumbel = function(t,d=2,theta){-2*log(theta) -t^(1/theta) + (1/theta-2) *log(t) + log(t^(1/theta)+theta-1)})
    gu <- apply(phi(u,theta),1,sum)
    term1 <- loggfunc(gu,d,theta)
    term2 <- apply(lnegphidash(u,theta),1,sum)
    out <- term1+term2
    if(!(logvalue))
      out <- exp(out)
  }
  out
}

#########################################################
#This method supplants the old fit.Archcopula2d() method of version 1.4
fit.AC <- function(Udata, name="clayton")
{
  negloglik <- function(x,data,name)
    {
      - sum(dcopula.AC(data, x, name, logvalue = TRUE))
    }
  lb <- switch(name,gumbel=1+1e-10,clayton=0)
  fit <- nlminb(2, negloglik, lower=lb, data=Udata, name=name)
  if (exists("is.R") && is.function(is.R) && is.R()) 
  {
    theta <- fit$par;
    if(fit$convergence == 0) converged=TRUE
    else                     converged=FALSE;
  } else
  {
    theta <- fit$parameters
    converged <- fit$message
  }
  hessianmatrix <- hessb(negloglik, theta, data=Udata, name=name);
  varcov <- solve(hessianmatrix);
  se <- sqrt(diag(varcov));
  ll.max <-  - fit$objective;
  out <- list(ll.max = ll.max, theta = theta, se=se, converged = 
              converged);
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

###############################################

dcopula.clayton <- function(u,theta,logvalue=FALSE){
  d <- dim(u)[2]
  if (d>2) stop("Clayton copula density only implemented for d=2")
  u1 <- u[,1]
  u2 <- u[,2]
  out <- log(1 + theta) + (-1 - theta)*log(u1*u2) + (-2 - 1/theta)*log(u1^( - theta) + u2^( - theta) - 1)
  if (!(logvalue))
    out <- exp(out)
  out
}

##############################################

dcopula.gumbel <- function(u,theta,logvalue=FALSE){
  d <- dim(u)[2]
  if (d>2) stop("Gumbel copula density only implemented for d=2")
  u1 <- u[,1]
  u2 <- u[,2]
  innerfunc <- function(x,y,theta){((-log(x))^theta + (-log(y))^theta)^(1/theta)}
  out <- -innerfunc(u1,u2,theta) -log(u1*u2) + (theta-1)*log(log(u1)*log(u2))+ log(theta-1+innerfunc(u1,u2,theta)) +(1-2*theta)*log(innerfunc(u1,u2,theta))
  if (!(logvalue))
    out <- exp(out)
  out
}

############################################
#This function has been deprecated by the new function fit.AC()
fit.Archcopula2d <- function(Udata,name)
{
  if (name=="clayton")
    negloglik <- function(x,data)
      {
       -sum(dcopula.clayton(data,x,logvalue=TRUE))
      }
  else if (name=="gumbel")
    negloglik <- function(x,data)
      {
       -sum(dcopula.gumbel(data,x,logvalue=TRUE))
      }
  else stop ("copula not implemented")
  lb <- switch(name,gumbel=1+1e-10,clayton=0)
  fit <- nlminb(2,negloglik, lower=lb, data=Udata)
  if (exists("is.R") && is.function(is.R) && is.R()) 
  {
    theta <- fit$par;
    if(fit$convergence == 0) converged=TRUE
    else                     converged=FALSE;
  } else
  {
    theta <- fit$parameters
    converged <- fit$message
  }
  hessianmatrix <- hessb(negloglik,theta,data=Udata)
  varcov <- solve(hessianmatrix)
  ses <- sqrt(diag(varcov))
  ll.max <- -fit$objective
  out <- list(ll.max=ll.max, theta = theta, ses=ses, converged = converged)
  out
}
