# QRMlib: version 1.4.2
# this file is a component of QRMlib 

# Copyright (C) 2005-07 Alexander McNeil
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

# Contact: Alexander J. McNeil:  a.j.mcneil@hw.ac.uk 
# Contact for R-language translation: Scott Ulman: scottulman@hotmail.com
# Note in the R-translations that TRUE has been substituted throughout the 
# code for T (and FALSE for F) when setting default parameter values as 
# specified in section 2.4 of an "Introduction to R" (R-intro.pdf).
# Otherwise "R CMD check -QRMlib" returned the following type of error when
# trying to run an example:
# > BiDensPlot(func=dmnorm,mu=c(0,0),Sigma=equicorr(2,-0.7))
#   Error in func(data, ...) : F used instead of FALSE
#   Execution halted
##########################################################################
#This module contains functions to build random samples and distributions for NORMAL MIXTURE models
#as described in Chapter 3 of QRM (especially section 3.2)
 
#Generate a random sample from the GIG (Generalized Inverse Gaussian) distribution.
# See p. 77 of QRM text.  Note this function calls a C-language function from QRSim.c
# to do the bulk of the work.  This method was evidently suggested by Atkinson in 1982.
# This function generates a random series for the Generalized Inverse Gaussian family.
# When this function is called from other functions in the QRMlib or from chapter scripts,
# its return value is assumed to be a vector of simulated values.
#This function was rewritten by SU on 11/2/2007 to eliminate all assign() statements
rGIG <- function(n, lambda, chi, psi, envplot = FALSE, messages = FALSE)
{
   if((chi < 0) | (psi < 0))
     stop("Invalid parameters for GIG with chi, psi both negative");
   if((chi == 0) & (lambda <= 0))
     stop("Invalid parameters for GIG with chi=0 and lambda <=0");
   if((psi == 0) & (lambda >= 0))
     stop("Invalid parameters for GIG with psi=0 and lambda>=0");
   if((chi == 0) & (lambda > 0))
     return(rgamma(n, shape = lambda, rate = (psi/2)));
   if((psi == 0) & (lambda < 0))
     return(1/rgamma(n, shape = ( - lambda), rate = (chi/2)));
   message <- NULL;
   if(abs(lambda) < 1)
     message <- paste(message, "Not necessarily efficient rejection method", "\n");
   neglambda <- FALSE;
   if(lambda < 0) 
   {
     neglambda = TRUE;
     lambda <- abs(lambda);
     tmp <- c(chi, psi);
     #switch values of chi and psi:
     chi <- tmp[2];
     psi <- tmp[1];
   }

   #Local function:
   efunc <- function(x, lambda, chi, psi)
   {
     (x^(lambda - 1)) * exp( - (chi/x + psi * x)/2);
   }
   #Local function:
   calcmode <- function(lambda, chi, psi)
   {
     if(psi > 0)
       return(((lambda - 1) + sqrt(chi * psi + (1 - lambda)^2))/psi)
     else if(psi == 0)
	return(chi/(2 - 2 * lambda))
     else stop("Problem in mode function")
   }

   #Local function (to be called only when lambda < 1 so only one theta parameter with which to optimize
   objectiveOneParam <- function(theta,lambda, chi, psi, themode)
   {
       if(lambda >= 1) stop("lambda must be less than 1 to call objectiveOneParam\n");
              
       if(length(theta)!= 1) stop("theta must be a single number in objectiveOneParam\n");
       if(theta <= 0)
         out <- NA
       else 
       {
         Delta1 <- (exp(themode * theta) - 1)/theta;
         Delta2 <- (2 * exp(( - themode * psi)/2))/psi;
         xL <- calcmode(lambda, chi, psi + 2 * theta);
         xH <- chi/(2 - 2 * lambda);
         #Call nested local functions:
         S1 <- efunc(xL, lambda, chi, psi + 2 * theta);
         S2 <- efunc(xH, lambda, chi, 0);
         out <- Delta1 * S1 + Delta2 * S2;
       } #end else
       out;
   } # end objectiveOneParam()

   #Local function (to be called only when lambda >= 1 so two theta parameters exist to optimize
   objectiveTwoParams <- function(theta, lambda, chi, psi, themode)
   {
       if(lambda < 1) stop("lambda must be at least equal to 1 to call objectiveTwoParam\n");
       if(length(theta)!= 2) stop("theta must be vector of length 2 in objectiveTwoParams\n");

       if((theta[1] <= 0) | (theta[2] <= 0))
         out <- NA
       else if((psi - 2 * theta[2]) < 0)
	  out <- NA
       else 
       {
         Delta1 <- (exp(themode * theta[1]) - 1)/theta[1];
         Delta2 <- exp( - themode * theta[2])/theta[2];
         xL <- calcmode(lambda, chi, psi + 2 * theta[1]);
         xH <- calcmode(lambda, chi, psi - 2 * theta[2]);
         S1 <- efunc(xL, lambda, chi, psi + 2 * theta[1]);
         S2 <- efunc(xH, lambda, chi, psi - 2 * theta[2]);
         out <- Delta1 * S1 + Delta2 * S2;
       } #end else
       out;
  } #end objectiveTwoParams() 

  #Call internal function:
  themode <- calcmode(lambda, chi, psi);

  #Initialize theta depending upon lambda value and call optim() using one of the 
  #two appropriate objective functions:
  if(lambda < 1) 
  {
     #there is only one parameter theta and we minimize relative to it: 
     lambdaLOW <- TRUE;
     theta <- 0.01;
     optimout <- optim(c(0.25), objectiveOneParam, gr=NULL, lambda=lambda,chi=chi,
         psi=psi, themode=themode, method="BFGS");
     spar <- optimout$par[1];
     ppar <- psi/2;
  }
  else
  {
     #there are two parameters (the spar and the ppar) and we minimize with respect to them.
     # spar will be extracted from theta[1] and ppar will be extracted from theta[2]
     #initialize theta[1] and theta[2]
     lambdaLOW <- FALSE;
     theta <- c(0.01, psi/4);
     optimout <- optim(c(.01,.25), objectiveTwoParams, gr=NULL, lambda=lambda,chi=chi,
         psi=psi, themode=themode, method="BFGS");
     spar <- optimout$par[1];
     ppar <- optimout$par[2];
   }

   if(optimout$convergence != 0)
   { 
     message <- paste(optimout$message, 
           "Problems finding optimal s and p (use option envplot for reassurance)","\n");
     print(message);
   } #end if no convergence
  
   #Calculate new values for xL, xH, etc using the parameters spar and ppar estimated
   xL <- calcmode(lambda, chi, psi + 2 * spar);
   xH <- calcmode(lambda, chi, psi - 2 * ppar);
   S1 <- efunc(xL, lambda, chi, psi + 2 * spar);
   S2 <- efunc(xH, lambda, chi, psi - 2 * ppar);
   Delta1 <- (exp(themode * spar) - 1)/spar;
   Delta2 <- exp( - themode * ppar)/ppar;
   k <- 1/((Delta1/S2) + (Delta2/S1));
   k1 <- k/S2;
   k2 <- k/S1;
   rpar <- k1 * Delta1;
   if(envplot) 
   {
     xdat <- seq(from = 0.01, to = themode * 20, length = 1000);
     envelope2 <- (xdat <= themode) * exp(spar * xdat) * S1 + 
			(xdat > themode) * exp( - ppar * xdat) * S2;
     envelope <- (xdat <= themode) * exp(spar * xdat) * k1 + (xdat >
			themode) * exp( - ppar * xdat) * k2;
     ydat <- efunc(xdat, lambda, chi, psi);
     yr <- range(ydat, envelope, envelope2);
     plot(xdat, ydat, ylim = yr, type = "l");
     abline(v = themode);
     lines(xdat, envelope, col = 2);
     lines(xdat, envelope2, lty = 2, col = 2);
   } #end if envplot
   
   #Initialize the xsim vector to all 0s.  
   xsim <- rep(0, n);
   #Initialize new to n:
   new = n;

   # SU 06/02/2006: made changes for consistency with R rather than S-Plus
   # The method call in C-language is 
   # void rgig(long *n, double *r, double *s, double *p, double *k1, double *k2, 
   # double *lambda, double *chi, double *psi, double *s1, double *s2, double *xsim)
   # This C-language function resets the values of both n (the 1st parameter) and xsim (the last parameter).
   # Note that when a .C() call is invoked in R or SPlus, ALL ARGUMENTS to the function are returned
   # as a LIST UNLESS you set one of the parameters up as a result by using notation like
   # res=as.integer(n) and then follow the function call with with the name of the individual list
   # element to return like ... PACKAGE="QRMlib")$res as the end of the function call.
   tmp <- .C("rgig",
            #the parameter n is changed by the function call:
		new = as.integer(n),  
		as.double(rpar),
		as.double(spar),
		as.double(ppar),
		as.double(k1),
		as.double(k2),
		as.double(lambda),
		as.double(chi),
		as.double(psi),
		as.double(S1), 
		as.double(S2),
            #the parameter xsim is changed by the function call
		xsim = as.double(xsim), 
            # Use the list variable $xsim to indicate we want tmp to be replaced by xsim
            # so we will return only the simulated values rather than the whole list
            PACKAGE="QRMlib") $xsim;

   efficiency <- n/new;
   
   message <- paste(message, "Efficiency", round(efficiency * 100, 1),"\n");
   if(messages)
     cat(message);

   #SU removed all 'assigned' parent environment variables:
   
   if(neglambda)
       return(1/tmp)  #tmp should now equal xsim due to listing $xsim after function call
   else return(tmp);  #tmp now equal to xsim due to listing $xsim after function call
}
###########################################################################

dsmghyp <- function(x, lambda, chi, psi, mu, Sigma, logvalue=FALSE)
{
  if ((psi==0) & (lambda <0)){
    nu <- chi
    out <- dmt(x,nu,mu,Sigma,logvalue)
  }
  else if (psi >0){
    d <- dim(x)[2]
    Q <- mahalanobis(x, mu, Sigma)
    log.top <- besselM3((lambda - d/2), sqrt(psi*(chi + Q)),logvalue=TRUE)
    log.bottom <- (d/2-lambda)*log(sqrt(psi*(chi + Q)))
    if (chi>0){
      log.const.top <- (-lambda/2)*log(psi*chi) +(d/2)*log(psi)
      log.const.bottom <- (d/2)*log(2 * pi) + besselM3(lambda, sqrt(chi * psi),logvalue=TRUE) + 0.5*log(det(Sigma))
    }
    else if (chi==0){
      log.const.top <- d*log(psi)/2 + (1-lambda)*log(2)
      log.const.bottom <- (d/2)*log(2 * pi) +log(gamma(lambda)) + 0.5*log(det(Sigma))
    }
    out <- log.const.top + log.top - log.const.bottom - log.bottom
  }
  if (!(logvalue)) out <- exp(out)
  out
}

#################################################################

dmghyp <- function(x, lambda, chi, psi, mu, Sigma, gamma,logvalue=FALSE)
{
  if (sum(abs(gamma))==0) out <- dsmghyp(x, lambda,chi,psi,mu,Sigma,logvalue=TRUE)
  else if ((psi==0) & (lambda <0)){
    nu <- chi
    d <- dim(x)[2]
    n <- dim(x)[1]
    Q <- mahalanobis(x, mu, Sigma)
    Offset <- t(gamma) %*% solve(Sigma) %*% gamma
    beta <- solve(Sigma) %*% gamma
    mu.matrix <- matrix(mu,nrow=n,ncol=d,byrow=TRUE)
    tilt <- as.vector((x-mu.matrix) %*% beta)
    interm <- sqrt((nu+Q)*Offset)
    log.top <- besselM3((nu+d)/2,interm,logvalue=TRUE) + tilt
    log.bottom <- (nu+d)*log(interm)/2
    log.const.top <- nu*log(nu)/2 +(nu+d)*log(Offset)/2
    log.const.bottom <- (d/2)*log(2*pi) + 0.5*log(det(Sigma)) +(nu/2-1)*log(2)+log(gamma(nu/2))
    out <- log.const.top + log.top - log.const.bottom - log.bottom
  }
  else if (psi >0){
    d <- dim(x)[2]
    n <- dim(x)[1]
    Q <- mahalanobis(x, mu, Sigma)
    Offset <- t(gamma) %*% solve(Sigma) %*% gamma
    beta <- solve(Sigma) %*% gamma
    mu.matrix <- matrix(mu,nrow=n,ncol=d,byrow=TRUE)
    tilt <- as.vector((x-mu.matrix) %*% beta)
    log.top <- besselM3((lambda - d/2), sqrt((psi+Offset)*(chi + Q)),logvalue=TRUE) + tilt
    log.bottom <- (d/2-lambda)*log(sqrt((psi+Offset)*(chi + Q)))
    if (chi>0){
      log.const.top <- (-lambda/2)*log(psi*chi) +(d/2)*log(psi) + (d/2-lambda)*log(1+Offset/psi)
      log.const.bottom <- (d/2)*log(2 * pi) + besselM3(lambda, sqrt(chi * psi),logvalue=TRUE) + 0.5*log(det(Sigma))
      out <- log.const.top + log.top - log.const.bottom - log.bottom
    }
    else if (chi==0){
      log.const.top <- d*log(psi)/2 + (1-lambda)*log(2) + (d/2-lambda)*log(1+Offset/psi)
      log.const.bottom <- (d/2)*log(2 * pi) +log(gamma(lambda)) + 0.5*log(det(Sigma))
      out <- log.const.top + log.top - log.const.bottom - log.bottom
    }
    else out <- NA
  }
  if (!(logvalue)) out <- exp(out)
  out
}

########################################################################

besselM3 <- function(lambda = 9/2, x = 2, logvalue = FALSE)
{
      lambda <- abs(lambda)
	integer.part <- as.integer(lambda)
      #SU: 6/29/2006: In S-Plus, modulo function is defined as as 
      #e1-floor(e1/e2)*e2 if e2!=0 and e1 otherwise (see Knuth, 1968, section 1.2.4)
      #R-help says 'x%%y' indicates 'x mod y' unless 'y == 0' where the result is 'NA' or 'NaN'
      #Hence we must alter the following S-Plus line:
	#remainder <- lambda %% integer.part
      #as follows to make the R-result consistent with the S-Plus result of 'e1 otherwise' where
      #e1 is lambda and e2 is integer.part
      if(integer.part == 0) 
         remainder <- lambda
      else 
        remainder <- lambda %% integer.part
 
	if(lambda == 0)
		intype = "zero"
	else if(remainder == 0)
		intype = "integer"
	else if((remainder == 0.5) & (integer.part == 0))
		intype = "half"
	else if(remainder == 0.5)
		intype = "halfinteger"
	else intype = "fractional"

      if(intype == "zero") 
      {
            # SU 06/08/2006: made changes for consistency with R rather than S-Plus
            # The method call in C-language is 
            # void besselM3z(double *x, long *n, double *err, double *y)
            # The following S-Plus code must be replaced
		#tmp <- .C("besselM3z",
		#	x = x,
		#	nx = length(x),
		#	err = rep(0, length(x)),
		#	y = rep(0, length(x)),
		#	CLASSES = c("double", "integer", "double",
		#		"double"))
		#out <- tmp$y
            #    if (logvalue==T) out <- log(out)
            #by the following R-code:
            nx <- length(x)
            err <- rep(0,nx)
            y <- rep(0,nx)
            tmp <- .C("besselM3z",
			     as.double(x),
			     as.integer(nx),
			     as.double(err),
                             y = as.double(y),
                             PACKAGE="QRMlib");
                  #You must use 'y = as.double(y)' if you want to use tmp$y as a list element.
                  #The LHS (y =) of the notation 'y = as.double(y)' defines an 'output' and 
                  #the RHS defines an 'input'. 
                  #If you want to use solely 'as.double(y)', you have no named 'output' parameter so you
                  #must set 'out <- tmp[[4]]' since 4 is the order-number of the y input parameter.
                  #Parameter orders are 'x,nx,err,y' = (1,2,3,4) so y is 4th parameter in list.
                  #The function name "besselM3z' is parameter 0.
                  #Hence you must refer to x as tmp[[0]], to nx as tmp[[1]] since they are not 'named'.
                  #If you had used 'x = as.double(x)', you could use tmp$x.

		out <- tmp$y;
             if (logvalue==TRUE) out <- log(out)
	}
	if(intype == "integer") {
            # SU 06/08/2006: made changes for consistency with R rather than S-Plus
            # The method call in C-language is 
            #void besselM3(long *lambda, double *x, long *n, double *err, double *y, long *logvalue)
		#tmp <- .C("besselM3",
		#	lambda = lambda,
		#	x = x,
		#	nx = length(x),
		#	err = rep(0, length(x)),
		#	y = rep(0, length(x)),
                #       logvalue = as.integer(logvalue),
		#	CLASSES = c("integer", "double", "integer", "double",
		#		"double", "integer"))
		#out <- tmp$y
            nx = length(x)
            err = rep(0, nx)
            y = rep(0,nx)
            tmp <- .C("besselM3",
			as.integer(lambda),
			as.double(x),
			as.integer(nx),
			as.double(err),
			y = as.double(y),
                        as.integer(logvalue),
                        PACKAGE="QRMlib");
                  #You must use 'y = as.double(y)' if you want to use tmp$y as a list element.
                  #The LHS (y =) of the notation 'y = as.double(y)' defines an 'output' and 
                  #the RHS defines an 'input'. 
                  #If you want to use solely 'as.double(y)', you have no named 'output' parameter so you
                  #must set 'out <- tmp[[5]]' since 5 is the order-number of the y input parameter.
                  #Parameter orders are 'lambda,x,nx,err,y,logvalue' = (1,2,3,4,5,6) so y is 
                  # 5th parameter in list. The function name "besselM3' is parameter 0.
                  #Hence you must refer to lambda as tmp[[0]], to x as tmp[[1]] since they are not 'named'.
                  #If you had used 'lambda = as.double(lambda)', you could use tmp$lambda
		out <- tmp$y
	}
	else if(intype == "half") {
		out <- sqrt(pi/2) * x^(-1/2) * exp( - x)
		if(logvalue)
			out <- log(out)
	}
	else if((intype == "fractional") | (intype == "halfinteger")) 
      {
            # SU 06/08/2006: made changes for consistency with R rather than S-Plus
            # The method call in C-language is 
            # void besselM3f(double *lambda, double *x, long *n, double *err, double *y, long *logvalue)
		#tmp <- .C("besselM3f",
		#	lambda = lambda,
		#	x = x,
		#	nx = length(x),
		#	err = rep(0, length(x)),
		#	y = rep(0, length(x)),
		#	logvalue = as.integer(logvalue),
		#	CLASSES = c("double", "double", "integer", "double",
		#		"double", "integer"))
		#out <- tmp$y
            nx = length(x)
            err = rep(0, nx)
            y = rep(0,nx)
            tmp <- .C("besselM3f",
			as.double(lambda),
			as.double(x),
			as.integer(nx),
			as.double(err),
			y = as.double(y),
                        as.integer(logvalue),
                        PACKAGE="QRMlib");
                  #You must use 'y = as.double(y)' if you want to use tmp$y as a list element.
                  #The LHS (y =) of the notation 'y = as.double(y)' defines an 'output' and 
                  #the RHS defines an 'input'. 
                  #If you want to use solely 'as.double(y)', you have no named 'output' parameter so you
                  #must set 'out <- tmp[[5]]' since 5 is the order-number of the y input parameter.
                  #Parameter orders are 'lambda,x,nx,err,y,logvalue' = (1,2,3,4,5,6) so y is 
                  # 5th parameter in list. The function name "besselM3f' is parameter 0.
                  #Hence you must refer to lambda as tmp[[0]], to x as tmp[[1]] since they are not 'named'.
                  #If you had used 'lambda = as.double(lambda)', you could use tmp$lambda.
		out <- tmp$y;

	}
	out
}

#############################################################################
#Simulate a multivariate generalized hyperbolic sample
rmghyp <- function(n,lambda,chi,psi,Sigma = equicorr(d, rho), mu = rep(0, d), gamma=rep(0,d),d = 2, rho = 0.7) 
{
     d <- dim(Sigma)[1.]
     W <- rGIG(n,lambda,chi,psi)
      m1 <- rmnorm(n, Sigma = Sigma)
      m2 <- matrix(rep(sqrt(W), d), ncol = d)
      offsetmatrix <- matrix(gamma, nrow = n, ncol = d, byrow = TRUE)
      mu.matrix <- matrix(mu, nrow = n, ncol = d, byrow = TRUE)
      return(m1 * m2 + offsetmatrix * m2^2 + mu.matrix)
}


##########################################################################################
#fit to a normalized hyperbolic.  See pp. 81-85 of QRM
fit.mNH <- function(data,symmetric=FALSE,case="NIG",kvalue=NA,nit=2000,tol=1e-10)
{
    if (is.matrix(data)==FALSE) data <- as.matrix(data);
    mix.pars <- switch(case,NIG=c(-0.5,1,1),hyp=c(1,1,1));
    optpars <- c(2,3);

    #optfunc() will be passed to MCECMupdate() below as an argument to be passed to optim().  
    #11/1/2007: added two additional parameters: "mixparam" will be mix.pars[1]and
    #"greeks" will be a list holding delta (vector), eta (vector), xi (scalar) respectively.
    #These parameters will be set (rather than globally assign()ed in MCECMupdate()prior to 
    #its call to optim() where it passes optfunc() and the extra parameters defined in this 
    #new version of optfunc(). CRITICALLY, the mixparam in fit.mNH() is mix.pars[1] whereas
    #it is mix.pars[3] in fit.mst()
    optfunc <- function(thepars, mixparams, greeks)
    {
      #parameters in order interpreted as lambda,chi,psi,delta,eta,xi:
      if(!is.list(greeks)) stop("greeks must be a list");
      if(!is.vector(mixparams)) stop("mixparams argument to optfunc() must be vector");
      MCECM.Qfunc(mixparams[1],thepars[1],thepars[2],greeks[[1]],greeks[[2]],greeks[[3]]);
    }

    #OLD METHOD NOW DEPRECATED BECAUSE IT IMPLICITLY RELIES ON ASSIGN()
    #optfunc <- function(thepars)
    #{
    #  MCECM.Qfunc(mix.pars.nl[1],thepars[1],thepars[2],greek.stats[[1]],greek.stats[[2]],greek.stats[[3]])
    #}

    n <- dim(data)[1];
    d <- dim(data)[2];
    Xbar <- apply(data,2,mean);
    mu <- Xbar;
    Sigma <- var(data);

    #Test kvalue passed in as parameter; it NA, set it via Sigma:
    if (is.na(kvalue)) kvalue <- det(Sigma);
    gamma <- rep(0,length(mu));
    scale.factor <- (det(Sigma)/kvalue)^(1/d);
    Sigma <- Sigma/scale.factor;
    i <- 0;
    #Use multivariate generalized hyperbolic distribution to sum the densities across data observations:
    ll <-  sum(dmghyp(data,mix.pars[1],mix.pars[2],mix.pars[3],mu,Sigma,gamma,logvalue=TRUE));
    closeness <- 100;
    while ((closeness > tol) & (i < nit))
    {
      i <- i+1;
      EMresult <- EMupdate(data,mix.pars,mu,Sigma,gamma,symmetric,scaling=TRUE,kvalue);
      mu <- EMresult$mu;
      Sigma <- EMresult$Sigma;
      gamma <- EMresult$gamma;

      MCECMresult <- MCECMupdate(data,mix.pars,mu,Sigma,gamma,optpars,optfunc,xieval=FALSE);
      mix.pars <- MCECMresult$mix.pars;
      conv <- MCECMresult$conv;
      conv.type <- MCECMresult$convtype;
      ll.old <- ll;

      ll <- sum(dmghyp(data,mix.pars[1],mix.pars[2],mix.pars[3],mu,Sigma,gamma,logvalue=TRUE));
      closeness <- abs((ll-ll.old)/ll.old);
    } #end while
    lambda <- mix.pars[1];
    chi <- mix.pars[2];
    psi <- mix.pars[3];
    mult <- det(Sigma)^(-1/d);
    Sigma <- Sigma*mult;
    chi <- chi/mult;
    psi <- psi*mult;
    gamma <- gamma*mult;
    mix.pars <- c(lambda,chi,psi);
    names(mix.pars) <- c("lambda","chi","psi");
    EW <- EGIG(lambda,chi,psi);
    EW2 <- EGIG(lambda,chi,psi,2);
    varW <- EW2-EW^2;
    beta <- as.vector(solve(Sigma) %*% gamma);
    mean <- as.numeric(mu+EW*gamma);
    covariance <- EW*Sigma + varW*outer(gamma,gamma);
    if (d >1)
      cor <- CovToCor(Sigma)
    else cor <- 1;
    delta <- as.numeric(sqrt(chi));
    alpha <- as.numeric(sqrt(psi+ t(beta) %*% Sigma %*% beta));
    alt.pars <- list(beta=beta,alpha=alpha,delta=delta);
    list(mix.pars=mix.pars, mu = mu, Sigma=Sigma, gamma=gamma, ll.max=ll,alt.pars=alt.pars,mean=mean,covariance=covariance,correlation=cor);
}
#######################################################################################

    MCECM.Qfunc <- function(lambda,chi,psi,delta,eta,xi){
      out <- NA
      if ((chi>0) & (psi>0)){
        n <- length(delta)
        term1 <- (lambda-1)*sum(xi)
        term2 <- -chi*sum(delta)/2
        term3 <- -psi*sum(eta)/2
        term4 <- -n*lambda*log(chi)/2+ n*lambda*log(psi)/2 -n*besselM3(lambda,sqrt(chi*psi),logvalue=TRUE)
        out <- -(term1+term2+term3+term4)
      }
      if ((psi==0) & (lambda <0)){
        n <- length(delta)
        nu <- chi
        term1 <- -n*nu*log(nu/2)/2
        term2 <- nu*(sum(xi)+sum(delta))/2
        term3 <- n*log(gamma(nu/2))
        out <- term1+term2+term3
      }
      if((chi==0) * (lambda>0)) stop("vg not implemented")
      out

    }

    #############################################################################
#Note that kvalue is used only if scaling=TRUE so you can pass the default (i.e.
#omit the parameter if you are NOT scaling.
  EMupdate <- function(data,mix.pars,mu,Sigma,gamma,symmetric,scaling=TRUE,kvalue=1)
  {
    d <- dim(data)[2]
    n <- dim(data)[1]
    lambda <- mix.pars[1]
    chi <- mix.pars[2]
    psi <- mix.pars[3]
    Q <- mahalanobis(data,mu,Sigma)
    Offset <- t(gamma) %*% solve(Sigma) %*% gamma
    delta <- EGIG(d/2-lambda,psi+Offset,Q+chi)
    delta.bar <- mean(delta)
    eta <- EGIG(lambda-d/2,Q+chi,psi+Offset)
    eta.bar <- mean(eta)
    delta.matrix <- matrix(delta,nrow=n,ncol=d,byrow=FALSE)
    if (symmetric==TRUE)
      gamma <- rep(0,d)
    else
      {
        Xbar <- apply(data,2,mean)
        Xbar.matrix <- matrix(Xbar,nrow=n,ncol=d,byrow=TRUE)
        Xbar.matrix <- Xbar.matrix-data
        gamma <- apply(delta.matrix*Xbar.matrix,2,sum)/(n*delta.bar*eta.bar-n)
      }
    mu <- (apply(delta.matrix*data,2,sum)/n - gamma)/delta.bar
    mu.matrix <- matrix(mu,nrow=n,ncol=d,byrow=TRUE)
    standardised <- data-mu.matrix
    tmp <- delta.matrix*standardised
    Sigma <- (t(tmp) %*% standardised)/n - outer(gamma,gamma)*eta.bar
    if (scaling){
      scale.factor <- (det(Sigma)/kvalue)^(1/d)
      Sigma <- Sigma/scale.factor
    }
    list(mu=mu,Sigma=Sigma,gamma=gamma)
  }

###################################################################
MCECMupdate <- function(data,mix.pars,mu,Sigma,gamma,optpars,optfunc,xieval=FALSE)
{
  d <- dim(data)[2];
  n <- dim(data)[1];
  lambda <- mix.pars[1];
  chi <- mix.pars[2];
  psi <- mix.pars[3];
  Q <- mahalanobis(data,mu,Sigma);
  Offset <- t(gamma) %*% solve(Sigma) %*% gamma;
  delta <- EGIG(d/2-lambda,psi+Offset,Q+chi);
  eta <- EGIG(lambda-d/2,Q+chi,psi+Offset);
  xi <- 0;
  if (xieval) xi <- ElogGIG(lambda-d/2,Q+chi,psi+Offset);
  #record the value of the initial 2nd parameter input into another variable
  thepars <- mix.pars[optpars];

  #11/1/2007: replace assign() statements by building parameters to pass to optim().  These represent
  #the extra parameters to be sent to optfunc() via optim():
  #We will send mix.pars vector as mixparams since fit.mst() and fit.mNH() use different
  #values from the vector.  E.g. mixparam.nl <- mix.pars[3] would be used by fit.mst() but
  #mix.pars[1] would be used by fit.mNH().
  #mixparam.nl < mix.pars[3];  
  #Build the "greeks" parameter:
  greek.stats <- list(delta=delta, eta=eta, xi=xi); 

  #SU: 06/07/2006. Trying to replace all S-Plus nlmin(f,x) calls by calls to R-function nlm(f,p)fails.
  # Use optim() instead.  The following are the commented-out S-Plus lines.
  #tmp <- nlmin(optfunc,thepars)
  #mix.pars[optpars] <- tmp$x
  #list(mix.pars=mix.pars,conv=tmp$converged,convtype=tmp$conv.type)

  #These are the replacement R-language lines (optfunc is an input function name):
  #We have ADDED parameters "mixparam" and "greeks" to optfunc() along with "thepars"
  #which are the original parameters to optimize over.  Hence we must add these new
  #parameters to the list passed to optim():
  #*******************************
  optimout <- optim(par=thepars,fn=optfunc,gr=NULL,mixparams=mix.pars, greeks=greek.stats, method="BFGS");
  #****************************************

  mix.pars[optpars] <- optimout$par
  if(optimout$convergence == 0)
  {
    conv = TRUE
    convtype = "Convergence succeeded via BFGS quasi-Newton method"
  }
  else
  {
    conv = FALSE 
    convtype = "Convergence failed via BFGS quasi-Newton method"
  }
  #SU: conv.type is a character string return from nlmin() in S-Plus.  It is replaced by
  # conv in R-language
  list(mix.pars=mix.pars,conv,convtype)
}


###################################
EGIG <- function(lambda,chi,psi,k=1){
  if ((chi[1]>0) & (psi[1]>0)){
    term1 <- k*log(chi/psi)/2
    term2 <- besselM3(lambda+k,sqrt(chi*psi),logvalue=TRUE)
    term3 <- besselM3(lambda,sqrt(chi*psi),logvalue=TRUE)
    out <- exp(term1+term2-term3)
  }
  else if ((chi[1]==0) & (lambda>0)){
    alpha <- lambda
    beta <- psi/2
    out <- gamma(k+alpha)/(gamma(alpha)*beta^k)
  }
  else if ((psi[1]==0) & (lambda <0)){
    alpha <- -lambda
    beta <- chi/2
    out <- (gamma(alpha-k)*beta^k)/gamma(alpha)
  }
  else stop("These GIG parameters are not allowed")
  out
}

#####################################################################

ElogGIG <- function(lambda,chi,psi){
 if ((chi[1]==0) & (lambda>0)){
    alpha <- lambda
    beta <- psi/2
    out <- psifunc(alpha)-log(beta)
  }
  else if ((psi[1]==0) & (lambda <0)){
    alpha <- -lambda
    beta <- chi/2
    out <- log(beta)-psifunc(alpha)
  }
 else {
stop("Log Moment of general GIG not implemented")
 }
  out
}

#########################################################################

psifunc <- function(x = 2,logvalue=FALSE)
{
      #The original S-Plus code must be replaced:
	#tmp <- .C("psifunc",
	#	x = x,
	#	nx = length(x),
	#	err = rep(0, length(x)),
	#	y = rep(0, length(x)),
	#	CLASSES = c("double", "integer", "double",
	#			"double"))
	#out <- tmp$y
      nx = length(x)
      err = rep(0,nx)
      y = rep(0,nx)
	tmp <- .C("psifunc",
		as.double(x),
		as.integer(nx),
		as.double(err),
		y = as.double(y),
                PACKAGE="QRMlib");
            #You must use 'y = as.double(y)' if you want to use tmp$y as a list element.
            #The LHS (y =) of the notation 'y = as.double(y)' defines an 'output' and 
            #the RHS defines an 'input'. 
            #If you want to use solely 'as.double(y)', you have no named 'output' parameter so you
            #must set 'out <- tmp[[4]]' since 4 is the order-number of the y input parameter.
            #Parameter orders are 'x,nx,err,y' = (1,2,3,4) so y is 4th parameter in list.
	out <- tmp$y;
	if(logvalue)
          out <- log(out)
        out
}


#########################################################################
#Fit data to multivariate student t distribution.  See descriptions pp. 81-85 QRM
#There is no class data.t.5d data set so the first parameter should NOT
#be set to default as data=data.t.5d:
#nit is number of iterations to try
fit.mst <- function(data,nit=2000,tol=1e-10)
{
    if (is.matrix(data)==FALSE) data <- as.matrix(data);
    mix.pars <- c(-4,8,0);
    optpars <- c(2);

    #optfunc() will be passed to MCECMupdate() below as an argument to be passed to optim().  
    #11/1/2007: added two additional parameters: "mixparams" will be vector mix.pars[]and
    #"greeks" will be a list holding delta (vector), eta (vector), xi (scalar) respectively.
    #These parameters will be set (rather than globally assign()ed in MCECMupdate()prior to 
    #its call to optim() where it passes optfunc() and the extra parameters defined in this 
    #new version of optfunc()
    #optfunc <- function(thepars, mixparam, greeks)
    optfunc <- function(thepars, mixparams, greeks)
    {
      #fit.mst() uses mixparams[3] whereas fit.mNH() uses mixparams[1]:
      #parameters in order to Qfunc() interpreted as lambda,chi,psi,delta,eta,xi:
      if(!is.list(greeks)) stop("greeks must be a list");
      if(!is.vector(mixparams)) stop("mixparams argument to optfunc() must be vector");
      MCECM.Qfunc(-thepars/2,thepars,mixparams[3],greeks[[1]],greeks[[2]],greeks[[3]]);

    }
    #Old Version now supplanted: 
    #optfunc <- function(thepar)
    #{
    #  MCECM.Qfunc(-thepar/2,thepar,mix.pars.nl[3],greek.stats[[1]],greek.stats[[2]],greek.stats[[3]])
    #}
    
    n <- dim(data)[1];
    d <- dim(data)[2];
    Xbar <- apply(data,2,mean);
    mu <- Xbar;
    Sigma <- var(data);
    gamma <- rep(0,d);
    i <- 0;
    #Use multivariate t distribution to sum the densities across data observations: 
    ll <-  sum(dmt(data,mix.pars[2],mu,Sigma,logvalue=TRUE));
    closeness <- 100;
    while ((closeness > tol) & (i < nit))
    {
      i <- i+1;
      #11/6/2007 SU: since you are not scaling in EMupdate, the kvalue is irrelevant and may be ignored,
      #i.e. don't pass anything and let the function set it to a default of 1:
      EMresult <- EMupdate(data,mix.pars,mu,Sigma,gamma,symmetric=TRUE,scaling=FALSE);  #,kvalue);
      mu <- EMresult$mu;
      Sigma <- EMresult$Sigma;
      gamma <- EMresult$gamma;
      #***************************************      
      MCECMresult <- MCECMupdate(data,mix.pars,mu,Sigma,gamma,optpars,optfunc,xieval=TRUE);
      #******************************************

      mix.pars <- MCECMresult$mix.pars;
      mix.pars[1] <- -mix.pars[2]/2;
      conv <- MCECMresult$conv;
      conv.type <- MCECMresult$convtype;
      ll.old <- ll;
      ll <-  sum(dmt(data,mix.pars[2],mu,Sigma,logvalue=TRUE));
      closeness <- abs((ll-ll.old)/ll.old);
    }  #end while
    lambda <- mix.pars[1];
    chi <- mix.pars[2];
    psi <- mix.pars[3];
    EW <- EGIG(lambda,chi,psi);
    EW2 <- EGIG(lambda,chi,psi,2);
    varW <- EW2-EW^2;
    Sigma <- symmetrize(Sigma);
    beta <- as.vector(solve(Sigma) %*% gamma);
    mean <- as.numeric(mu+EW*gamma);
    covariance <- EW*Sigma + varW*outer(gamma,gamma);
    if (d >1)
      cor <- CovToCor(Sigma)
    else cor <- 1;
    nu <- chi;
    list(nu=nu, mu = mu, Sigma=Sigma, gamma=gamma,ll.max=ll,mean=mean,covariance=covariance,correlation=cor);
}

#######################################################################################
