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

# Contact: Alexander J. McNeil:  a.j.mcneil@hw.ac.uk */
# Contact for R-language version: Scott Ulman  scottulman@hotmail.com
# Note in the R-translations that TRUE has been substituted throughout the 
# code for T (and FALSE for F) when setting default parameter values as 
# specified in section 2.4 of an "Introduction to R" (R-intro.pdf).
# Otherwise "R CMD check -QRMlib" returned the following type of error when
# trying to run an example:
# > BiDensPlot(func=dmnorm,mu=c(0,0),Sigma=equicorr(2,-0.7))
#   Error in func(data, ...) : F used instead of FALSE
#   Execution halted
#########################################################
fit.st <- function(data)
{
  #10/4/2007: Fixed for R-2.6.0. In R-2.5.0 and prior, as.numeric() could process timeSeries objects directly.  
  #Beginning with R-2.6.0 you must test for timeSeries and pass only the @Data slot to as.numeric():
  if (exists("is.R") && is.function(is.R) && is.R()) 
  {
    if(is.timeSeries(data)) data <- data@Data; 
  } 

  mu <- mean(data); #sample mean
  m2 <- mean((data-mu)^2); #autocovariance function
  m4 <- mean((data-mu)^4);
  nu <- 4 + (6*m2^2)/(m4-3*m2^2);
  sigma <- sqrt((nu-2)*m2/nu);
  theta <- c(nu,mu,sigma);
  
  #10/5/2007: added 2nd parameter to negloglik() to avoid using assign()
  negloglik <- function(theta, datavector)
  {   
     #dt() is the density function for the student-t distribution with the first parameter
     # being a vector of quantiles (q) and the second parameter df being the degrees of freedom.
     #A value of the density is returned for each value in the data set and these values are summed.
      - sum(log(dt((datavector - theta[2])/abs(theta[3]), df = abs(theta[1]))) - log(abs(theta[3])));
  }
  #10/5/2007: added 3rd parameter to be passed to negloglik() to avoid using assign()
  optimfit <- optim(theta, fn=negloglik, gr=NULL, datavector = data, method="BFGS");
  par.ests <- optimfit$par;
  if(optimfit$convergence == 0) #if code is 0 we are OK; 
    converged <- TRUE
  else
    converged <- FALSE;

  # replace parameter estimates with their absolute values:
  par.ests[1] <- abs(par.ests[1]);
  par.ests[2] <- abs(par.ests[2]);
  #10/5/2007: added 3rd parameter to hessb() to be passed to negloglik() to avoid using assign()
  nItheta <- hessb(negloglik,par.ests, datavector=data);
  
  asymp.cov <- solve(nItheta);
  #10/5/2007: added 2nd parameter to be passed to negloglik() to avoid using assign()

  loglh.max <- -negloglik(par.ests, datavector=data);
  par.ses <- sqrt(diag(asymp.cov));
  names(par.ests) <- c("nu","mu","sigma");
  names(par.ses) <- names(par.ests);
  dimnames(asymp.cov) <- list(names(par.ests),names(par.ests));

  #Return a list telling whether convergence occurred, the parameter estimates, the std error of
  #parameter estimates, the asymptotic covariance matrix, and the value of the maximize loglikelihood
  list(converged=converged,par.ests=par.ests, par.ses=par.ses, asymp.cov=asymp.cov, ll.max=loglh.max);      
}

#####################################################################################
#Create a mean-variance normal mixture of random variables called the generalized hyperbolic
#It is not necessarily elliptical but its univariate version will be.  See p. 78 in QRM.
# This is the GH model.
rghyp <- function(n, lambda, chi, psi, mu=0, gamma=0)
{
      #generate a series of random Generalized Inverse Gaussian variables
      # see p. 77 of QRM text
	W <- rGIG(n, lambda, chi, psi);
      # Generate a similar random sequence of standard normals: 
	Z <- rnorm(n);
      #Mix the two distributions using equation 3.25 (p. 77) but with gamma possibly 0 or a scalar
	sqrt(W) * Z + mu + gamma * W;
}

#####################################################################################
#Generate a random series for Generalized Hyperbolic with the so-called 
# alpha-delta-capDelta-beta parameterization. See p.80 in QRM book.
rghypB <- function(n, lambda, delta, alpha, beta=0, mu=0)
{
  rghyp(n,lambda,delta^2,alpha^2-beta^2,mu,beta);
}

#################################################################################
#generate a density for the generalized hyperbolic function
dghyp <- function(x, lambda, chi, psi, mu=0, gamma=0, logvalue=FALSE)
{
	Q <- (x - mu)^2;
	infunc <- sqrt((chi + Q) * (psi + gamma^2));
	top <- besselM3((lambda - 1/2), infunc, logvalue = TRUE);
	bottom <- (1/2 - lambda) * log(infunc);
	tilt <- gamma * (x - mu);
	const.top <- (-lambda/2) * log(psi * chi) + lambda*log(psi) + (1/2-lambda)*log(psi + gamma^2);
	const.bottom <- log(2 * pi)/2 + besselM3(lambda, sqrt(chi * psi), logvalue = TRUE);
	out <- (const.top + top + tilt) - (const.bottom + bottom);
        if (!(logvalue)) out <- exp(out);
	out;
}

############################################################################
#Get the density for the so-called alpha-delta-capDelta-beta parameterization. 
#See p.80 in QRM book.
dghypB <- function(x, lambda, delta, alpha, beta=0, mu=0, logvalue=FALSE)
{
  dghyp(x,lambda,delta^2,alpha^2-beta^2,mu,beta,logvalue);
}

####################################################################
#Fit data to univariate generalized hyperbolic.  See p. 80 in QRM.  
#The default case is NIG requiring lambda = -1/2.  The NIG has a slightly heavier
# tail than the generalized hyperbolic set by case="hyp"  
fit.NH <- function(data, case="NIG", symmetric = FALSE, se=FALSE)
{
  #10/4/2007: Fixed for R-2.6.0. In R-2.5.0 and prior, as.numeric() could process timeSeries objects directly.  
  #Beginning with R-2.6.0 you must test for timeSeries and pass only the @Data slot to as.numeric():
  if (exists("is.R") && is.function(is.R) && is.R()) 
  {
    if(is.timeSeries(data)) data <- data@Data; 
  } 

  mu <- mean(data);
  
  #The R-function var() does not have an unbiased=F alternative. 
  #In S-Plus, unbiased=F implies the sum will be divided by N (the series length) rather than 
  # by N-1 which is always used by R. The value will be biased and hence smaller.
  #Hence we need to multiply the R-var() by (N-1)/N to get the unbiased=F (biased)S-Plus-var()to match
  # McNeil's code.  Hence we comment out the following S-Plus code from McNeil.
  #variance <- var(data,unbiased=F) # S-Plus call
  #and replace it with the following longer version:
  N = length(data);
  if(N==1) stop ("only one observation in data sent to fit.nH");
  variance <- (N-1)*var(data)/N; # creates the 'biased' (smaller) version which is divided by N
  #Use the S-Plus version of kurtosis (properly modified) since there is a discrepancy between
  # the R and S-Plus versions (the R "excess" version should match the S-Plus "moment" method but
  # it doesn't and we want to match McNeil's S-Plus results here.
  #Comment out the R-version:
  #ekurt <- kurtosis(data,method="moment")  
  #Substitute our analogue version built from S-Plus type code:
  ekurt <- kurtosisSPlus(data,method="moment");

  #Added test for parameter from case.  In S-Plus, the return from a switch() is 
  #the value of the expression that is picked, or NULL if no expression is picked. 
  lambda <- switch(case,NIG=-1/2,hyp=1); #returns -1/2 if 'case' is "NIG" and 1 if 'case' is "hyp"
  if(is.null(lambda)) stop("case must be 'NIG' or 'hyp' in call to fit.NH");
  
  alpha <- sqrt(3/(variance*ekurt));
  delta <- alpha*variance;
  chi <- delta^2;
  psi <- alpha^2;
  gamma <- 0;
  #set local theta to have either 3 or 4 parameters
  if (symmetric == TRUE)
  {
    theta <- c(chi,psi,mu)
  }
  else
    theta <- c(chi,psi,mu,gamma);
    
   
  #Negloglik() will be called implicitly from optim() where it is passed as the 2nd parameter.
  #It will also be called implicitly from hessb().
  #10/5/2007: added multiple parameters to negloglik to avoid using assign()
  #These new parameters input to negloglik() must also be passed to optim() and hessb() below
  negloglik <- function(theta, datavector, symmIn, lambdaIn)
  {
      #gammaLocal is a local variable set to either 0 or theta[4] 
      if(symmIn) #effectively set final parameter to 0
        gammaLocal <- 0
      else gammaLocal <- theta[4];
 
      negloglikSum <-
      #Minus sum() will be returned by negloglik() by default since this is the last step in function     
      - sum(dghyp(x= datavector, lambda = lambdaIn, chi=abs(theta[1]),
                  psi = abs(theta[2]), mu = theta[3], gamma
                  = gammaLocal, logvalue=TRUE))
       
      return(negloglikSum)
  }  #end negloglik()

  #10/5/2007: named parameters 1-3 and passed additional paramters 4-6 to be passed to negloglik()
  optimfit <- optim(par = theta, fn=negloglik, gr=NULL, datavector= data, symmIn = symmetric, lambdaIn= lambda, method="BFGS");
  par.ests <- optimfit$par; 
  if(optimfit$convergence == 0) #if code is 0 we are OK; 
    converged <- TRUE
  else
    converged <- FALSE;
  #replace first two parameter estimators by their absolute values:
  par.ests[1] <- abs(par.ests[1]);
  par.ests[2] <- abs(par.ests[2]);
  if(symmetric)
    names(par.ests) <- c("chi","psi", "mu")
  else names(par.ests) <- c("chi","psi", "mu", "gamma");

  #Set the loglikelihood maximum value:
  ll.max <-  - negloglik(par.ests, datavector= data, symmIn = symmetric, lambdaIn= lambda);
  if(se) 
  {
    #10/5/2007 Must also pass additional parameters passed to new negloglik() to hessb() which 
    #calculates the hessian matrix of negloglik() function:
    hessmatrix <- hessb(negloglik,par.ests,datavector= data, symmIn = symmetric, lambdaIn= lambda);
    vcmatrix <- solve(hessmatrix);
    par.ses <- sqrt(diag(vcmatrix));
    names(par.ses) <- names(par.ests);
    dimnames(vcmatrix) <- list(names(par.ests), names(par.ests));
  }
  else 
  {
    par.ses <- NA;
    vcmatrix <- NA;
  }
  chi <- par.ests[1];
  psi <- par.ests[2];
  mu <- par.ests[3];
  if(!(symmetric))
    gamma <- par.ests[4];
  alt.pars <- c(sqrt(chi), sqrt(psi + gamma^2),gamma,mu);
  names(alt.pars) <- c("delta", "alpha","beta","mu");
  
  list(converged = converged, case=case,symmetric=symmetric,par.ests = par.ests, par.ses = par.ses,
       alt.pars = alt.pars, vcmatrix=vcmatrix, ll.max= ll.max);
}

#################################################################################
#SU: Added function to resolve conflict between R and S-Plus.
#This method parallels the calculation of kurtosis in S-Plus.  The values calculated
#by R and S-Plus differ when we use the call kurtosis(x, method="moment") which causes
#serious consequences in the fit.NH() function call.  Hence we introduce the S-Plus
#version used by McNeil here to call from fit.NH().
#Note the "moment" method in S-Plus subtracts three from (sum(x^4)/n)/(sum(x^2)/n)^2.
# In R, the "excess" method subtracts three and hence SHOULD be equivalent to the "moment"
# method of S-Plus.  However, I found when using a vector of 3998 values that the value
# calculated in R was
# kurtosis(data.NIG.splus,method="excess")
# [1] 1.347410
#whereas the value reported in S-Plus was
# kurtosis(data.NIG.splus, method= "moment")
#kurtosis(data.NIG,method="moment")
# [1] 1.349585
#Hence I decided to introduce the following S-Plus code modified to R:
#where "moment" approximates the "excess" method of R
kurtosisSPlus <- function(x, na.rm = FALSE, method = "fisher")
{
   method <- char.expand(method, c("fisher", "moment"), stop(
		"argument 'method' must match either \"fisher\" or \"moment\""
		))
    #Replace these S-Plus lines of code
	#if(na.rm) {
	#	wnas <- which.na(x)
	#	if(length(wnas))
	#		x <- x[ - wnas]
	#}
	#else if(anyMissing(x))
	#	return(NA)
    #with these R-language lines of code
    #if we are told to remove NA items, then use the series with omitted values
    if(na.rm) 
    {
       x <- na.omit(x)
    }
    else if(length(x) != length(na.omit(x)))
       return(NA)
    
	n <- as.double(length(x))
	if(method == "fisher" && n < 4)
		return(NA)
	x <- x - mean(x)
	if(method == "moment")
		(sum(x^4)/n)/(sum(x^2)/n)^2 - 3
	else ((n + 1) * (n - 1) * ((sum(x^4)/n)/(sum(x^2)/n)^2 -
			(3 * (n - 1))/(n + 1)))/((n - 2) * (n -3))
}
##########################################################
