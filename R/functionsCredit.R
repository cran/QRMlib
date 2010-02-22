# QRMlib: version 1.4.5
# this file is a component of QRMlib 

# Copyright (C) 2005-07 Alexander McNeil
# R-language additions Copyright (C) 2006-2010 Scott Ulman

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
# R-language contact:  Scott Ulman: scottulman@hotmail.com
# Note in the R-translations that TRUE has been substituted throughout the 
# code for T (and FALSE for F) when setting default parameter values as 
# specified in section 2.4 of an "Introduction to R" (R-intro.pdf).
# Otherwise "R CMD check -QRMlib" returned the following type of error when
# trying to run an example:
# > BiDensPlot(func=dmnorm,mu=c(0,0),Sigma=equicorr(2,-0.7))
#   Error in func(data, ...) : F used instead of FALSE
#   Execution halted
#########################################################
momest <- function(data, trials, limit = 10.)
{
	out <- rep(NA, limit)
	for(k in 1.:limit) {
		term <- 1.
		for(j in 1.:k)
			term <- (term * (data - (j - 1.)))/(trials - (j - 1.))
		out[k] <- mean(term)
	}
	out
}

#######################################################

cal.probitnorm <-function(pi1 = 0.1837, pi2 = 0.0413)
{
       #In R, the library which contains the function pmvnorm() is mvtnorm which is a dependency
       #of QRMlib. Hence it is automatically loaded and we don't need to reload it.
       #library(mvtnorm);
	rootfunc <- function(x, k1, k2)
	{
             #Major difference between S-Plus and R in this call:
		#k2 - pmvnorm(c(qnorm(k1), qnorm(k1)), rho = x); #S-Plus way
             #R-language way: must set correlation matrix and provide lower bound for each of the 
             #pair of bivariate variables:
             k2 - pmvnorm(lower=-Inf, upper=c(qnorm(k1), qnorm(k1)), corr = equicorr(2,x));
	}
	out <- uniroot(rootfunc, c(0.0001, 0.999), k1 = pi1, k2 = pi2)
	rho <- out$root
       #In R, we no longer need to unload the mvtnorm package since it is a dependency of QRMlib:
       #detach("package:mvtnorm");

	mu <- qnorm(pi1)/sqrt(1 - rho)
	sigma <- sqrt(rho)/sqrt(1 - rho)
	c(mu = mu, sigma = sigma, rho.asset = rho)
}

###############################################

 cal.beta <- function(pi1 = 0.1837, pi2 = 0.0413)
{
	a <- (pi1 * (pi1 - pi2))/(pi2 - pi1^2.)	
	b <- ((1. - pi1) * (pi1 - pi2))/(pi2 - pi1^2.)
	c(a = a, b = b)
}

#####################################################

 cal.claytonmix <- function(pi1 = 0.1837, pi2 = 0.0413)
{
	rootfunc <- function(theta, k1, k2)
	{
		log(2 * k1^( - theta) - 1) + theta * log(k2)
	}
	out <- uniroot(rootfunc, c(1e-010, 1), k1 = pi1, k2 = pi2)
	theta <- out$root
	c(pi = pi1, theta = theta)
}

######################################################

 pprobitnorm <- function(q, mu, sigma)
{
	pnorm((qnorm(q) - mu)/sigma)
}

####################################################
 pclaytonmix <- function(q, pi, theta)
{
	1 - pgamma((( - log(q))/(pi^( - theta) - 1)), (1/theta))
}

########################################################
#switched notation so x passed as 1st argument rather than q which goes with pprobitnorm()
dprobitnorm <- function(x, mu, sigma)
{
	(dnorm((qnorm(x) - mu)/sigma))/(dnorm(qnorm(x)) * sigma)
}

########################################################
#switched notation so x passed as 1st argument rather than q which goes with pclaytonmix()
dclaytonmix <- function(x, pi, theta)
{
	dgamma((( - log(x))/(pi^( - theta) - 1)), (1/theta))/(x * (pi^( - theta) - 1))
}

########################################################

rprobitnorm <- function(n, mu, sigma)
{
	pnorm(rnorm(n, mu, sigma))
}

########################################################

rlogitnorm <- function(n, mu, sigma)
{
	(1. + exp( - rnorm(n, mu, sigma)))^(-1.)
}

#######################################################

rclaytonmix <- function(n, pi, theta)
{
	exp( - (pi^( - theta) - 1) * rgamma(n, 1/theta))
}

#############################################################

rtcopulamix <- function(n, pi, rho.asset, nu)
{
	W <- nu/rchisq(n, nu)
	THETA <- rnorm(n)
	pnorm((qt(pi, nu)/sqrt(W) - THETA * sqrt(rho.asset))/sqrt(1 - rho.asset))
}

##############################################################

rbinomial.mixture <- function(n = 1000, m = 100, model = "probitnorm", ...)
{
	mixdist <- eval(parse(text = paste("r", model, sep = "")))
	Q <- mixdist(n, ...)
	rbinom(n, m, Q)
}

####################################################################

# Fits a binomial distribution to data for a single group
 fit.binomial <- function(M, m)
{
	phat <- sum(M)/sum(m)
	loglik <- function(p, M, m)
	{
		logb(p) * sum(M) + logb(1. - p) * (sum(m) - sum(M))
	}
	llmax <- loglik(phat, M, m)
	pse <- sqrt((phat * (1 - phat))/length(M))
	pi <- phat
	pi2 <- phat^2
	rhoY <- 0
	list(par.ests = phat, par.ses = pse, maxloglik = llmax, pi = pi, pi2 = pi2, rhoY = rhoY)
}

#########################################################
# Fits a binomial-beta distribution to data for a single group; e.g.
#fits default and obligors to binomial distribution
fit.binomialBeta <- function(M, m, startvals = c(2, 2), ses = FALSE)
{
#ARGUMENTS:
#  M        vector of number of successes (defaults)
#  m        vector of number of trials (obligors)-vector 
#           M and M will have equal length which will probably be the number of 
#           credit classes/ratings from a Ratings Agency
# startvals starting point for parameter estimations
# ses       Indicates whether standard errors should be calculated

   #10/5/2007: added parameters 2-3 to function to replace assign()
   negloglik <- function(theta, defaults, trials)
   {
     length(trials) * lbeta(theta[1]^2, theta[2]^2) - sum(lbeta(theta[1]^2 + defaults, 
          theta[2]^2 + trials - defaults));
   }
   #10/5/2007: added parameters 2-3 to truenegloglik() to avoid assign(); 
   truenegloglik <- function(theta, defaults, trials)
   {
      length(trials) * lbeta(theta[1], theta[2]) - sum(lbeta(theta[1] + 
          defaults, theta[2] + trials - defaults));
   }
   #10/5/2007: removed assign() settings
 
   theta <- sqrt(startvals);
  
   #10/5/2007: added parameters 2-3 to function to replace assign()
   optimout <- optim(theta,negloglik, defaults=M, trials=m, method="BFGS");
   par.ests <- optimout$par;
   if(optimout$convergence == 0)
   {
      conv = TRUE;
      convtype = "Convergence succeeded via BFGS quasi-Newton method";
   }
   else
   {
      conv = FALSE; 
      convtype = "Convergence failed via BFGS quasi-Newton method";
   }

   #10/5/2007: added parameters 2 and 3 which must be passed to negloglik()    
   loglik <-  - negloglik(par.ests, defaults=M, trials=m);
   par.ests <- par.ests^2;
   
	
   if(ses) 
   {
      #10/5/2007: added new parameters 3-4 to hessb; they must both be passed to truenegloglik():
      hessmatrix <- hessb(truenegloglik, par.ests, defaults=M, trials=m);
      varcov <- solve(hessmatrix);
      par.ses <- sqrt(diag(varcov));
   }
	
   else par.ses <- rep(NA, 2);

   pi <- par.ests[1]/sum(par.ests);
   pi2 <- (par.ests[1] + 1)/(sum(par.ests) + 1) * pi;
   rhoY <- (pi2 - pi^2)/(pi - pi^2);
       
  list(par.ests = par.ests, par.ses = par.ses, maxloglik = loglik, converged = conv, 
       details = convtype, pi = pi, pi2 = pi2, rhoY = rhoY);
}


####################################################################
#fit.binomialProbitnorm() differs from fit.binomialLogitnorm() only in the link function.  Here the 
#link funtion is pnorm() where in the logit version, it is 1/(1-exp(-z))).
#The argument M is a vector like 'Bdefaults' (the number of defaults) and m is a vector like 
#'Bobligors' (the number of obligors).  This function calls optim(...method="L-BFGS-B") so you may need to
#reset lower and upper parameter estimators if convergence occurs at an endpoint of either limit and run
#another test.
fit.binomialProbitnorm <- function(M, m, startvals = c(-1, 0.5),lowerParamLimits = c(-3.0, 0.02),
  upperParamLimits = c(1,0.9))
{
   #ARGUMENTS:
   #  M        vector of number of successes (defaults)
   #  m        vector of number of trials (obligors)-vector 
   #probitnormal is the one-factor KMV/CreditMetrics model. Link function for probitnorm is the normal.
   #'link' is the name of a function to be passed via 'link.fn' to the 'integrate()' method in 'negloglik()' after 
   #'link' is assigned() to 'link.fn'.  pnorm() is the probability distribution for a normal random variable.


   link <- pnorm;

   #This local function name will be passed to negloglik() as the 3rd parameter. The integrand will be
   #calculated as a function of p and has five other parameters (kk,n ,mu, sigma, lfunc where lfunc is
   #a link function passed in.
   integrand <- function(p, kk, n, mu, sigma, lfunc)
   {
      exp(kk * log(lfunc(mu + sigma * qnorm(p))) + (n - kk) * log(1 - lfunc(mu + sigma * qnorm(p))));
   }
   
   #Renamed this to integrand2; it is called directly from the code of this function and not indirectly 
   #through an internal function.  It is a function of p and has only four parameters. lfunc is name 
   #of link function passed in as one of the parameters. 
   #It assumes probability of exactly kk obligors defaulting:
   integrand2 <- function(p, mu, sigma, kk, lfunc)
   {
      exp(kk * log(lfunc(mu + sigma * qnorm(p))));
   }

   #10/9/2007: removed assign() statements.
   #Combine the defaults M in column 1 and the number of obligors m in column 2 of a matrix:
   MLdata <- cbind(M,m); 

   #10/9/2007: added parameters which replace assign(). We pass integrand which integrates over p
   #with 5 parameters including a link function:
   negloglik <- function(theta, defaultData, FNIntegrate=integrand, FNlink=link)
   {
      #10/9/2007: added tests to make sure some parameters are functions:
      if(!is.function(FNIntegrate)) stop("parameter FNIntegrate must be a function");
      if(!is.function(FNlink)) stop("parameter FNlink must be a function");

      n <- dim(defaultData)[1]; #determine the number of rows in defaultData
      #numeric() creates a vector of length n with all values equal to 0.  The sum of the elements of this vector 
      #will be the return value of negloglik()
      out <- numeric(n); 
      for(i in 1:n) 
      {
         #R version of integrate() returns list with $value rather than the S-Plus $integral.  
         #5-parameter version integrating over p: 
         out[i] <-  - logb(integrate(FNIntegrate, lower = 0., upper = 1.,  kk = defaultData[i, 1], 
                      n = defaultData[i, 2], mu = theta[1], sigma = theta[2]^2, lfunc = FNlink)$value)
      }
      sum(out);
   }
 
   theta <- c(startvals[1], sqrt(startvals[2]))
   #SU: Replace S-Plus nlmin(f,x) calls by calls to R-function optim().
   #out <- nlmin(negloglik, theta, max.fcal = 1000., max.iter = 1000.) #S-Plus
   #par.ests <- out$x
   #These are the replacement R-language lines ('negloglik' is an input function name):
   #Set lower and upper bounds on parameters: startvals are c(-1, 0.5). Hence use L-BFGS-B method rather than BFGS method.
   #10/9/2007: added parameters to optim() which need to be passed to negloglik:
   optimout <- optim(theta,negloglik,defaultData=MLdata, FNIntegrate=integrand, FNlink=link,
           method="L-BFGS-B", lower=lowerParamLimits ,upper=upperParamLimits)
   par.ests <- optimout$par

   if(optimout$convergence == 0)
   {
      conv = TRUE
      convtype = "Convergence succeeded via L-BFGS-B quasi-Newton method"
   }
   else
   {
      conv = FALSE 
      convtype = "Convergence failed via L-BFGS-B quasi-Newton method"
   }
   #10/9/2007: added further parameters to negloglik:
   loglik <-  - negloglik(par.ests, defaultData=MLdata, FNIntegrate=integrand, FNlink=link);
   par.ests[2] <- par.ests[2]^2;

   #integrate() performs NUMERICAL INTEGRATION. See pp. 354-5 in QRM book for a description (Examples 8.12
   #and 8.13. Integrate to get first moment (probability of precisely one defaulting obligor).  There is no
   #n value passed here
   tmp1 <- integrate(f=integrand2, lower = 0, upper = 1, mu = par.ests[1], sigma = par.ests[2], kk = 1, lfunc = link);
   #In R, the first list item in the return list from function 'integrate()' is called 'value'
   #rather than 'integral' so we replace the following S-Plus line with its R equivalent:
   pi <- tmp1$value  #R-language equivalent
   #Integrate to get 2nd noncentral moment: joint probability of exactly two obligors defaulting.
   tmp2 <- integrate(f=integrand2, lower = 0, upper = 1,  mu = par.ests[1], sigma = par.ests[2], kk = 2, lfunc = link)
   pi2 <- tmp2$value;
   rhoY <- (pi2 - pi^2)/(pi - pi^2); #See equation 8.22 on p. 345 of QRM Book.
   list(par.ests = par.ests, maxloglik = loglik, converged = conv, details = convtype, pi = pi, pi2 = pi2, rhoY = rhoY);
}

###########################################################################
#fit.binomialLogitnorm() differs from fit.binomialProbitnorm() only in the link function.  Here the 
#link funtion is pnorm() where in the logit version, it is 1/(1-exp(-z))).
#The argument M is a vector like 'Bdefaults' (the number of defaults) and m is a vector like 
#'Bobligors' (the number of obligors).  This function calls optim(...method="L-BFGS-B") so you may need to
#reset lower and upper parameter estimators if convergence occurs at an endpoint of either limit and run
#another test.
fit.binomialLogitnorm <- function(M, m, startvals = c(-1, 0.5),lowerParamLimits = c(-5.0, 0.02),
  upperParamLimits = c(1,0.9))
{
   #ARGUMENTS:
   #  M                      vector of number of successes (defaults)
   #  m                      vector of number of trials (obligors) 
   # startvals               values at which to start parameter estimation process
   # lowerParamLimits        feasible limits for parameter estimates-rerun if outputs equal these limits
   # upperParamLimits        feasible limits for parameter estimates-rerun if outputs equal these limits
   # link                    The link function here is the so-called logistic function
   
   link <- function(z)
   {
	1/(1 + exp( - z))
   }

   #This local function name will be passed to negloglik() as the 3rd parameter. The integrand will be
   #calculated as a function of p and has five other parameters (kk,n ,mu, sigma, lfunc where lfunc is
   #a link function passed in.
   integrand <- function(p, kk, n, mu, sigma, lfunc)
   {
      exp(kk * log(lfunc(mu + sigma * qnorm(p))) + (n - kk) * log(1 - lfunc(mu + sigma * qnorm(p))));
   }
   
   #Renamed this to integrand2; it is called directly from the code of this function and not indirectly 
   #through an internal function.  It is a function of p and has only four parameters. lfunc is name 
   #of link function passed in as one of the parameters. 
   #It assumes probability of exactly kk obligors defaulting:
   integrand2 <- function(p, mu, sigma, kk, lfunc)
   {
      exp(kk * log(lfunc(mu + sigma * qnorm(p))));
   }

   #10/9/2007: removed assign() statements.
   #Combine the defaults M in column 1 and the number of obligors m in column 2 of a matrix:
   MLdata <- cbind(M,m); 

   #10/9/2007: added parameters which replace assign(). We pass integrand which integrates over p
   #with 5 parameters including a link function:
   negloglik <- function(theta, defaultData, FNIntegrate=integrand, FNlink=link)
   {
      #10/9/2007: added tests to make sure some parameters are functions:
      if(!is.function(FNIntegrate)) stop("parameter FNIntegrate must be a function");
      if(!is.function(FNlink)) stop("parameter FNlink must be a function");

      n <- dim(defaultData)[1]; #determine the number of rows in defaultData
      #numeric() creates a vector of length n with all values equal to 0.  The sum of the elements of this vector 
      #will be the return value of negloglik()
      out <- numeric(n); 
      for(i in 1:n) 
      {
         #R version of integrate() returns list with $value rather than the S-Plus $integral.  
         #5-parameter version integrating over p: 
         out[i] <-  - logb(integrate(FNIntegrate, lower = 0., upper = 1.,  kk = defaultData[i, 1], 
                      n = defaultData[i, 2], mu = theta[1], sigma = theta[2]^2, lfunc = FNlink)$value)
      }
      sum(out);
   }
 
   theta <- c(startvals[1], sqrt(startvals[2]))
   #SU: Replace S-Plus nlmin(f,x) calls by calls to R-function optim().
   #out <- nlmin(negloglik, theta, max.fcal = 1000., max.iter = 1000.) #S-Plus
   #par.ests <- out$x
   #These are the replacement R-language lines ('negloglik' is an input function name):
   #Set lower and upper bounds on parameters: startvals are c(-1, 0.5). Hence use L-BFGS-B method rather than BFGS method.
   #10/9/2007: added parameters to optim() which need to be passed to negloglik:
   optimout <- optim(theta,negloglik,defaultData=MLdata, FNIntegrate=integrand, FNlink=link,
           method="L-BFGS-B", lower=lowerParamLimits ,upper=upperParamLimits)
   par.ests <- optimout$par

   if(optimout$convergence == 0)
   {
      conv = TRUE
      convtype = "Convergence succeeded via L-BFGS-B quasi-Newton method"
   }
   else
   {
      conv = FALSE 
      convtype = "Convergence failed via L-BFGS-B quasi-Newton method"
   }
   #10/9/2007: added further parameters to negloglik:
   loglik <-  - negloglik(par.ests, defaultData=MLdata, FNIntegrate=integrand, FNlink=link);
   par.ests[2] <- par.ests[2]^2;

   #integrate() performs NUMERICAL INTEGRATION. See pp. 354-5 in QRM book for a description (Examples 8.12
   #and 8.13. Integrate to get first moment (probability of precisely one defaulting obligor).  There is no
   #n value passed here
   tmp1 <- integrate(f=integrand2, lower = 0, upper = 1, mu = par.ests[1], sigma = par.ests[2], kk = 1, lfunc = link);
   #In R, the first list item in the return list from function 'integrate()' is called 'value'
   #rather than 'integral' so we replace the following S-Plus line with its R equivalent:
   pi <- tmp1$value  #R-language equivalent
   #Integrate to get 2nd noncentral moment: joint probability of exactly two obligors defaulting.
   tmp2 <- integrate(f=integrand2, lower = 0, upper = 1,  mu = par.ests[1], sigma = par.ests[2], kk = 2, lfunc = link)
   pi2 <- tmp2$value;
   rhoY <- (pi2 - pi^2)/(pi - pi^2); #See equation 8.22 on p. 345 of QRM Book.
   list(par.ests = par.ests, maxloglik = loglik, converged = conv, details = convtype, pi = pi, pi2 = pi2, rhoY = rhoY);
}

#################################################

 lbeta <- function(a, b)
{
	lgamma(a) + lgamma(b) - lgamma(a + b)
}


######################################################################
