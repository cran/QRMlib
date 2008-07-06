# QRMlib: version 1.4.3
# modifies QRMlib 1.4.2 to correct parameter order in fit.GEV 
# this file is a component of QRMlib 

# Copyright (C) 2005-08 Alexander McNeil 
# R-language additions Copyright (C) 2006-2008 Scott Ulman

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
# R-language contact: Scott Ulman : scottulman@hotmail.com
# Note in the R-translations that TRUE has been substituted throughout the 
# code for T (and FALSE for F) when setting default parameter values as 
# specified in section 2.4 of an "Introduction to R" (R-intro.pdf).
# Otherwise "R CMD check -QRMlib" returns the following type of error when
# trying to run an example:
# > BiDensPlot(func=dmnorm,mu=c(0,0),Sigma=equicorr(2,-0.7))
#   Error in func(data, ...) : F used instead of FALSE
#   Execution halted
###################################################
#11/7/2007: changed first parameter name from x to q so help file will know q 
#is a quantile.  Hence changed transform to q = (q-mu)/sigma
pGumbel <- function(q, mu = 0, sigma = 1)
{
#Standard Gumbel has mu=0 and sigma=1. Sigma must be positive.
    if(sigma <= 0) stop("Scale parameter sigma must be strictly positive!");
    q <- (q - mu)/sigma;
    exp(-exp(-q));
}

###################################################
#To use Gumbel rather than standard Gumbel we must pass scale and location parameters
#as arguments; 
qGumbel <- function(p, mu = 0, sigma = 1)
{
    if(sigma <= 0) stop("Scale parameter sigma must be strictly positive!");
    #Make sure p parameter is vector of probabilities from cdf. Use logical AND on (0,1) which will force proper range:
    if(length(p[(p > 0) & (p < 1)]) < length(p)) stop("p parameter does not represent probabilities");
    mu + sigma*(-log(-log(p)));
}

#######################################################
dGumbel <- function(x, mu = 0, sigma = 1, logvalue=FALSE)
  {
    if(sigma <= 0) stop("Scale parameter sigma must be strictly positive!");
    q <- (x - mu)/sigma; 
    out <- -q-exp(-q)-log(sigma);
    if(!(logvalue)) out <- exp(out);
    out;
  }
#########################################################

rGumbel <- function(n, mu=0, sigma=1)
{
    if(sigma <= 0) stop("Scale parameter sigma must be strictly positive!");
    U <- runif(n);
    qGumbel(U,mu,sigma);
}

#########################################################

pGEV <- function(q, xi, mu = 0, sigma = 1)
{
  x <- (q-mu)/sigma
  #use the new pGumbel which passes mu and sigma and pass q rather than x:
  if (xi==0) 
  {
        return(out <- pGumbel(q, mu, sigma));
  }
  else out <- exp( - (1 + xi* x)^(-1/xi))
  if (xi  > 0) out[1+xi*x  <= 0] <- 0
  if (xi < 0) out[1+xi*x <=0] <- 1
  out
}


#######################################################

qGEV <- function(p, xi, mu = 0, sigma = 1)
{
  #use the new qGumbel which passes mu and sigma and return result immediately:
  if (xi==0) 
  {
       return(out <- qGumbel(p, mu, sigma))
  }
  else
    {
      inner <- ((-log(p))^(-xi)-1)/xi
      if (xi> 0) out <- pmax(inner,-1/xi)
      if (xi <0) out <- pmin(inner,1/abs(xi))
    }
  mu +sigma*out
}

#######################################################

dGEV <- function(x, xi, mu = 0, sigma = 1, logvalue=FALSE)
  {
    xx <- (x-mu)/sigma
    #use the new dGumbel which passes mu and sigma:
    #if (xi==0) out <- dGumbel(xx,logvalue=TRUE)-log(sigma)
    if (xi==0) {
       return(out <- dGumbel(x, mu, sigma, logvalue));
    }

    else
      { out <- rep(-Inf,length(x))
        out[1+xi*xx>0] <- (-1/xi-1)*log(1+xi*xx[1+xi*xx>0]) - (1+xi*xx[1+xi*xx>0])^(-1/xi) -log(sigma)
      }
    if (!(logvalue))
      out <- exp(out)
    out
  }

######################################################

rGEV <- function(n, xi, mu=0, sigma=1)
  {
    U <- runif(n)
    qGEV(U,xi,mu,sigma)
  }

#######################################################

fit.GEV <- function(maxima)
  {
    sigma0 <- sqrt((6. * var(maxima))/pi)
    mu0 <- mean(maxima) - 0.57722 * sigma0
    xi0 <- 0.1
    theta <- c(xi0, mu0, sigma0)
    
    #10/5/2007: removed assign() for maxima.nl
    #10/5/2007: passed additional parameter to avoid using assign()
    negloglik <- function(theta, maxvalue)
    {
      -sum(dGEV(maxvalue,theta[1],theta[2],abs(theta[3]),logvalue=TRUE));        
    }

    #10/5/2007: passed additional 4th parameter which will be passed to negloglik
    optimfit <- optim(theta, fn=negloglik, gr=NULL, maxvalue=maxima, method="BFGS");
    par.ests <- optimfit$par;
    if(optimfit$convergence == 0) #if code is 0 we are OK; 
      converged <- TRUE
    else
      converged <- FALSE;
    #replace 3rd parameter estimate with its absolute value:
    par.ests[3] <- abs(par.ests[3]);
    #10/5/2007: added final parameter which must be passed to negloglik
    fisher <- hessb(negloglik, par.ests, maxvalue=maxima);
    varcov <- solve(fisher);

    par.ses <- sqrt(diag(varcov));
    out <- list(par.ests = par.ests, par.ses = par.ses, varcov = varcov, converged = 
    #10/5/2007: passed additional 2nd parameter to -negloglik()
    converged, llmax = -negloglik(par.ests,maxima));
    #06/30/2008: the order of the names is incorrect; note that the parameter vector theta 
    #(passed to the optimization function) has the order c(xi0, mu0, sigma0); the order of 
    #the return parameters should be identical;
    #hence we need to make the following correction to maintain the same order
    #names(out$par.ests) <- c("xi", "sigma", "mu");
    #names(out$par.ses) <- c("xi", "sigma", "mu");
    names(out$par.ests) <- c("xi", "mu", "sigma");
    names(out$par.ses) <- c("xi", "mu", "sigma");
    out;
  }

##############################################################################

pGPD <- function(q, xi, beta=1){
  x <- q/beta
  if (xi==0)
    out <- pexp(x)
  else
    out <- 1 -(1+xi*x)^(-1/xi)
  out[x<0] <- 0
  if(xi<0)
    out[x > -1/xi] <- 1
  out 
}

###############################################################################

qGPD <- function(p, xi, beta=1)
{
    if (xi==0)
      out <- qexp(p)
    else
      {
        inner <- (1/xi) * ((1-p)^(-xi)-1)
        out <- pmax(inner,0)
        if (xi<0)
          out <- pmin(inner,1/abs(xi))
      }
    beta*out
}

#########################################################################

rGPD <- function(n, xi, beta=1)
{
    U <- runif(n);
    #Bugfix on 08/21/2006: must call qGPD(), not gGPD()
    #gGPD(U,xi,beta)
    qGPD(U,xi,beta);
}

#########################################################################

dGPD <- function(x, xi, beta=1, logvalue=FALSE)
  {
    xx <- x/beta
    if (xi==0)
      out <- log(dexp(xx))-log(beta)
    else{
      out <- rep(-Inf,length(x))
      cond <- (xx>0)
      if (xi<0) cond <- cond & (xx < 1/abs(xi))
      out[cond] <- (-1/xi-1)*log(1+xi*xx[cond])-log(beta)
    }
    if(!(logvalue))
      out <- exp(out)
    out
}

########################################################################
#Find a threshold for GPD estimation
findthreshold <- function(data, ne)
{
#Parameters Required:
#   data - a data set (vector or numeric) whose threshold is to be determined.  
#          Should be a vector, not a matrix, dataframe, or multiple timeSeries
#          If using a matrix, pass matname[,n] to pass nth column
#          If using a df, pass dfname[["colname"]] or dfname[[n]] or dfname$colname or 
#          dfname[ , "colname"] or dfname[ ,n] where n is col number
#          If using a timeSeries, pass "as.vector(tS@Data[,n]" to pass nth column of 
#          timeSeries data
#   ne   - scalar or vector number of exceesses to use.  
#  Return value: If passing a scalar ne, return will be the lowest data value from the vector
#  such that ne values in the vector exceed it.  If passing a vector as ne, the return value
#  will be a vector containing the threshold value corresponding to each number of exceedances
#  in the ne vector.
  #10/4/2007: Fixed for R-2.6.0. In R-2.5.0 and prior, as.numeric() could process timeSeries objects directly.  
  #Beginning with R-2.6.0 you must test for timeSeries and pass only the @Data slot to as.numeric():
  if (exists("is.R") && is.function(is.R) && is.R()) 
  {
     if(is.timeSeries(data)) data <- as.vector(data@Data); 
     if(!is.vector(data)) stop("data input to findthreshold() must be a vector or timeSeries with only one data column");
     if(length(data) < ne) stop("data length less than ne (number of exceedances");
  } 
  data <- rev(sort(as.numeric(data)));
  thresholds <- unique(data);
  indices <- match(data[ne], thresholds);
  indices <- pmin(indices + 1., length(thresholds));
  thresholds[indices];
}
##########################################################################
#This method uses numerical derivatives and the optim() function for the MLE parameters.
#Alternatively, you may use fit.GPDb() which uses the actual calculated derivative function
#in place of numerical derivatives and uses nlminb() rather than optim().
fit.GPD <- function(data, threshold = NA, nextremes = NA, method = "ml", information = "observed")
{
    if(is.na(nextremes) & is.na(threshold))
      stop("Enter either a threshold or the number of upper extremes");
	
    if(!is.na(nextremes) & !is.na(threshold))
      stop("Enter EITHER a threshold or the number of upper extremes");

    #10/4/2007: Fixed for R-2.6.0. In R-2.5.0 and prior, as.numeric() could process timeSeries objects directly.  
    #Beginning with R-2.6.0 you must test for timeSeries and pass only the @Data slot to as.numeric():
    if (exists("is.R") && is.function(is.R) && is.R()) 
    {
       if(is.timeSeries(data)) data <- data@Data; 
    } 

    data <- as.numeric(data);
    n <- length(data);
	
    if(!is.na(nextremes))
    {
      #10/4/2007: Substitute QRM findthreshold() for S-Plus findthresh()
      # R complains: no visible global function for 'findthresh' despite fact that R will ALWAYS EXECUTE if statement
      #and S-plus will execute the else statement
      #if (exists("is.R") && is.function(is.R) && is.R()){threshold <- findthreshold(data, nextremes);}
      #else {threshold <- findthresh(data, nextremes);} #Splus
      threshold <- findthreshold(data, nextremes);
    }
    exceedances <- data[data > threshold];
    excess <- exceedances - threshold;
    Nu <- length(excess);
    if(method == "ml") 
    {
       xbar <- mean(excess);
       a0 <- xbar;
       gamma <- -0.35;
       delta <- 0.;
       pvec <- ((1.:Nu) + delta)/(Nu + delta);
       a1 <- mean(sort(excess) * (1. - pvec))
       xi <- 2. - a0/(a0 - 2. * a1);
       beta <- (2. * a0 * a1)/(a0 - 2. * a1);
       par.ests <- c(xi,beta);

       #10/5/2007: added 2nd parameter to internal negloglik to remove assign():
       negloglik <- function(theta, excesses)
       {
         -sum(dGPD(excesses,theta[1],abs(theta[2]),logvalue=TRUE))
       }

       #added parameters to pass to negloglik in 3rd position to remove assign(): the is 
       #the new 
       optimfit <- optim(par.ests, fn=negloglik, gr=NULL, excesses=excess, method="BFGS");
       par.ests <- optimfit$par;
       #replace the 2nd parameter with its absolute value:
       par.ests[2] <- abs(par.ests[2]);
       if(optimfit$convergence == 0) #if code is 0 we are OK; 
         converged <- TRUE
       else
         converged <- FALSE;

       #10/5/2007: added 2nd parameter to negloglik to remove assign(): 
       ll.max <- -negloglik(optimfit$par, excesses=excess); 
  
       if(information == "observed") 
       {
         #10/5/2007: passed 3rd parameter to hessb(); it is the excesses parameter to negloglik()
         fisher <- hessb(negloglik, optimfit$par, excesses = excess);
         varcov <- solve(fisher);
       }
       if(information == "expected") 
       {
         one <- (1. + par.ests[1.])^2./Nu;
         two <- (2. * (1. + par.ests[1.]) * par.ests[2.]^2.)/Nu;
         cov <-  - ((1. + par.ests[1.]) * par.ests[2.])/Nu;
         varcov <- matrix(c(one, cov, cov, two), 2.);
       }
       
    }

    if(method == "pwm") 
    {
      xbar <- mean(excess);
      a0 <- xbar;
      gamma <- -0.35;
      delta <- 0.;
      pvec <- ((1.:Nu) + delta)/(Nu + delta);
      a1 <- mean(sort(excess) * (1. - pvec));
      xi <- 2. - a0/(a0 - 2. * a1);
      beta <- (2. * a0 * a1)/(a0 - 2. * a1);
      par.ests <- c(xi, beta);
      denom <- Nu * (1. - 2. * xi) * (3. - 2. * xi);
      if(xi > 0.5) 
      {
        denom <- NA;
        warning("Asymptotic standard errors not available for PWM Method when xi > 0.5");
      }
      one <- (1. - xi) * (1. - xi + 2. * xi^2.) * (2. - xi)^2.;
      two <- (7. - 18. * xi + 11. * xi^2. - 2. * xi^3.) * beta^2.;
      cov <- beta * (2. - xi) * (2. - 6. * xi + 7. * xi^2. - 2. *xi^3.);
      varcov <- matrix(c(one, cov, cov, two), 2.)/denom;
      information <- "expected";
      converged <- NA;
      ll.max <- NA;
    }
    par.ses <- sqrt(diag(varcov));
    p.less.thresh <- 1. - Nu/n;
	
    out <- list(n = length(data), data = exceedances, threshold = threshold,
		p.less.thresh = p.less.thresh, n.exceed = Nu, method = method,
		par.ests = par.ests, par.ses = par.ses, varcov = varcov, 
		information = information, converged = converged, ll.max = 
		ll.max);
    names(out$par.ests) <- c("xi", "beta");
    names(out$par.ses) <- c("xi", "beta");
    
    out;
}
######################################################################################
#This function parallels fit.GPD() except it uses nlminb() for the optimization and utilizes 
#the actual derivative function rather than calculating numerical first derivatives.
#Pass the data slot of a time series (tS@Data) in lieu of timeSeries.
fit.GPDb <- function(data, threshold = NA, nextremes = NA, method = "ml", information = "observed")
{
   if(is.na(nextremes) & is.na(threshold))
     stop("Enter either a threshold or the number of upper extremes");
		
   if(!is.na(nextremes) & !is.na(threshold))
     stop("Enter EITHER a threshold or the number of upper extremes");

    #10/4/2007: Fixed for R-2.6.0. In R-2.5.0 and prior, as.numeric() could process timeSeries objects directly.  
    #Beginning with R-2.6.0 you must test for timeSeries and pass only the @Data slot to as.numeric():
    if (exists("is.R") && is.function(is.R) && is.R()) 
    {
       if(is.timeSeries(data)) data <- data@Data; 
    } 

    data <- as.numeric(data);
    n <- length(data);

    if(!is.na(nextremes))
    {
      #10/4/2007: Substitute QRM findthreshold() for S-Plus findthresh().  
      # R complains: no visible global function for 'findthresh' despite fact that R will ALWAYS EXECUTE if statement
      #and S-plus will execute the else statement
      #if (exists("is.R") && is.function(is.R) && is.R()){threshold <- findthreshold(data, nextremes);}  #R-language
      #else {threshold <- findthresh(data, nextremes);} #Splus
      threshold <- findthreshold(data, nextremes);
    }
		
    exceedances <- data[data > threshold];
    excess <- exceedances - threshold;
    Nu <- length(excess);
    if(method == "ml") 
    {
       xbar <- mean(excess);
       a0 <- xbar;
       gamma <- -0.35;
       delta <- 0.;
       pvec <- ((1.:Nu) + delta)/(Nu + delta);
       a1 <- mean(sort(excess) * (1. - pvec));
       xi <- 2. - a0/(a0 - 2. * a1);
       beta <- (2. * a0 * a1)/(a0 - 2. * a1);
       par.ests <- c(xi,beta);
          
       #SU: COMMENT: this (negloglik()) function is called from nlminb()(a nonlinear minimization
       # subject to Box constraints).  It appears the R-function for nlminb() is very similar to 
       #the S-Plus function.  The Ydata parameter to nlminb() is an argument passed both to the 
       #objective (negloglik()) function (the 2nd parameter passed to nlminb() and the derivative
       #function.  
       negloglik <- function(theta,Ydata)
       {
          -sum(dGPD(Ydata,theta[1],abs(theta[2]),logvalue=TRUE))
       }
       #Similarly, this deriv() function is called  from nlminb()in lieu of calculating numerical
       #derivatives  The Ydata parameter is also an input to nlminb() which calls this function. 
       deriv <- function(theta,Ydata)
       {
          xi <- theta[1];
          beta <- theta[2];
          term1 <- sum(Ydata/(beta+xi*Ydata));
          term2 <- sum(log(1+xi*Ydata/beta));
          d1 <- -term2*xi^(-2)+(1+1/xi)*term1;
          d2 <- (length(Ydata)-(xi+1)*term1)/beta;
          c(d1,d2);
       }
       options(warn=-2);
       #Ydata is the 2nd parameter passed to both negloglik() and gradient() when nlminb()
       #calls them; the first parameter sent to both functions is the first parameter to nlminb(),
       #i.e. the 'start=par.ests' parameter value.
       fit <- nlminb(start=par.ests, objective=negloglik,gradient=deriv, Ydata=excess);
       #In R, the fit parameters returned by list from nlminb() are called $par
       #par.ests <- fit$parameters
       #Reset par.ests to the fitted values:          
       par.ests <- fit$par
       #The 'beta' parameter must be greater than 0 so convert via absolute value
       par.ests[2] <- abs(par.ests[2])
       #R has an output parameter named 'convergence' which is 0 if convergence occurs
       #converged <- fit$message
       if(fit$convergence == 0) #if code is 0 we are OK 
          converged <- TRUE
        else
           converged <- FALSE;

       ll.max <- -fit$objective;
       if(information == "observed") 
       {
          fisher <- hessb(negloglik, par.ests, ep=0.0001, Ydata=excess);
          varcov <- solve(fisher);
       }
       if(information == "expected") 
       {
          one <- (1. + par.ests[1.])^2./Nu;
          two <- (2. * (1. + par.ests[1.]) * par.ests[2.]^2.)/Nu;
          cov <-  - ((1. + par.ests[1.]) * par.ests[2.])/Nu;
          varcov <- matrix(c(one, cov, cov, two), 2.);
       }
    }
    if(method == "pwm") 
    {
	xbar <- mean(excess);
       a0 <- xbar;
       gamma <- -0.35;
       delta <- 0.;
       pvec <- ((1.:Nu) + delta)/(Nu + delta);
       a1 <- mean(sort(excess) * (1. - pvec));
       xi <- 2. - a0/(a0 - 2. * a1);
       beta <- (2. * a0 * a1)/(a0 - 2. * a1);
       par.ests <- c(xi, beta);
       denom <- Nu * (1. - 2. * xi) * (3. - 2. * xi);
       if(xi > 0.5) 
       {
           denom <- NA;
           warning("Asymptotic standard errors not available for PWM Method when xi > 0.5");
       }
       one <- (1. - xi) * (1. - xi + 2. * xi^2.) * (2. - xi)^2.;
       two <- (7. - 18. * xi + 11. * xi^2. - 2. * xi^3.) * beta^2.;
       cov <- beta * (2. - xi) * (2. - 6. * xi + 7. * xi^2. - 2. *xi^3.);
       varcov <- matrix(c(one, cov, cov, two), 2.)/denom;
       information <- "expected";
       converged <- NA;
	ll.max <- NA;
    }
    par.ses <- sqrt(diag(varcov));
    p.less.thresh <- 1. - Nu/n;
    out <- list(n = length(data), data = exceedances, threshold = threshold,
		p.less.thresh = p.less.thresh, n.exceed = Nu, method = method,
		par.ests = par.ests, par.ses = par.ses, varcov = varcov, 
		information = information, converged = converged, ll.max = 
		ll.max);
    names(out$par.ests) <- c("xi", "beta");
    names(out$par.ses) <- c("xi", "beta");
    out;
}

############################################################################################
#This method will be called from showRM() but you may call it yourself directly IF you
#have called fit.GPD() first since the first parameter is the return value from fit.GPD
plotTail <- function(object, extend=2,fineness=1000,...)
{
#Parameters Required:
#   object - the return value from a call to fit.GPD which is a list like this:
#            list(n = length(data), data = exceedances, threshold = threshold,
#		p.less.thresh = p.less.thresh, n.exceed = Nu, method = method,
#		par.ests = par.ests, par.ses = par.ses, varcov = varcov, 
#		information = information, converged = converged, ll.max = 
#		ll.max)

   data <- as.numeric(object$data);
   threshold <- object$threshold;
   xi <- object$par.ests[names(object$par.ests) == "xi"];
   beta <- object$par.ests[names(object$par.ests) == "beta"];
   xpoints <- sort(data);
   ypoints <- ppoints(sort(data));
   xmax <- max(xpoints)*extend;
   prob <- object$p.less.thresh;
   ypoints <- (1- prob) * (1 - ypoints);       
   x <- threshold+qGPD((0:(fineness-1))/fineness, xi, beta);
   x <- pmin(x,xmax);
   y <- pGPD(x-threshold, xi,beta);
   y <- (1- prob) * (1 - y);        
   plot(xpoints, ypoints, xlim=range(threshold,xmax), ylim = range(ypoints, y), xlab = "x (on log scale)", ylab = "1-F(x) (on log scale)", log="xy",...)
        lines(x,y);
   NULL;
}

############################################################################
# Show Risk-Measure Estimates on Tailplot
#Parameters Required:
#   object - the return value from a call to fit.GPD which is a list like this:
#            list(n = length(data), data = exceedances, threshold = threshold,
#		p.less.thresh = p.less.thresh, n.exceed = Nu, method = method,
#		par.ests = par.ests, par.ses = par.ses, varcov = varcov, 
#		information = information, converged = converged, ll.max = 
#		ll.max)
#   alpha - the probability level (a high value like .99)
# Parameters - Optional
#   RM - the risk measure, either VaR or ES
#   extend - how far to extend picture; x-axis extends to this value times the largest observation 
#   ci.p - confidence level for confidence interval 
#   like.num - number of evaluations of profile likelihood 
showRM <- function(object, alpha, RM="VaR", extend=2, ci.p = 0.95, like.num = 50.)
{	
      threshold <- object$threshold;
      par.ests <- object$par.ests;
      xihat <- par.ests[names(par.ests) == "xi"];
      betahat <- par.ests[names(par.ests) == "beta"];
      p.less.thresh <- object$p.less.thresh;
      a <- (1-alpha)/(1 - p.less.thresh);
      quant <- threshold+betahat*(a^(-xihat)-1)/xihat;
      es <- quant/(1-xihat) + (betahat-xihat*threshold)/(1-xihat);
      #in the following switch, if RM is "VaR:, the 2nd argument quant is selected; otherwise if "ES", then es is selected;
      #otherwise an error occurs.
      point.est <- switch(RM,VaR=quant,ES=es);
      plotTail(object,extend=2);
      abline(v=point.est);
      #stretch the x-axis out by multiplying the largest data item by the 4th parameter, which defaults to 2
      xmax <- max(object$data)*extend;
      #10/5/2007: eliminated assign() by adding parameters to local function called internally from showRM() function
      parloglik <- function(theta, excessesIn, xpiIn, aIn, uIn, RMIn)
      {
         xi <- theta;
         if (RMIn=="VaR")
           beta <- xi*(xpiIn-uIn)/(aIn^(-xi)-1);
         if (RMIn=="ES")
           beta <- ((1 - xi) * (xpiIn - uIn))/(((aIn^( - xi) - 1)/xi) +1);
         #Must arbitrarily set 'out' to 1.0e17 (a large value) rather than NA.  
         #Otherwise we get the error message "Error in optim(xihat,...)Initial value in 'vmmin'
         # is not finite." Hence replace NA with 1.0e17. We cannot use Inf which gives same error.
         #Furthermore, if beta == 0, the dGPD() will overflow to NaN in 'else' clause and the same 
         #optim() error with infinite 'vmmin' will occur.  Hence change the original 'if(beta<0)' 
         #to 'if(beta <= 0)' and out to 1.0e17 from NA.
         if (beta <= 0) out <- 1.0e17 #NA
         else
            out <- -sum(dGPD(excessesIn,xi,beta,logvalue=TRUE));
         out;
       }

       parmax <- NULL;
       start <- switch(RM,VaR=threshold,ES=quant);
       xp <- exp(seq(from = log(start), to = log(xmax), length = like.num));

       for(i in 1.:length(xp)) 
       {
         #10/5/2007: added parameteres to optim which must be passed to parloglik in positions 3-7
         # These parameters eliminate the need for any assign() statements.
          optimfit2 <- optim(xihat, parloglik, excessesIn=(object$data - threshold), xpiIn=xp[i], 
                       aIn=a, uIn=threshold, RMIn=RM, method="BFGS");
          #build a matrix by adding the parameter rows from each iteration of -parloglik(): 
          parmax <- rbind(parmax, -parloglik(optimfit2$par, excessesIn=(object$data - threshold), 
                       xpiIn=xp[i], aIn=a, uIn=threshold, RMIn=RM));
        }
        #set 'overallmax' to the max from the call to fit.GPD:
        overallmax <-  object$ll.max;
        crit <- overallmax - qchisq(0.999, 1)/2.;
        cond <- parmax > crit;
        #Get reduced size vectors by copying only those elements which satisfy 'parmax > crit'
        #into a vector of the same name; hence xp[] and parmax[] will have fewer than 50 elements:
        xp <- xp[cond];
        parmax <- parmax[cond];  
        #The following 'par()' cmd says the "next high-level plot() command should NOT  
        #clean the frame before drawing as if it were on a new device" according to HELP  
        par(new = TRUE);
        plot(xp, parmax, type = "n", xlab = "", ylab = "", axes = FALSE,
			ylim = range(overallmax, crit), xlim=range(threshold,xmax), log = "x");
        #Add an axis to the current plot:
        #10/5/2007: Fixed for R-2.6.0. R uses "tick=" as parameter rather than "ticks=". Added if/else.
        if (exists("is.R") && is.function(is.R) && is.R()) 
          {axis(4., at = overallmax - qchisq(c(0.95, 0.99), 1.)/2., labels=c("95", "99"), tick = TRUE)}
        else {axis(4., at = overallmax - qchisq(c(0.95, 0.99), 1.)/2., labels=c("95", "99"), ticks = TRUE)};
        aalpha <- qchisq(ci.p, 1.);
        abline(h = overallmax - aalpha/2, lty = 2, col = 2);
        #Introduce a new condition to replace parmax > crit:
        cond <- !is.na(xp) & !is.na(parmax);
        #perform cubic spline interpolation of given points
        #xp and parmax have already been reduced in size. Here they are potentially reduced again:
        #any xp, parmax pair containing valid numbers (both not NA) is included)
        smth <- spline(xp[cond], parmax[cond], n = 200.);
        #Add connected line segments to plot: smth is a vector containing coordinates to join
        lines(smth, lty = 2., col = 2.);
        ci <- smth$x[smth$y > overallmax - aalpha/2.];
        out <- c(min(ci), point.est, max(ci));
        names(out) <- c("Lower CI", "Estimate", "Upper CI");
	out
}
################################################################################
MEplot <- function(data, omit = 3., labels = TRUE, ...)
{
  #10/4/2007: Fixed for R-2.6.0. In R-2.5.0 and prior, as.numeric() could process timeSeries objects directly.  
  #Beginning with R-2.6.0 you must test for timeSeries and pass only the @Data slot to as.numeric():
  if(is.timeSeries(data)) data <- data@Data;   
  data <- as.numeric(data);
  n <- length(data);
  #The rank() function "averages ties", so if ranks 4 and 5 have equal values, then they are both listed
  #as 4.5, a decimal.  Hence rank() returns non-integer values. Function myrank() does the same thing as rank except 
  #it returns INTEGER VALUES for ranks.  TIES are all given the HIGHEST integer value. Hence if ranked index values
  #3,4 are tied, EACH receives a 4 so there are is 3 index listed nor are there any decimal values.
  myrank <- function(x, na.last = TRUE)
  {
     #sort.list() will give the index values for the list sorted from smaller to bigger.
     #Applying sort.list() to sort.list() will rearrange to give the index in the 1st sort.list  
     #Example.  Suppose 1st sort list returns
     # [1] 13  9 25  1 30 14  3 29  4 20 18 27  2  8 16 26 10 12 21 19  5 23 11  7  6
     #[26] 15 24 28 22 17  so the lowest index is in position 4.  Then sort.list(sort.list())
     #will show 4 as its first value since the 4th value in 1st sort list is 1 (the lowest value in the original list.   
     ranks <- sort.list(sort.list(x, na.last = na.last));
     #If the USER OVERRIDES the na.last=TRUE PARAMETER, THIS TEST IS NEEDED; OTHERWISE NO.
     if(is.na(na.last))
        #Is there an R-function named is.orderable()? No.  In SPlus, is.orderable() returns !is.na(x) 
        x <- x[!is.na(x)];  #SPlus uses x <- x[is.orderable(x)]
     for(i in unique(x[duplicated(x)])) 
     {
       which <- x == i & !is.na(x);
       ranks[which] <- max(ranks[which]);
     }
     ranks;
  }
  data <- sort(data);
  #unique() returns a vector with duplicates removed:
  n.excess <- unique(floor(length(data) - myrank(data)));
  points <- unique(data);
  nl <- length(points);
  n.excess <- n.excess[ - nl];
  points <- points[ - nl];
  excess <- cumsum(rev(data))[n.excess] - n.excess * points;
  y <- excess/n.excess;
  plot(points[1.:(nl - omit)], y[1.:(nl - omit)], xlab = "", ylab = "",...);
  if(labels)
    title(xlab = "Threshold", ylab = "Mean Excess");
}

##############################################################
#calculates risk measures like VaR and expected shortfall based on a 
#generalized Pareto model fitted to losses over a high threshold
#Output - matrix with quantile and shortfall estimates for each probability level 
#Parameters: 
#    out - Results of a GPD fit to excesses over high thresholds, i.e the list return from
#          gpd.fit()
#    p - vector of probability measures for risk levels (e.g. c(.99,.995)  
# Return Value - matrix with quantile and shortfall estimates for each probability level in p
RiskMeasures <- function(out, p)
{
	u <- out$threshold
	par.ests <- out$par.ests
	xihat <- par.ests[names(par.ests) == "xi"]
	betahat <- par.ests[names(par.ests) == "beta"]
	p.less.thresh <- out$p.less.thresh
	lambda <- 1./(1. - p.less.thresh)
	quant <- function(pp, xi, beta, u, lambda)
	{
		a <- lambda * (1. - pp)
		u + (beta * (a^( - xi) - 1.))/xi
	}
	short <- function(pp, xi, beta, u, lambda)
	{
		a <- lambda * (1. - pp)
		q <- u + (beta * (a^( - xi) - 1.))/xi
		(q * (1. + (beta - xi * u)/q))/(1. - xi)
	}
	q <- quant(p, xihat, betahat, u, lambda)
	es <- short(p, xihat, betahat, u, lambda)
	cbind(p, quantile = q, sfall = es)
}

###################################################
#create a plot showing how the estimate of GPD shape varies with threshold or number of extremes. 
xiplot <- function(data, models = 30., start = 15., end = 500., reverse = TRUE, ci = 
	0.95, auto.scale = TRUE, labels = TRUE, table = FALSE, ...)
{
       #10/4/2007: Fixed for R-2.6.0. In R-2.5.0 and prior, as.numeric() could process timeSeries objects directly.  
       #Beginning with R-2.6.0 you must test for timeSeries and pass only the @Data slot to as.numeric():
       if (exists("is.R") && is.function(is.R) && is.R()) 
       {
         if(is.timeSeries(data)) data <- data@Data; 
       } 
	data <- as.numeric(data);
	qq <- 0.;
	if(ci)
		qq <- qnorm(1. - (1. - ci)/2.);
	x <- trunc(seq(from = min(end, length(data)), to = start, length = 
		models));
	gpd.dummy <- function(nex, data)
	{
		out <- fit.GPD(data = data, nex = nex, information = "expected");
		c(out$threshold, out$par.ests[1.], out$par.ses[1.]);
	}
	mat <- apply(as.matrix(x), 1., gpd.dummy, data = data);
	mat <- rbind(mat, x);
	dimnames(mat) <- list(c("threshold", "shape", "se", "exceedances"),
		NULL);
	thresh <- mat[1.,  ];
	y <- mat[2.,  ];
	yrange <- range(y);
	if(ci) {
		u <- y + mat[3.,  ] * qq;
		l <- y - mat[3.,  ] * qq;
		yrange <- range(y, u, l);
	}
	index <- x;
	if(reverse)
		index <-  - x;
	if(auto.scale)
		plot(index, y, ylim = yrange, type = "l", xlab = "", ylab = "",
			axes = FALSE, ...)
	else plot(index, y, type = "l", xlab = "", ylab = "", axes = FALSE, ...);
       #10/4/2007: Fixed for R-2.6.0. R uses "tick=" as parameter rather than "ticks=". Added if/else.
       if (exists("is.R") && is.function(is.R) && is.R()) 
       {
	  axis(1., at = index, labels = paste(x), tick = FALSE);
	  axis(2.);
	  axis(3., at = index, labels = paste(format(signif(thresh, 3.))), tick = FALSE);
       }
       else
       {
         axis(1., at = index, lab = paste(x), ticks = FALSE);
         axis(2.);
         axis(3., at = index, lab = paste(format(signif(thresh, 3.))), ticks = FALSE);
       }
	box();
	if(ci) {
		lines(index, u, lty = 2., col = 2.);
		lines(index, l, lty = 2., col = 2.);
	}
	if(labels) {
		labely <- "Shape (xi)";
		if(ci)
			labely <- paste(labely, " (CI, p = ", ci, ")", sep = "");
		title(xlab = "Exceedances", ylab = labely);
		mtext("Threshold", side = 3., line = 3.);
	}
	if(table)
		print(mat);
	NULL;
}

###################################################
#This code originally appeared in EVIS library for S-Plus by Alexander McNeil,
#the main author of the S-Plus code for QRMlib. (It did not exist in QRMlib.)
#It was ported to evir library in R by Alexander Stephenson.  It was added here
#by Scott Ulman only to make it easier for QRMlib users who now won't have to 
#attach library(evir) to build the Hill plot.
hillPlot <- function (data, option = c("alpha", "xi", "quantile"), start = 15, 
    end = NA, reverse = FALSE, p = NA, ci = 0.95, auto.scale = TRUE, 
    labels = TRUE, ...) 
{
    #10/4/2007: Fixed for R-2.6.0. In R-2.5.0 and prior, as.numeric() could process timeSeries objects directly.  
    #Beginning with R-2.6.0 you must test for timeSeries and pass only the @Data slot to as.numeric():
    if (exists("is.R") && is.function(is.R) && is.R()) 
    {
      if(is.timeSeries(data)) data <- data@Data; 
    } 

    data <- as.numeric(data);
    ordered <- rev(sort(data));
    ordered <- ordered[ordered > 0];
    n <- length(ordered);
    option <- match.arg(option);
    if ((option == "quantile") && (is.na(p))) 
        stop("Input a value for the probability p");
    if ((option == "quantile") && (p < 1 - start/n)) 
    {
        cat("Graph may look strange !! \n\n");
        cat(paste("Suggestion 1: Increase `p' above", format(signif(1 - 
            start/n, 5)), "\n"));
        cat(paste("Suggestion 2: Increase `start' above ", ceiling(length(data) * 
            (1 - p)), "\n"));
    }
    k <- 1:n;
    loggs <- logb(ordered);
    avesumlog <- cumsum(loggs)/(1:n);
    xihat <- c(NA, (avesumlog - loggs)[2:n]);
    alphahat <- 1/xihat;
    y <- switch(option, alpha = alphahat, xi = xihat, quantile = ordered * 
        ((n * (1 - p))/k)^(-1/alphahat));
    ses <- y/sqrt(k);
    if (is.na(end)) 
        end <- n;
    x <- trunc(seq(from = min(end, length(data)), to = start));
    y <- y[x];
    ylabel <- option;
    yrange <- range(y);
    if (ci && (option != "quantile")) 
    {
        qq <- qnorm(1 - (1 - ci)/2);
        u <- y + ses[x] * qq;
        l <- y - ses[x] * qq;
        ylabel <- paste(ylabel, " (CI, p =", ci, ")", sep = "");
        yrange <- range(u, l);
    }
    if (option == "quantile") 
        ylabel <- paste("Quantile, p =", p);
    index <- x;
    if (reverse) 
        index <- -x;
    if (auto.scale) 
        plot(index, y, ylim = yrange, type = "l", xlab = "", 
            ylab = "", axes = FALSE, ...)
    else plot(index, y, type = "l", xlab = "", ylab = "", axes = FALSE, 
        ...);
    axis(1, at = index, lab = paste(x), tick = FALSE);
    axis(2);
    #10/4/2007 Substitute QRM findthreshold() for S-Plus findthresh()
    #R complains: no visible global function for 'findthresh' despite fact that R will ALWAYS EXECUTE if statement
    #and S-plus will execute the else statement
    #if (exists("is.R") && is.function(is.R) && is.R()){threshold <- findthreshold(data, x);}
    #else {threshold <- findthresh(data, x);} #Splus
    threshold <- findthreshold(data, x);

 
    axis(3, at = index, lab = paste(format(signif(threshold, 
        3))), tick = FALSE);
    box();
    if (ci && (option != "quantile")) {
        lines(index, u, lty = 2, col = 2)
        lines(index, l, lty = 2, col = 2)
    }
    if (labels) {
        title(xlab = "Order Statistics", ylab = ylabel)
        mtext("Threshold", side = 3, line = 3)
    }
    invisible(list(x = index, y = y));
}

###################################################
#This function can be used to create Figures 7.4(c) and 7.5(c) in QRM (pp. 281-2).
#It plots the fitted GPD distribution for excesses against the empirical distribution
#function for excesses from the data.
plotFittedGPDvsEmpiricalExcesses <- function(data, threshold=NA, nextremes=NA)
{
  if(is.na(nextremes) & is.na(threshold))
     stop("Enter either a threshold or the number of upper extremes");
  #10/4/2007: Fixed for R-2.6.0. In R-2.5.0 and prior, as.numeric() could process timeSeries objects directly.  
  #Beginning with R-2.6.0 you must test for timeSeries and pass only the @Data slot to as.numeric():
  if (exists("is.R") && is.function(is.R) && is.R()) 
  {
    if(is.timeSeries(data)) data <- data@Data; 
  } 

  mod <- fit.GPD(data, threshold, nextremes);
  
  #We need threshold for plot if we don't have it:
  if(!is.na(nextremes))
    #10/4/2007 Substitute QRM findthreshold() for S-Plus findthresh()
    # R complains: no visible global function for 'findthresh' despite fact that R will ALWAYS EXECUTE if statement
    #and S-plus will execute the else statement
    #if (exists("is.R") && is.function(is.R) && is.R()){threshold <- findthreshold(data, nextremes);}
    #else {threshold <- findthresh(data, nextremes);} #Splus
    threshold <- findthreshold(as.vector(data), nextremes);


  #The values mod$data gives the exceedance values which can be plotted as excesses
  #by subtracting 'threshold' from each value.  We can get the empirical cdf of these excesses.
  pECDF <- edf(mod$data);

  #Get the max data loss:
  if(is.timeSeries(data)) 
    maxVal <- max(data@Data)
  else  
    maxVal <- max(data);

  #Get the quantile vectors for the excesses over the threshold from the GPD fit:
  quantVector <- seq(threshold, maxVal, 0.25);
  pG <- pGPD(quantVector-threshold, mod$par.ests["xi"],mod$par.ests["beta"]); 
  #Plot the two graphs:
  split.screen(c(1,1));
  plot(quantVector,pG, type="l", log="x", xlab="x (log scale)", ylab="Fu(x-u)");
  #Add the ECDF data to the plot.  Indicate we DO NOT START A NEW SCREEN:
  screen(1,new=FALSE);
  plot(mod$data,pECDF,log="x", pch=19, xlab="x (log scale)", ylab="Fu(x-u)", col="blue");
  close.screen(all=TRUE);
}
###################################################



