# S-Plus script developed by Professor Alexander McNeil, A.J.McNeil@hw.ac.uk
# R-version adapted by Scott Ulman (scottulman@hotmail.com)
# QRMlib 1.4.2
# This free script using QRMLib is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

######Load the QRMlib and spdata.raw data set##################
#QRMlib.pdf is a help file for the functions used by QRMlib.  It is available at
#...\Program Files\R\R-2.6.0\library\QRMlib\Docs
#If you have created the QRMBook workspace and .Rprofile  as described in QRMlib.pdf
#topics 'QRMBook-workspace' and 'profileLoadLibrary', then you may comment out the
#following line:
library(QRMlib);
#if you have previously opened the spdata.raw data set  AND saved 
#the workspace, you may comment out the following line:
data(spdata.raw);
#################################################

#MODEL RISK See especially Section 8.4.6 on p. 364 of QRM book
spdata.raw;
attach(spdata.raw);

# Select S&P B rating data for analysis.  
# Estimating default probabilities and correlations by simple means
#momest() is an internal function in functionsCredit.R to calculate moment estimators 
#for default probabilities. The first parameter input is a vector containing the number
#of defaults in each time period; the 2nd parameter input is a vector containing the
#number of credits in the group during the time period.
momest(Bdefaults,Bobligors);
#The values calculated from momest(Bdefaults,Bobligors) are the parameter estimates
#shown in Table 8.6, p.365 of QRM book under the model column labeled 'B'
#The first value returned is the probability of a single default.
pi.B <- momest(Bdefaults, Bobligors)[1];  #one obligor defaulting pi = .04896
#the second value returned is the probability of a joint default probability for two firms.
pi2.B <- momest(Bdefaults, Bobligors)[2];  #two obligors defaulting jointly pi2 = .0031265
rhoY.B <- (pi2.B-pi.B^2)/(pi.B-pi.B^2); #correlation rhoy= .015665

# Model Risk Experiment
# Calibrate a 1-Factor Creditmetrics (probitnormal) model to pi.B and pi2.B for all models:
#The following values are shown in Table 8.6, column B, row labeled 'Probit-normal'.  In other
#words, we find the probitnorm mu and sigma values which give same probabilities as momest()
probitnorm.pars <- cal.probitnorm(pi.B,pi2.B);
probitnorm.pars;

# Calibrate beta-mixing model ( ~ CreditRisk+) to pi.B and pi2.B
#The following values are shown in Table 8.6, QRM book p. 365, column B, row labeled 'Beta'
beta.pars <- cal.beta(pi.B,pi2.B);
beta.pars;

# Calibrate Calyton copula model to pi.B and pi2.B
#The following values are shown in Table 8.6, column B, QRM book p. 365,row labeled 'Clayton'.
#Note they match the 'All models' estimators.
claytonmix.pars <- cal.claytonmix(pi.B,pi2.B);
claytonmix.pars;


# Calculate cumulative probabilities of Q on a unit interval. Plot picture of tail of three mixing distributions
# This picture essentially shows large sample asymptotics
#Build 1000 equally-spaced value on unit interval as multiples of .000999; discard all values except those below 0.25
#because we want to look at the tail, i.e. Q > 0.25 via the tail function [1 - P(Q <= 0.25)]
q <- (1:1000)/1001;
q <- q[q<0.25];
#pprobitnorm() gives cumulative probabilites for random variable q defined on unit interval
#such that the probit transform of q has a normal distribution with mu=probitnorm.pars[1]
#and sigma=probitnorm.pars[2]
p.probitnorm <- pprobitnorm(q,probitnorm.pars[1],probitnorm.pars[2])
#pbeta() gives cumulative probabilites for random variable q defined on unit interval
#where the probitnorm model parameters are used as the shape parameters for a beta distribution: 
p.beta <- pbeta(q, beta.pars[1], beta.pars[2])
#pclaytonmix() gives cumulative probabilites for random variable q defined on unit interval
#which gives an exchangeable Bernoulli mixture model equivalent to a Clayton Copula model:
p.claytonmix <- pclaytonmix(q,claytonmix.pars[1],claytonmix.pars[2])

#The following is very similar to Figure 8.3, p. 364 in QRM. It would be the same if we had used q<0.35.
#'range' returns a vector containing the minimum and maximum of all the given arguments.
scale <- range((1-p.probitnorm),(1-p.beta),(1-p.claytonmix))
plot(q, (1 - p.probitnorm), type = "l", log = "y", xlab = "q", ylab = "P(Q>q)",ylim=scale, lty="solid")
lines(q, (1 - p.beta), col = 2, lty="dashed")
lines(q, (1 - p.claytonmix), col = 3, lty="dotted")
abline(h = 0.01)
legend(0.05, 1e-4, c("Probit-normal", "Beta", "Clayton-Mixture"), lty=c(1,2,3),col = (1:3))

# We could also look at mixing densities. Remember that density values for continuous variables may
#exceed 1 since they give an approximation for the change in the cdf value as we change the x value.
#Hence if the cdf increases by 0.2 as we increase x from .1 to .2, the density should be about 2.0 (dF(x)/dx).
d.probitnorm <- dprobitnorm(q,probitnorm.pars[1],probitnorm.pars[2]);
d.beta <- dbeta(q, beta.pars[1], beta.pars[2]);
d.claytonmix <- dclaytonmix(q,claytonmix.pars[1],claytonmix.pars[2]);

plot(q, d.probitnorm, type = "l", xlab = "q", ylab = "f(q)", main="Density functions");
lines(q, d.beta, col = 2, lty="dashed");
lines(q, d.claytonmix, col = 3, lty="dotted");
legend(0.15, 15, c("Probit-normal", "Beta", "Clayton"), lty=c(1,2,3),col = (1:3));


# Simulation of 1000 years of default data; 500 firms in each year
# All models parameterised to pi.B and pi2.B from the S&P500 data over 20 years
# Models: CreditRisk+, CreditMetrics, Clayton copula
n <- 1000;
m <- rep(500,n);
#m represents a vector containing the "number of coin flips or firms in each year"
#The return value from rbinomial.mixture() is a vector with the "number of successes" in each year (from
#the 500 flips each year)
M.beta <- rbinomial.mixture(n,m,"beta",shape1=beta.pars[1],shape2=beta.pars[2]);
M.probitnorm <- rbinomial.mixture(n,m,"probitnorm",mu=probitnorm.pars[1],sigma=probitnorm.pars[2]);
M.claytonmix <- rbinomial.mixture(n,m,"claytonmix",pi=claytonmix.pars[1],theta=claytonmix.pars[2]);
#Test the moment estimators for each distribution.  momest() returns a vector of size 10 giving 
#successively the probability for one firm default, joint probability of two firms defaulting, joint
#probability of three firms defaulting,...joint probability of ten firms defaulting.
momest(M.beta,m);
momest(M.probitnorm,m);
momest(M.claytonmix,m);

# Simulate two Models with the same Asset Correlation
# Models: CreditMetrics and t latent variable model
#probitnorm.pars contains the parameter estimates for the 1-factor CreditMetrics probitnorm model.
#probitnorm.pars[3] contains the estimated correlation coefficient, i.e. the rho.asset.
probitnorm.pars;
M.tcopulamix <- rbinomial.mixture(n,m,"tcopulamix",pi=pi.B,rho.asset=probitnorm.pars[3],nu=10);
momest(M.tcopulamix,m);
pvals <- c(0.5,0.95,0.99,0.995,0.999);
cbind(quantile(M.probitnorm,pvals),quantile(M.tcopulamix,pvals));

