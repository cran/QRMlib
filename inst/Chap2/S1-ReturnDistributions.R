# S-Plus script developed by Professor Alexander McNeil, mcneil@math.ethz.ch
# R-version adapted by Scott Ulman (scottulman@hotmail.com)
# This free script using QRMLib is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

######Load the QRMlib and DJ data set##################
#QRMlib.pdf is a help file for the functions used by QRMlib.  It is available at
#...\Program Files\R\R-2.2.1\library\QRMlib\Docs
#If you have created the QRMBook workspace and .Rprofile  as described in QRMlib.pdf
#topics 'QRMBook-workspace' and 'profileLoadLibrary', then you may comment out the
#following line:
library(QRMlib);
#if you have previously opened the DJ data set (the Dow Jones for 30 stocks) AND saved 
#the workspace, you may comment out the following line:
data(DJ);
#################################################

# EXPLORING RETURN DISTRIBUTIONS

# Simulation and fitting of return distributions
n <- 4000
#Generate 4000 random data points from a student-T distribution with mu=3, sigma=4, df=6
data.st <- 3+4*rt(n,df=6)
#Try to fit a student-T to the generated random data.  Do mu, sigma, nu=df match the 
#parameters passed to the random-number generator?
out <- fit.st(data.st)

#Generate random data from generalized hyperbolic distributions (NIG requires lambda = -0.5)
data.NIG <- rghyp(n,lambda=-0.5,chi=2,psi=2,mu=0,gamma=0)
out <- fit.NH(data.NIG,symmetric=T,se=T)

data.NIG <- rghyp(n,lambda=-0.5,chi=2,psi=2,mu=0.5,gamma=0.3)
out <- fit.NH(data.NIG,se=T,symmetric=F)

data.hyp <- rghyp(n,lambda=1,chi=2,psi=2,mu=0,gamma=0)
out <- fit.NH(data.hyp,case="hyp",symmetric=T,se=T)

data.hyp <- rghyp(n,lambda=1,chi=2,psi=2,mu=0.5,gamma=0.3)
out <- fit.NH(data.hyp,case="hyp",se=T,symmetric=F)

# Working in R.  
#Make returns from timeSeries (the default is log-returns). Ret.DJ is a timeSeries class.
Ret.DJ <- mk.returns(DJ)
#In R, cutSeries() method works by selecting data only between the 'to' and 'from' dates. Hence 
#we will use the remaining (cut) data from 1993-01-01 to 2000-12-31. MUST use PRIOR DAY on 'from' 
#In version 240.10068, fCalendar uses cut() rather than cutSeries() to select a subset from timeseries:
DJ30dailyTS <- cut(Ret.DJ, from="1992-12-31", to="2000-12-31");
#We cannot use the R-function 'aggregate()' since it uses a 'ts' class rather than a 'timeSeries'
#Hence call the new function aggregateWeeklySeries from functionsUtility.R
DJ30weeklyTS <- aggregateWeeklySeries(DJ30dailyTS);
 
# make matrices from timeSeries and multiply by 100 (avoids numerical problems).  Hence we are looking at
#log-returns in percentage form and throwing away dates
DJ30daily <- 100*seriesData(DJ30dailyTS)
DJ30weekly <- 100*seriesData(DJ30weeklyTS)
dim(DJ30daily)
dim(DJ30weekly)


#Extract only the Microsoft returns as 'rseries'; remember this is a vector and not a timeSeries
rseries <- DJ30daily[,"MSFT"]
#Try to fit normal, student-T, NIG, hyp:
mod.gauss <- fit.norm(rseries)
mod.t <- fit.st(rseries)
#The default case for fit.NH(() is NIG requiring lambda = -1/2.
mod.NIG <- fit.NH(rseries)
#The alterntive case is 'hyp' where lambda = 1:
mod.hyp <- fit.NH(rseries,case="hyp")

#Build a character vector holding the maximum likelihood value from each of the normal
# student-T, NIG, and hyperbolic cases using the MSFT data:
c(mod.gauss$ll.max,mod.t$ll.max,mod.NIG$ll.max,mod.hyp$ll.max)

# make a graph 
#divide area from minimum data value in MSFT series to maximum data value into 100 equally spaced points:
xvals <- seq(from=min(rseries),to=max(rseries),length=100)
yvals.gauss <- dnorm(xvals,mean=mod.gauss$mu,sd=sqrt(mod.gauss$Sigma[1,1]))
yvals.t <- dt((xvals-mod.t$par.ests[2])/mod.t$par.ests[3],df=mod.t$par.ests[1])/mod.t$par.ests[3]
yvals.NIG <- dghyp(xvals,lambda=-1/2,chi=mod.NIG$par.ests[1],psi=mod.NIG$par.ests[2],mu=mod.NIG$par.ests[3],gamma=mod.NIG$par.ests[4])
hist(rseries,nclass=30,prob=T,ylim=range(yvals.gauss,yvals.t,yvals.NIG), main="Histogram of MSFT")
lines(xvals,yvals.gauss,col=3)
lines(xvals,yvals.t,col=4)
lines(xvals,yvals.NIG,col=5)

qqnorm(rseries)
QQplot(rseries,ref="student",df=mod.t$par.ests[1])

# Make function applying multiple models to datasets
# The following function fit.models()will be called lower in the code by being a
# parameter to the function call
#    results <- t(apply(DJ30weekly,2,fit.models))
#where t() is the transpose matrix function and apply() implicitly calls fit.models()
# which is passed as an argument.
fit.models <- function(data){
  tmp0 <- fit.norm(data)
  tmp <- fit.st(data)
  tmp2 <- fit.NH(data,case="NIG")
  out <- c(tmp0$ll.max,tmp$par.ests[1],tmp$par.ses[1],tmp$converged,tmp$ll.max,tmp2$par.ests[1],tmp2$par.ests[2],tmp2$converged,tmp2$ll.max)
  names(out) <- c("llmax.gauss","nu","nu.se","conv.t","llmax.t","chi","psi","conv.nig","llmax.nig")
  out
}

# Apply fit.models to all stocks and check convergence
#Error Tests:
dimnames(DJ30weekly)[[2]] #gives names of all columns--([[1]] would give rownames
dimnames(DJ30weekly)[[2]][1] #gives name of first column in DJ30weekly
#test <- apply(DJ30weekly[,1:2],2,fit.models) #use the first two columns of data (stock returns)
#test <- apply(DJ30weekly[,1:5],2,fit.models) #use the first five columns of data (stock returns)
#t() transposes the matrix so the stock names which were column names now appear as row names.
#apply() applies the function 'fit.models()' to each item in the DJ30weekly
results <- t(apply(DJ30weekly,2,fit.models))
results <- data.frame(results)
results

#The following attaches 'results' to the datasurface, allowing us to fetch llmax.nig and llmax.t
attach(results)
#Test whether the NIG distribution outperforms the t.  Note values returned from negloglikelihood
#are negative, so the less negative, the better. Each of these is a vector for the 30 stocks so the
#result will be a vector of the differences.  If all are positive, then NIG outperforms t for all variables.
llmax.nig-llmax.t

# Pick INTEL where NIG clearly outperforms t (llmax.nig-llmax.t = 3.0106)
rseries <- DJ30weekly[,"INTC"]

#Fit the Intel weekly series using gauss,t, NIG (requiring lambda = -1/2), 
#hyp (requiring lambda = 1):
mod.gauss <- fit.norm(rseries)
mod.t <- fit.st(rseries)
mod.NIG <- fit.NH(rseries)
mod.hyp <- fit.NH(rseries,case="hyp")

#Visually compare the loglikelihoods and observe that mod.NIG is the 
#largest since all are negative
c(mod.gauss$ll.max,mod.t$ll.max,mod.NIG$ll.max,mod.hyp$ll.max)
diff <- mod.NIG$ll.max - mod.t$ll.max

# make graph
xvals <- seq(from=min(rseries),to=max(rseries),length=100)
yvals.gauss <- dnorm(xvals,mean=mod.gauss$mu,sd=sqrt(mod.gauss$Sigma[1,1]))
yvals.t <- dt((xvals-mod.t$par.ests[2])/mod.t$par.ests[3],df=mod.t$par.ests[1])/mod.t$par.ests[3]
yvals.NIG <- dghyp(xvals,lambda=-1/2,chi=mod.NIG$par.ests[1],psi=mod.NIG$par.ests[2],mu=mod.NIG$par.ests[3],gamma=mod.NIG$par.ests[4])
hist(rseries,nclass=30,prob=T,ylim=range(yvals.gauss,yvals.t,yvals.NIG))
lines(xvals,yvals.gauss,col=3)
lines(xvals,yvals.t,col=4)
lines(xvals,yvals.NIG,col=5)

