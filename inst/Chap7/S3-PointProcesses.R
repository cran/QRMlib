# S-Plus script developed by Professor Alexander McNeil, A.J.McNeil@hw.ac.uk
# R-version adapted by Scott Ulman (scottulman@hotmail.com)
# QRMlib 1.4.2
# This free script using QRMLib is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

######Load the QRMlib and sp500 index data set##################
#QRMlib.pdf is a help file for the functions used by QRMlib.  It is available at
#...\Program Files\R\R-2.6.0\library\QRMlib\Docs
#If you have created the QRMBook workspace and .Rprofile  as described in QRMlib.pdf
#topics 'QRMBook-workspace' and 'profileLoadLibrary', then you may comment out the
#following line:
library(QRMlib);
#if you have previously opened the DJ data set (the Dow Jones for 30 stocks) AND saved 
#the workspace, you may comment out the following line:
data(sp500);
#################################################

plot(sp500, type="l");
sp500.nreturns <- -mk.returns(sp500);
#In R, seriesPositions() extracts the time slot from a timeSeries so change S-Plus code. Additionally, R's timeDate()
#requires a format string. Remember in R, a small 'y' indicates a two-digit year. Finally, the & indicates 'vectorized and'.
#window <- (positions(sp500.nreturns) > timeDate("12/31/95")) & (positions(sp500.nreturns) < timeDate("01/01/04"))
window <- (seriesPositions(sp500.nreturns) > timeDate("12/31/1995",format = "%m/%d/%Y")) & 
        (seriesPositions(sp500.nreturns) < timeDate("01/01/2004",format = "%m/%d/%Y"));
sp500.nreturns <- sp500.nreturns[window];
head(sp500.nreturns);
tail(sp500.nreturns);
plot(sp500.nreturns,type="l",main="S&P500 1/1/1996-12/31/2003")
 

#The following functions are contained in functionsHawkes.R. 
#Plot the 100 largest exceedances. Create an MPP (marked point process) class 
tmp <- extremalPP(sp500.nreturns,ne=100);
par(mfrow=c(2,1));
#Specifically call plot.MPP(). Since tmp is class MPP, it should plot via plot(tmp) as well.
plot.MPP(tmp);

#Be sure to graph with plot.PP instead of plot.MPP:
tmp2 <- unmark(tmp);
plot.PP(tmp2);
par(mfrow=c(1,1));

#tmp$times represents the Julian data counters for the highest 100 exceedances. Get the
#c() joins the starttime and endtime to the times vector, creating with 102 points. diff()
#then calculates the differences between each successive julian date count. Note that the
#diff()function has a default lag of 1 and takes first differences (diff = 1). Hence this
#calculates a vector of gaps or number of periods between the 100 largest exceedances
gaps <- diff(c(tmp$starttime,tmp$times,tmp$endtime));
plot(gaps,type="l");
#Plot the gaps against the exponential which should be the density if the exceedance times are
#independently distributed (i.e. not clustering).  Should get a straight line if not clustering.
QQplot(gaps,ref="exp");

#This is known as the POT-PP model (where PP means Point Process) rather than the POT-GPD model.
#Remember, we used fit.GPD before.
#fit the MARKED POINT process model first:
mod1 <- fit.POT(tmp);
#Now fit the UNMARKED POINT PROCESS model:
mod1b <- fit.POT(tmp2);
mod1;
mod1b;

#Now try to fit the SELF-EXCITING point processes:
mod2a <- fit.sePP(tmp,mark.influence=FALSE,std.errs=TRUE);
mod2b <- fit.sePP(tmp,mark.influence=TRUE,std.errs=TRUE);
mod2c <- fit.sePP(tmp,model="ETAS",mark.influence=FALSE,std.errs=TRUE);
mod2d <- fit.sePP(tmp,model="ETAS",mark.influence=TRUE,std.errs=TRUE);
par(mfrow=c(3,1));
plot(tmp2);
plot(tmp);
plot(mod2a);

stationary.sePP(mod2a);
stationary.sePP(mod2b);
stationary.sePP(mod2c);
stationary.sePP(mod2d);

#Now try to fit using the SELF-EXCITING MARKED PP:
mod3a <- fit.seMPP(tmp,mark.influence=FALSE,std.errs=TRUE);
mod3b <- fit.seMPP(tmp,mark.influence=TRUE,std.errs=TRUE);
mod3c <- fit.seMPP(tmp,model="ETAS",mark.influence=FALSE,std.errs=TRUE);
mod3d <- fit.seMPP(tmp,model="ETAS",mark.influence=TRUE,std.errs=TRUE);

mod4a <- fit.seMPP(tmp,mark.influence=FALSE,predictable=TRUE,std.errs=TRUE);
mod4b <- fit.seMPP(tmp,mark.influence=TRUE,predictable=TRUE,std.errs=TRUE);
mod4c <- fit.seMPP(tmp,model="ETAS",mark.influence=FALSE,predictable=TRUE,std.errs=TRUE);
mod4d <- fit.seMPP(tmp,model="ETAS",mark.influence=TRUE,predictable=TRUE,std.errs=TRUE);

stationary.sePP(mod3a);
stationary.sePP(mod3b);
stationary.sePP(mod3c);
stationary.sePP(mod3d);
stationary.sePP(mod4a);
