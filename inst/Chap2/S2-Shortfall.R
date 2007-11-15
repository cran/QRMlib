# S-Plus script developed by Professor Alexander McNeil, A.J.McNeil@hw.ac.uk
# R-version adapted by Scott Ulman (scottulman@hotmail.com)
# QRMlib 1.4.2
# This free script using QRMLib is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

######Load the QRMlib ##################
#QRMlib.pdf is a help file for the functions used by QRMlib.  It is available at
#...\Program Files\R\R-2.6.0\library\QRMlib\Docs
#If you have created the QRMBook workspace and .Rprofile  as described in QRMlib.pdf
#topics 'QRMBook-workspace' and 'profileLoadLibrary', then you may comment out the
#following line:
library(QRMlib);

# SHORTFALL TO QUANTILE RATIOS
#Set up the quantile probabilities
p <- c(0.90,0.95,0.975,0.99,0.995,0.999,0.9999,0.99999,0.999999)
#set up the alpha-type quantile probabilities for upper tail. This works in R but NOT in S-Plus
alpha <- c(0.10, 0.05, 0.025, 0.01, 0.005, 0.001, 0.0001, 0.00001, 0.000001)
sigma <- 0.2*10000/sqrt(250)
#Get a corresponding VaR vector which equals the quantile values calculated by inverting 
# normal probability fn and evaluating at quantile 
VaR.normal <- qnorm(p,sd=sigma)
#the following will not work in S-Plus since there is no lower-tail parameter
VaR.normal.ut <- qnorm(alpha, sd=sigma, lower.tail=FALSE)
ES.normal <- ESnorm(p,sd=sigma)
Ratio.normal <- ES.normal/VaR.normal
#Now look at VaR for student t with 4 degrees of freedom:
VaR.t4 <- qst(p,4,sigma=sigma,scale=TRUE)
ES.t4 <- ESst(p,4,sigma=sigma,scale=TRUE)
Ratio.t4 <- ES.t4/VaR.t4
#Display the comparisons in a matrix (table)
cbind(p,VaR.normal,VaR.t4,ES.normal,ES.t4,Ratio.normal,Ratio.t4)

