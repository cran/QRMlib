# S-Plus script developed by Professor Alexander McNeil, mcneil@math.ethz.ch
# R-version adapted by Scott Ulman (scottulman@hotmail.com)
# This free script using QRMLib is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

######Load the QRMlib and spdata.raw data set##################
#QRMlib.pdf is a help file for the functions used by QRMlib.  It is available at
#...\Program Files\R\R-2.2.1\library\QRMlib\Docs
#If you have created the QRMBook workspace and .Rprofile  as described in QRMlib.pdf
#topics 'QRMBook-workspace' and 'profileLoadLibrary', then you may comment out the
#following line:
library(QRMlib);
#if you have previously opened the spdata.raw data set  AND saved 
#the workspace, you may comment out the following line:
data(spdata.raw);
#################################################


# Exploratory Analyses of Default Data
#Show the raw credit data
spdata.raw;
#Attach the data set so you can use the data variable names without prepending the dataset name.
attach(spdata.raw);

# Calculate default rates.  Bobligors is the number of issuers with a B credit rating.
#Bdefaults is a vector containing the number of defaults in the group with B credit rating.
BdefaultRate <- Bdefaults/Bobligors;
BBdefaultRate <- BBdefaults/BBobligors;
BBBdefaultRate <- BBBdefaults/BBBobligors;
AdefaultRate <- Adefaults/Aobligors;
CCCdefaultRate <- CCCdefaults/CCCobligors;

# Plot default rates

plot(year,CCCdefaultRate,xlab="Year",ylab="Rate",type="l");
lines(year,BdefaultRate,col=2);
lines(year,BBdefaultRate,col=3);
lines(year,BBBdefaultRate,col=4);
lines(year,AdefaultRate,col=5);

