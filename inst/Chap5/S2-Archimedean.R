# S-Plus script developed by Professor Alexander McNeil, A.J.McNeil@hw.ac.uk
# R-version adapted by Scott Ulman (scottulman@hotmail.com)
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
#################################################


# Simulating Archimedean copulas with six variables:
data.gumbelcopula <- rcopula.gumbel(1000,theta=2, d=6)
data.claytoncopula <- rcopula.clayton(1000, theta=1, d=6)
data.frankcopula <- rcopula.frank(1000, theta=5, d=6)
pairs(data.gumbelcopula)
hist(data.gumbelcopula[,3])
pairs(data.claytoncopula)
pairs(data.frankcopula)

# Tests of multivariate normality.  Generate values via quantile function
#for normal marginals 
data.metagumbel <- apply(data.gumbelcopula,2,qnorm)
#produce QQ-plot assuming normal distribution; use only 1 column of data from matrix:
qqnorm(data.metagumbel[,1])
jointnormalTest(data.metagumbel)

# Asymmetric Gumbel
data.AGumbel <- rcopula.AGumbel(10000,theta=4,alpha=c(0.95,0.7))
plot(data.AGumbel)
data.AGumbel <- rcopula.AGumbel(5000,theta=2,c(0.1,0.7,0.8,0.9))
pairs(data.AGumbel)
hist(data.AGumbel[,4])

# Gumbel 2 groups
data <- rcopula.Gumbel2Gp(n=3000,gpsizes=c(3,4),theta=c(2,3,5))
pairs(data)

# Gumbel Nested
data <- rcopula.GumbelNested(n=3000,theta=1:3)
pairs(data)

