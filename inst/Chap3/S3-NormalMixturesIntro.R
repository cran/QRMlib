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
#################################################


### Normal distribution: visualization, simulation, estimation
BiDensPlot(func=dmnorm,mu=c(0,0),Sigma=equicorr(2,-0.7))
ndata <- rmnorm(2000,rho=0.7,d=3)
fit.norm(ndata)

### t distribution: visualization, simulation, estimation
BiDensPlot(func=dmt,xpts=c(-4,4),ypts=c(-4,4),mu=c(0,0),Sigma=equicorr(2,-0.7),nu=4)
tdata <- rmt(2000,df=4,rho=0.7,d=3)
mod1 <- fit.mst(tdata)
mod0 <- fit.norm(tdata)

### (Multivariate generalized) Hyperbolic distribution: visualization with PERSPECTIVE or CONTOUR plots
par(mfrow=c(2,2))
ll <- c(-4,4)
BiDensPlot(func=dmghyp,xpts=ll,ypts=ll,mu=c(0,0),Sigma=equicorr(2,-0.7),lambda=1,chi=1,psi=1,gamma=c(0,0))
BiDensPlot(func=dmghyp,type="contour",xpts=ll,ypts=ll,mu=c(0,0),Sigma=equicorr(2,-0.7),lambda=1,chi=1,psi=1,gamma=c(0,0))
BiDensPlot(func=dmghyp,xpts=ll,ypts=ll,mu=c(0,0),Sigma=equicorr(2,-0.7),lambda=1,chi=1,psi=1,gamma=c(0.5,-0.5))
BiDensPlot(func=dmghyp,type="contour",xpts=ll,ypts=ll,mu=c(0,0),Sigma=equicorr(2,-0.7),lambda=1,chi=1,psi=1,gamma=c(0.5,-0.5))
par(mfrow=c(1,1))

### NIG distribution: visualization.  (Normal Inverse Gaussian with lambda = -0.5)
par(mfrow=c(2,2))
ll <- c(-2,2)
BiDensPlot(func=dmghyp,xpts=ll,ypts=ll,mu=c(0,0),Sigma=equicorr(2,-0.7),lambda=-0.5,chi=1,psi=1,gamma=c(0,0))
BiDensPlot(func=dmghyp,type="contour",xpts=ll,ypts=ll,mu=c(0,0),Sigma=equicorr(2,-0.7),lambda=-0.5,chi=1,psi=1,gamma=c(0,0))
BiDensPlot(func=dmghyp,xpts=ll,ypts=ll,mu=c(0,0),Sigma=equicorr(2,-0.7),lambda=-0.5,chi=1,psi=1,gamma=c(0.5,-0.5))
BiDensPlot(func=dmghyp,type="contour",xpts=ll,ypts=ll,mu=c(0,0),Sigma=equicorr(2,-0.7),lambda=-0.5,chi=1,psi=1,gamma=c(0.5,-0.5))
par(mfrow=c(1,1))

### GH (Generalized Hyperbolic) distributions: simulation
d <- 5;  #use 5 dimensions
n <- 2000;
P <- equicorr(d,0.7);
scale.factor <- (det(P))^(1/d);
Sigma <- P/scale.factor;
mu <- rep(0,d);
gamma <- (1:d)/10;
nu <- 4;
set.seed(13)
#See p. 80 in QRM for a description of these special distributions.  Note the .5d indicates we have 5 variables.
#multivariate whose marginals are 1-dimensional HYPERBOLIC distributions with lambda = 1 
data.hyp.5d <- rmghyp(n,lambda=1,chi=1,psi=1,Sigma=Sigma,mu=mu,gamma=gamma);
#NIG (normal inverse Gaussian) with lambda = -0.5 is  which is different from GIG (generalized inverse Gaussian)
data.nig.5d <- rmghyp(n,lambda=-0.5,chi=1,psi=1,Sigma=Sigma,mu=mu,gamma=gamma);
# multidimensional SKEWED (asymmetric) T distribution with lambda=-nu/1, chi=nu,psi=0:
data.t.5d <- rmghyp(n,lambda=(-nu/2),chi=nu,psi=0,Sigma=P,mu=mu,gamma=gamma);
#Full GENERALIZED HYPERBOLIC mean-variance mixture distribution:
data.gh.5d <- rmghyp(n,lambda=-1,chi=1,psi=1,d=5,Sigma=Sigma,mu=mu,gamma=gamma);
#VARIANCE-GAMMA distribution with lambda >0 and chi = 0. Also known as generalized Laplace or generalized Bessel.
data.vg.5d <- rmghyp(n,lambda=2,chi=0,psi=1,d=5,Sigma=Sigma,mu=mu,gamma=gamma);


###  GH distributions: fitting (NIG and hyp only)
mod1 <- fit.mNH(data.nig.5d,symmetric=FALSE,case="NIG");
mod2 <- fit.mNH(data.hyp.5d,symmetric=FALSE,case="hyp");

