#### Example script to simulate data, fit a model, and analyze results ####
# R version 4.0.1, Rtools 40
library(SPIM)
library(coda) # version 0.19-3
library(MCMCglmm) # version 2.29

#### Set working directory here ####
setwd("C:/Users/Robbie/Desktop/dissertation/results_and_writing/dissertation/Ecology_SuppInfo_Code_Emmet_et_al_Revised/")
source("mcmc.Pack.known.R")
source("sim.Pack.R")

#### Simulate data ####
set.seed(18)
Np=20 # Number of groups
zeta=12.5 # Average group size parameter for zero-truncated Poisson distribution 
lam0=0.034 # baseline group encounter rate
sigma=4.50 # group encounter scale parameter
theta=0.3 # individual detection rate
K=31 # number of occasions
spacing=2/3
X=expand.grid(1:11,1:12)*spacing*sigma # camera grid
buff=sigma*3 # buffer around camera grid
data=sim.Pack(Np=Np,zeta=zeta,lam0=lam0,sigma=sigma,theta=theta,K=K,X=X,buff=buff)

#### Run model ####
niter=500 # number of iterations
nburn=0 # amount of burnin
nthin=1 # thin rate
Z=60 # group superpopulation size
mincluster=50 # the minimum number of total individuals (observed and augmented) that should be in each group to ensure proper group sizes can be proposed
W=3600 # individual superpopulation size - set to at least Z*mincluster
# Initial values
inits=list(zeta=zeta,lam0=lam0,sigma=sigma,theta=theta,psi1=0.5,psi=0.5)
# Proposal distribution parameters (standard deviations of normal distributions centered at the current value of each)
proppars=list(zeta=1.3,lam0=0.004,sigma=0.25,theta=0.04,s=2.0,ncluster=1,MoveOneProb=0.5)




storeLatent=list(ID=FALSE,s=FALSE,z=FALSE)#set s to true to set proposal parameters
plotLatent=FALSE #turn off for simulations
J=nrow(X) # Number of detectors
trapop=matrix(1,nrow=J,ncol=K) # Trap operation matrix - detectors by occasions
yprime=apply(data$y.I.obs,c(1,2),sum) # Capture histories summed over occasions
yprime=array(yprime,dim=c(dim(yprime),1))
# Since our models used one occasion, we use the number of original simulated occasions (31) in a trap operation matrix
trapopprime=rowSums(trapop) # The trap operation file is now Jx1-the number of occasions in which each detector was operational
trapopprime=matrix(trapopprime,nrow=length(trapopprime),ncol=1)
data$trapop=trapopprime
data$y.I.obs=yprime
#if you want to initTrue, this adds group encounters ("SiteVisits" in the output) to the data
yprimefull=apply(data$y.I,c(1,2),sum)
yprimefull=array(yprimefull,dim=c(dim(yprimefull),1))
yPprime=apply(data$v,c(1,2),sum) #yPprime is the sum of group encounters for each detector and group over all occasions
yPprime=array(yPprime,dim=c(dim(yPprime),1))
data$y.I=yprimefull
data$v=yPprime

chain1=mcmc.Pack.known(data,niter=niter,nthin=nthin,nburn=nburn,W=W,Z=Z,inits=inits,proppars=proppars,
                                mincluster=mincluster,storeLatent=storeLatent,
                                plotLatent=plotLatent)

# Basic summary statistics, acceptance rate, etc.
# Npack is the realized number of groups, Nind is the realized number of individuals, and SiteVisits is the total sum of all group encounters
chain1mcmc <- mcmc(chain1$out)
posterior.mode(chain1mcmc)
HPDinterval(chain1mcmc)
1-rejectionRate(chain1mcmc)
