library(snow)
library(doSNOW)
library(foreach)
cores=25 #try 20-30?
reps=25
cl.tmp = makeCluster(rep("localhost",cores), type="SOCK")
registerDoSNOW(cl.tmp)
clusterSetupRNG(cl.tmp, seed=rep(18,6)) # originally rep(18,6)
out=foreach(i=1:100) %dopar% {
  #setwd("C:/Users/ba378/Dropbox/Robbie/Ben Analyses")
  setwd("D:/Workspace/emmetr/newest_october20_awd/Ben Analyses/")
  source("sim.Pack.R")
  source("mcmc.Pack.Final.R")
  source("mcmc.Pack.known.R") # want to look at both known and unknown probably
  # Rcpp::sourceCpp("siteup3D.cpp")
  library(SPIM)
  library(jagsUI)
  library(coda)
  Np=20  # CHANGE BY SCENARIO
  lambda.P=12.5 # up to 16-17 to get 250 individuals with Np=15 # CHANGE BY SCENARIO
  lam0.P=0.034
  sigma.P=4.50
  buff=3*sigma.P
  lambda.I=0.3 # Vary this by simulation # CHANGE BY SCENARIO
  scenfolder <- paste("D:/Workspace/emmetr/newest_october20_awd/Ben Analyses/sims/np", Np, "/lamp", lambda.P, "/lami", lambda.I, sep="")
  #set.seed(18) # delete later
  #data=sim.Pack(Np=Np,lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I,K=K,X=X,buff=buff)
  load(file=paste(scenfolder, "/rep", i, "/simdata.RData",sep="")) # if not tuning parms
  
  niter=30000 #just running 250 to set proppars # try 30k for simulations
  nburn=0
  nthin=1 # burnin and thin after the fact maybe?
  nind<-length(unique(data$packID[1:dim(data$y.I.obs)[1]])) # this is technically npack here - pack is ind
  K<-1
  ntraps<-dim(data$y.I.obs)[2]
  ymat_ind <- apply(data$y.I.obs, c(1,2), sum)
  ymat_group <- aggregate(ymat_ind, by=list(packID=data$packID[1:dim(data$y.I.obs)[1]]), sum)
  ymat_group <- as.matrix(ymat_group[,-1])
  M=60
  nz<-M-nind
  
  #create augmented array
  Yaug <- array(0, dim=c(M,ntraps))
  Yaug[1:nind,]<-ymat_group
  y<-Yaug
  
  #center the coordinates of the trap matrix
  X=data$X
  
  #set up the state-space
  
  Xl=min(X[,1]) - buff
  Xu=max(X[,1]) + buff
  Yl=min(X[,2]) - buff
  Yu=max(X[,2]) + buff
  areaX=(Xl-Xu)*(Yl-Yu)
  
  # Set up trapop
  trapop<-rep(31, ntraps)
  
  
  #get mean activity centers for observed bears; create initial values for remaining s
  Sin<-matrix(NA, ncol=2, nrow=M)
  
  for(i in 1:nind){
    active_coords <- (y[i,]>0)*X
    active_coords <- active_coords[active_coords[,2]!=0,]
    if(is.null(dim(active_coords))){
      active_coords <- array(active_coords, dim=c(1,2))
    }
    Sin[i,]<- colMeans(active_coords)
  }
  Sin[(nind+1):M,]<-cbind(runif(nz, Xl, Xu), runif(nz,Yl,Yu))
  
  dataIntoJAGS<-list(y=y,M=M, J=ntraps, Xl=Xl, Yl=Yl, Xu=Xu, Yu=Yu, X=X, area=areaX, trapop=trapop)
  parameters<-c('psi','lambda0','N', 'D', 'sigma')
  #set.seed(18)
  inits =  function() {list(z=c(rep(1,nind), rbinom(nz,1,0.5)),psi=runif(1), s=Sin, 
                            sigma=runif(1,0,15),lambda0=1) }
  scr_group <- tryCatch(jagsUI::jags(dataIntoJAGS, inits=inits, parameters.to.save = parameters, "D:/Workspace/emmetr/newest_october20_awd/SCR0_group.txt", n.chains=3, n.adapt=1000, n.iter=niter, n.burnin=nburn, n.thin = 1, parallel=FALSE), error=function(e){scr_group <- NA})
  
  
  # #comment these out when simulating
  # plot(mcmc(out$out))
  # 1-rejectionRate(mcmc(out$out))
  # min(1-rejectionRate(mcmc(out$sout[,,1])))
  
  return(scr_group) #originally out=out #save data set, too
}
stopCluster(cl.tmp)
#save(out,file="namethissimulationfilesomething")

# Process and save simulations HERE
# Make sure to note any that throw errors

# Get mcmc.list
# # TODO add IDout
# library(coda)
# bigout <- lapply(out, function(x){
#   #ifelse(any(is.na(list(x$chain1, x$chain2, x$chain3))), NA,mcmc.list(mcmc(x$chain1$out), mcmc(x$chain2$out), mcmc(x$chain3$out)))
#   if(any(is.na(list(x$chain1, x$chain2, x$chain3)))){
#     NA
#   } else {
#     mcmc.list(mcmc(x$chain1$out), mcmc(x$chain2$out), mcmc(x$chain3$out))
#   }
#   })
# bigidout <- lapply(out, function(x){
#   #ifelse(any(is.na(list(x$chain1, x$chain2, x$chain3))), NA,mcmc.list(mcmc(x$chain1$IDout[,1:dim(x$data$y.I.obs)[1]]), mcmc(x$chain2$IDout[,1:dim(x$data$y.I.obs)[1]]), mcmc(x$chain3$IDout[,1:dim(x$data$y.I.obs)[1]])))
#   if(any(is.na(list(x$chain1, x$chain2, x$chain3)))){
#     NA
#   } else {
#     mcmc.list(mcmc(x$chain1$IDout[,1:dim(x$data$y.I.obs)[1]]), mcmc(x$chain2$IDout[,1:dim(x$data$y.I.obs)[1]]), mcmc(x$chain3$IDout[,1:dim(x$data$y.I.obs)[1]]))
#   }
# })
setwd("D:/Workspace/emmetr/newest_october20_awd/Ben Analyses/sims/")
for(countone in c(20)){
  #dir.create(paste("np", countone,sep=""))
  for(counttwo in c(12.5)){
    #dir.create(paste("np", countone, "/lamp", counttwo, sep=""))
    for(countthree in c(0.3)){
      #dir.create(paste("np", countone, "/lamp", counttwo, "/lami", countthree, sep=""))
      for(countfour in 1:100){ # change as you add more reps
        tempout <- out[[countfour]] # right? -30? Think
        save(tempout,file=paste("np", countone, "/lamp", counttwo, "/lami", countthree, "/rep", countfour, "/basic_group_scr.RData", sep=""))
        #tempidout <- bigidout[[countfour]]
        #save(tempidout,file=paste("np", countone, "/lamp", counttwo, "/lami", countthree, "/rep", countfour, "/idsimresults.RData", sep=""))
        
      }
    }
  }
}

# Known-ID model

# Tuning
library(coda)
lapply(out, function(x){1-rejectionRate(mcmc(x$chain1$out))})
lapply(out, function(x){min(1-rejectionRate(mcmc(x$chain1$sout[,,1])))})

