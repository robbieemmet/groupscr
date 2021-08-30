# This 
#match 132 traps from real data, but in regular array 11 x 12
#Simulate from estimated parameters.
#Assume area is fixed, match this regular grid to 6931 km^2 from data set
spacing=2/3 #trap spacing in sigma units
nrep=10
store=rep(NA,nrep)

# Read in stuff
setwd("D:/Workspace/emmetr/newest_october20_awd/Ben Analyses/")
source("sim.Pack.R")
#source("mcmc.Pack.Final.R")
#source("mcmc.Pack.known.R") # want to look at both known and unknown probably
# Rcpp::sourceCpp("siteup3D.cpp")
library(SPIM)

# Set up scenario vectors
packabund <- c(10, 20)
packsize <- c(12.5, 5) # maybe add in later?
lamI <- c(0.3)
nsimdata <- 100

# Set up folder structure
setwd("D:/Workspace/emmetr/newest_october20_awd/Ben Analyses/sims/")
for(countone in c(10,20)){
  dir.create(paste("np", countone,sep=""))
  for(counttwo in c(5,12.5)){
    dir.create(paste("np", countone, "/lamp", counttwo, sep=""))
    for(countthree in c(0.3)){
      dir.create(paste("np", countone, "/lamp", counttwo, "/lami", countthree, sep=""))
      for(countfour in 1:100){ # change as you add more reps
        dir.create(paste("np", countone, "/lamp", counttwo, "/lami", countthree, "/rep", countfour, sep=""))
        
      }
    }
  }
}

# Simulate data
set.seed(18) # if you need to add more, change seed for new data
for(countone in c(10,20)){
  Np=countone # vary this by simulation # up to 20 to get 250 individuals with lambda.P=12.5 - originally 15 to get 2014 density
  for(counttwo in c(5,12.5)){
    lambda.P=counttwo
  for(countthree in c(0.3)){
    #lambda.P=12.5 # up to 16-17 to get 250 individuals with Np=15
    lam0.P=0.034
    sigma.P=4.50
    lambda.I=countthree # Vary this by simulation
    K=31 # why this?
    spacing=2/3 # originally 1.2
    X=expand.grid(1:11,1:12)*spacing*sigma.P 
    buff=sigma.P*3
    #set.seed(18) # delete later
    for(countfour in 1:100){
    data=sim.Pack(Np=Np,lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I,K=K,X=X,buff=buff)
    save(data, file=paste("np", countone, "/lamp", counttwo, "/lami", countthree, "/rep", countfour, "/simdata.RData", sep=""))
    } # countfour
  } # countthree
  } #
} # add progress bar

# # Intermediate - look at simulated data pack recaps, ind recaps, npack, nind, etc.
data_summaries <- data.frame()
for(countone in c(10,20)){
  Np=countone # vary this by simulation # up to 20 to get 250 individuals with lambda.P=12.5 - originally 15 to get 2014 density
  for(counttwo in c(5,12.5)){
    lambda.P=counttwo
    for(countthree in c(0.3)){
      #lambda.P=12.5 # up to 16-17 to get 250 individuals with Np=15
      #lam0.P=0.034
      #sigma.P=4.50
      #lambda.I=countthree # Vary this by simulation
      #K=31 # why this?
      #spacing=2/3 # originally 1.2
      #X=expand.grid(1:11,1:12)*spacing*sigma.P
      #buff=sigma.P*3
      #set.seed(18) # delete later
      for(countfour in 1:100){
        #data=sim.Pack(Np=Np,lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I,K=K,X=X,buff=buff)
        load(file=paste("np", countone, "/lamp", counttwo, "/lami", countthree, "/rep", countfour, "/simdata.RData", sep=""))
# Possible correlates of model performance - nind, some recapture statistics, group overlap, individual-group recapture ratio, etc.
nind <- nrow(data$y.I.obs)
npack <- length(unique(data$packID[1:dim(data$y.I.obs)[1]]))

# Figure out some measure of pack ID success
packID.true <- data$packID[1:dim(data$y.I.obs)[1]]

# Co-detection proportion and group overlap calculations
total_det_sites <- sum(apply(data$y.I.obs,2,sum)>0)

pack_by_trap <- aggregate(apply(data$y.I.obs,c(1,2),sum), by=list(packID=packID.true), sum)
pack_by_trap <- pack_by_trap[,-1]

# Pack recaptures and individual recaptures
inds_captured <- apply(data$y.I.obs, c(1,2), sum)>0
ind_recaptures <- mean(apply(inds_captured, 1, sum))

pack_recaptures <- mean(apply(pack_by_trap>0, 1, sum))

pack_site_visits <- sum(data$y.P)/dim(data$y.P)[1]
data_summaries <- rbind(data_summaries, c(nind, npack, ind_recaptures, pack_recaptures, pack_site_visits))

} # countfour
} # countthree
} #
} # add progress bar

# Brief analysis
names(data_summaries) <- c("nind", "npack", "ind_recaptures","pack_recaptures", "pack_site_visits")

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
  Np=20  # CHANGE BY SCENARIO
  lambda.P=12.5 # up to 16-17 to get 250 individuals with Np=15 # CHANGE BY SCENARIO
  lam0.P=0.034
  sigma.P=4.50
  lambda.I=0.3 # Vary this by simulation # CHANGE BY SCENARIO
  scenfolder <- paste("D:/Workspace/emmetr/newest_october20_awd/Ben Analyses/sims/np", Np, "/lamp", lambda.P, "/lami", lambda.I, sep="")
  K=31 # why this?
  spacing=2/3
  X=expand.grid(1:11,1:12)*spacing*sigma.P
  buff=sigma.P*3
  #set.seed(18) # delete later
  #data=sim.Pack(Np=Np,lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I,K=K,X=X,buff=buff)
  load(file=paste(scenfolder, "/rep", i, "/simdata.RData",sep="")) # if not tuning parms
  # par(mfrow=c(1,1),ask=FALSE)
  # xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  # ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  # plot(NA,pch=16,xlim=xlim,ylim=ylim,xlab="X",ylab="Y")
  # visiti=which(rowSums(data$y.P)>0)
  # y.P2D=apply(data$y.P,c(1,2),sum)
  # for(i in visiti){
  #   trapcap=which(y.P2D[i,]>0)
  #   for(j in trapcap){
  #     lines(x=c(data$s[i,1],X[j,1]),y=c(data$s[i,2],X[j,2]),lty=3)
  #   }
  # }
  # capi=which(rowSums(data$y.I)>0)
  # y.I2D=apply(data$y.I,c(1,2),sum)
  # for(i in capi){
  #   trapcap=which(y.I2D[i,]>0)
  #   for(j in trapcap){
  #     ID=data$packID[i]
  #     lines(x=c(data$s[ID,1],X[j,1]),y=c(data$s[ID,2],X[j,2]),lty=1,col="#999999",lwd=2)
  #   }
  # }
  # points(data$s,cex=2,col="#56B4E9",pch=16)
  # text(data$s[,1],data$s[,2],label=as.character(data$NperPack),cex=0.75)
  # points(X,pch=4)
  
  niter=30000 #just running 250 to set proppars # try 30k for simulations
  nburn=0
  nthin=1
  M2=60 #Not sure how high this will need to be. Originally 40 # check if it needs to be higher
  if(Np==10 & lambda.P==5){
  M1=2200
  mincluster=30 # originally 30
  
  inits=list(lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I,psi1=0.5,psi2=0.5)
  proppars=list(lambda.P=0.9,lam0.P=0.004,sigma.P=0.25,lambda.I=0.06,s=1.8,ncluster=1,MoveOneProb=0.5)
  }
  
  if(Np==10 & lambda.P==12.5){
    M1=3600
    mincluster=50 # originally 30
    
    inits=list(lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I,psi1=0.5,psi2=0.5)
    proppars=list(lambda.P=1.3,lam0.P=0.004,sigma.P=0.25,lambda.I=0.06,s=1.5,ncluster=1,MoveOneProb=0.5)
  }
  
  if(Np==20 & lambda.P==5){
    M1=2200
    mincluster=30 # originally 30
    
    inits=list(lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I,psi1=0.5,psi2=0.5)
    proppars=list(lambda.P=0.9,lam0.P=0.004,sigma.P=0.25,lambda.I=0.06,s=2.0,ncluster=1,MoveOneProb=0.5)
  }
  
  if(Np==20 & lambda.P==12.5){
    M1=3600
    mincluster=50 # originally 30
    
    inits=list(lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I,psi1=0.5,psi2=0.5)
    proppars=list(lambda.P=1.3,lam0.P=0.004,sigma.P=0.25,lambda.I=0.04,s=2.0,ncluster=1,MoveOneProb=0.5)
  }
  
  
  
  #initTrue=TRUE #want to initTrue when setting proppars
  storeLatent=list(ID=FALSE,s=FALSE,z=FALSE)#set s to true to set proppars
  plotLatent=FALSE #turn off for simulations
  J=nrow(X)
  trapop=matrix(1,nrow=J,ncol=K)
  yprime=apply(data$y.I.obs,c(1,2),sum)
  yprime=array(yprime,dim=c(dim(yprime),1))
  trapopprime=rowSums(trapop)
  trapopprime=matrix(trapopprime,nrow=length(trapopprime),ncol=1)
  data$trapop=trapopprime
  data$y.I.obs=yprime
  #if you want to initTrue
  yprimefull=apply(data$y.I,c(1,2),sum)
  yprimefull=array(yprimefull,dim=c(dim(yprimefull),1))
  yPprime=apply(data$y.P,c(1,2),sum)
  yPprime=array(yPprime,dim=c(dim(yPprime),1))
  data$y.I=yprimefull
  data$y.P=yPprime
  
  # add pack IDs of detected individuals to data for known packs
  #data$packID.obs <- as.numeric(as.factor(data$packID[unique(which(data$y.I>0,arr.ind = T)[,1])])) # TODO is this OK?

  
  
  a=Sys.time()
  #set.seed(18) #make repeatable
  # add tryCatch to this in a sensible way
  chain1=tryCatch(mcmc.Pack.known(data,niter=niter,nthin=nthin,nburn=nburn,M1=M1,M2=M2,inits=inits,proppars=proppars,
                mincluster=mincluster,storeLatent=storeLatent,
                plotLatent=plotLatent), error=function(e){chain1 <- NA})
  chain2=tryCatch(mcmc.Pack.known(data,niter=niter,nthin=nthin,nburn=nburn,M1=M1,M2=M2,inits=inits,proppars=proppars,
                            mincluster=mincluster,storeLatent=storeLatent,
                            plotLatent=plotLatent), error=function(e){chain2 <- NA})
  chain3=tryCatch(mcmc.Pack.known(data,niter=niter,nthin=nthin,nburn=nburn,M1=M1,M2=M2,inits=inits,proppars=proppars,
                            mincluster=mincluster,storeLatent=storeLatent,
                            plotLatent=plotLatent), error=function(e){chain3 <- NA})
  b=Sys.time()
  # b-a
  
  # #comment these out when simulating
  # plot(mcmc(out$out))
  # 1-rejectionRate(mcmc(out$out))
  # min(1-rejectionRate(mcmc(out$sout[,,1])))
  
  return(list(chain1=chain1,chain2=chain2,chain3=chain3,time=b-a,data=data)) #originally out=out #save data set, too
}
stopCluster(cl.tmp)
#save(out,file="namethissimulationfilesomething")

# Process and save simulations HERE
# Make sure to note any that throw errors

# Get mcmc.list
# TODO add IDout
library(coda)
bigout <- lapply(out, function(x){
  #ifelse(any(is.na(list(x$chain1, x$chain2, x$chain3))), NA,mcmc.list(mcmc(x$chain1$out), mcmc(x$chain2$out), mcmc(x$chain3$out)))
  if(any(is.na(list(x$chain1, x$chain2, x$chain3)))){
    NA
  } else {
    mcmc.list(mcmc(x$chain1$out), mcmc(x$chain2$out), mcmc(x$chain3$out))
  }
  })
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
        tempout <- bigout[[countfour]] # right? -30? Think
        save(tempout,file=paste("np", countone, "/lamp", counttwo, "/lami", countthree, "/rep", countfour, "/knownidsimresults.RData", sep=""))
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

