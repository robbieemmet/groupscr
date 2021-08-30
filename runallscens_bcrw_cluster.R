# This 
#match 132 traps from real data, but in regular array 11 x 12
#Simulate from estimated parameters.
#Assume area is fixed, match this regular grid to 6931 km^2 from data set
# This time we simulate biased correlated random walk movement
spacing=1.2 #trap spacing in sigma units
nrep=10
store=rep(NA,nrep)

# Read in stuff
setwd("D:/Workspace/emmetr/newest_october20_awd/Ben Analyses/")
source("sim.Pack.R")
#source("mcmc.Pack.Final.R")
#source("mcmc.Pack.known.R") # want to look at both known and unknown probably
# Rcpp::sourceCpp("siteup3D.cpp")
library(SPIM)
library(momentuHMM) # dev version as of July 2019 on RLE's computer
library(raster)
library(ctmcmove)
library(VGAM)

# Set up scenario vectors
packabund <- c(10, 20)
packsize <- c(12.5, 5) # maybe add in later?
lamI <- c(0.3)
nsimdata <- 100

# Set up folder structure
setwd("D:/Workspace/emmetr/newest_october20_awd/Ben Analyses/sims/")
for(countone in c(10,20)){
  #dir.create(paste("np", countone,sep=""))
  for(counttwo in c(5,12.5)){
    #dir.create(paste("np", countone, "/lamp", counttwo, sep=""))
    for(countthree in c(0.3)){
      for(countattract in c(-1, -0.1)){
      dir.create(paste("np", countone, "/lamp", counttwo, "/lami", countthree, "/attr", abs(countattract), sep=""))
      for(countfour in 1:100){ # change as you add more reps
        dir.create(paste("np", countone, "/lamp", counttwo, "/lami", countthree, "/attr", abs(countattract), "/rep", countfour, sep=""))
        #print("hi")
        # To potentially add after countthree loop - attraction/repulsion scenarios
      }
      }
    }
  }
}

# Simulate data
set.seed(18)
library(snow)
library(doSNOW)
library(foreach)
for(countone in c(10,20)){
  Np=countone # vary this by simulation # up to 20 to get 250 individuals with lambda.P=12.5 - originally 15 to get 2014 density
  for(counttwo in c(5,12.5)){
    lambda.P=counttwo
  for(countthree in c(0.3)){
    #lambda.P=12.5 # up to 16-17 to get 250 individuals with Np=15
    lam0.P=1
    sigma.P=4.50
    lambda.I=countthree # Vary this by simulation
    K=31 # why this?
    spacing=2/3
    X=expand.grid(1:11,1:12)*spacing*sigma.P 
    J <- nrow(X)
    buff=sigma.P*3
    #set.seed(18) # delete later
    # Add loop for repulsion parameters?
    for(countattract in c(-1, -0.1)){
    cores=25 #try 20-30?
    reps=25
    cl.tmp.hmm = makeCluster(rep("localhost",cores), type="SOCK")
    registerDoSNOW(cl.tmp.hmm)
    clusterSetupRNG(cl.tmp.hmm, seed=rep(18,6))
    outhmm <- foreach(countfour=1:100) %dopar% { # foreach this fella
    
      #### Groups with variable/dynamic group size
      ## 1/6/2020


      ## ----------------- SCR with explicit movement process ----------------
      ## Group Dynamic Movement - variable group size
      library(SPIM)
      library(momentuHMM) # dev version as of July 2019 on RLE's computer
      library(raster)
      library(ctmcmove)
      library(VGAM)
      
        skip_to_next <- FALSE
        #N <- 80     ## Abundance - originally 250
        #M <- N+100     ## Superpopulation size - see line 256 (approx) instead
        Ngroups <- Np  ## Number groups - group 25 got lost somewhere down below - originally 25
        lambda0 <- 2  ## Local baseline encounter rate. Lots of photos if near camera
        sigma <- 1 ## SD of stationary distribution
        rho <- 0.9    ## Correlation parameter
        sigmaDet <- 0.45 ## Scale parameter of local detection function - originally 0.005
        ylim.S <- range(X[,2])+c(-buff, buff)
        xlim.S <- range(X[,1])+c(-buff, buff)
        ## Activity centers
        s <- cbind(runif(Ngroups, xlim.S[1], xlim.S[2]),
                   runif(Ngroups, ylim.S[1], ylim.S[2]))

        # create potential surface to keep animals within state space
        
        ncols <- 11 # nrows in secr grid
        nrows <- 12 # ncols in secr grid
        potSurface <- raster(ncol=11*5,nrow=12*5,xmn=xlim.S[1]-diff(xlim.S)*2,xmx=xlim.S[2]+diff(xlim.S)*2,ymn=ylim.S[1]-diff(ylim.S)*2,ymx=ylim.S[2]+diff(ylim.S)*2, vals=rep(0,(11*5)*12*5))
        crs(potSurface) <- CRS("+proj=utm +zone=16N +datum=WGS84 +units=km +no_defs")
        # Set values to 0 for middle cells (which are?), otherwise 1
        potSurface<-setValues(potSurface,values=1-c(rep(rep(0,ncols*5),nrows*2),rep(c(rep(0,ncols*2),rep(1,ncols),rep(0,ncols*2)),nrows),rep(rep(0,ncols*5),nrows*2)))
        potSurface[potSurface>0] <- NA
        potSurface <- raster::distance(potSurface)*1000
        potSurfaceRast <- ctmcmove::rast.grad(potSurface)[c("rast.grad.x","rast.grad.y")]
        
        # Should I do absolute value to make sure gradient is in right direction?
        #potSurfaceRast$rast.grad.x <- abs(potSurfaceRast$rast.grad.x)
        #potSurfaceRast$rast.grad.y <- abs(potSurfaceRast$rast.grad.y)
        
        ## Locations at each occasion (hour)
        
        #simulate a group centroid as a single-state biased correlated random walk (N=Ngroups)
        
        mu0 <- list()
        for(i in 1:Ngroups){
          mu0[[i]] <- s[i,]
        }
        
        centroidData <- tryCatch(simData(nbAnimals=Ngroups,nbStates=1,
                                obsPerAnimal=K,
                                dist=list(ctr="rw_mvnorm2"),
                                DM=list(ctr=list(mean.x=~crw(ctr.x_tm1)+rast.grad.x, #t minus 1
                                                 mean.y=~crw(ctr.y_tm1)+rast.grad.y, #t minus 1
                                                 sigma.x=~1,
                                                 sigma.xy=~1,
                                                 sigma.y=~1)),
                                Par=list(ctr=c(1,
                                               0.3, #amount of correlation
                                               -0.5,
                                               1,
                                               0.3, #amount of correlation
                                               -0.5,
                                               log(sigma^2),
                                               0,
                                               log(sigma^2))),
                                mvnCoords="ctr",
                                initialPosition = mu0,
                                spatialCovs = potSurfaceRast), error=function(e){skip_to_next <- TRUE})
        #plot(centroidData,dataNames=c("ctr.x","ctr.y"),compact=TRUE)
        if(skip_to_next){next}
        
        # create centroid data frame
        cD <- data.frame(x = centroidData$ctr.x, y = centroidData$ctr.y)
        
        #simulate the movement of a group of K individuals as a random walk (relative to the centroid)
        # This is the part where you add in random group sizes - individual abundance is not fixed
        #meanGroupSize <- N/Ngroups
        #groupSizes<-sample((meanGroupSize-2):(meanGroupSize+3), size=Ngroups, replace=T) #random sample around mean group size (upper bound higher to ensure enough N are generated)
        
        groupSizes <- rzapois(Ngroups, counttwo)
        
        #population initial positions for individuals in each group
        mu0_2 <- list()
        for (i in 1:Ngroups) {
          individual.starts <- list()
          for(j in 1:groupSizes[i]){
            individual.starts[[j]] <- c(rnorm(1, mu0[[i]][1], sigma), rnorm(1, mu0[[i]][2], sigma)) #initialize from initial position of each centroid track
          }
          mu0_2[[i]] <- individual.starts
        }
        
        groupData = data.frame()
        for (i in 1:Ngroups){
          grpDat <- tryCatch(simData(nbAnimals=groupSizes[i],nbStates=1,
                            obsPerAnimal=K,
                            dist=list(mu="rw_mvnorm2"),
                            covs=centroidData[centroidData$ID==i,c("ctr.x","ctr.y")],
                            DM=list(mu=list(mean.x=~I(mu.x_tm1 - ctr.x) + crw(mu.x_tm1) + rast.grad.x, #t minus 1
                                            mean.y=~I(mu.y_tm1 - ctr.y) + crw(mu.x_tm1) + rast.grad.y, #t minus 1
                                            sigma.x=~1,
                                            sigma.xy=~1,
                                            sigma.y=~1)),
                            Par=list(mu=c(1,
                                          countattract, #beta for centroid attraction (x-direction)
                                          0.3, #beta for correlation coefficient (x-direction)
                                          -0.5,
                                          1,
                                          countattract, #beta for centroid attraction (y-direction)
                                          0.3, #beta for correlation coefficient (y-direction)
                                          -0.5,
                                          log(sigma^2),
                                          0,
                                          log(sigma^2))),
                            mvnCoords="mu",
                            centroids = list(centroid = cD),
                            initialPosition = mu0_2[[i]],
                            spatialCovs = potSurfaceRast), error=function(e){skip_to_next <- TRUE})
          if(!skip_to_next){
          grpDat$Group <- i
          grpDat$id <- grpDat$ID
          grpDat$ID <- paste(grpDat$Group, grpDat$id, sep="_")
          groupData <- rbind(groupData, grpDat)
          } else {
            break
          }
        }
        if(skip_to_next){next}
        #plot(groupData,dataNames=c("mu.x","mu.y"),compact=TRUE)
        
        # #Check group movements against their centroid path
        # ggplot() + geom_path(data=groupData, aes(x=mu.x, y=mu.y, color=id, group=ID)) + 
        #    theme_classic() + geom_path(data=centroidData, aes(x=ctr.x, y=ctr.y, group=ID), size=2)
        grpID=unique(groupData$ID)
        groupData$indID=pmatch(groupData$ID, grpID, duplicates.ok = TRUE)
        
        # Ensure N=250 - scratch that
        #groupData <- groupData[1:(N*K),] # Will this ensure that Ngroups is the same?
        #length(unique(groupData$ID))
        
        # #Check group movements against their centroid path
        #ggplot() + geom_path(data=groupData, aes(x=mu.x, y=mu.y, color=id, group=ID)) +
        # theme_classic() + geom_path(data=centroidData, aes(x=ctr.x, y=ctr.y, group=ID), size=1.5) +
        #geom_point(data=traps, aes(x=x/1000,y=y/1000),size=1.5)
        
        #ggplot() + geom_path(data=groupData, aes(x=mu.x, y=mu.y, color=factor(Group), group=ID)) +
         #  theme_classic() + geom_path(data=centroidData, aes(x=ctr.x, y=ctr.y, group=ID), size=1.5)+
        # geom_point(data=X, x=Var1, y=Var2)
        # ^ Add camera locations?
        M <- 2200 # large number...
        uall<-array(NA,dim=c(M,K,2))
        d <- lambda <- yall <- array(NA,dim=c(M,J,K))
        Nind <- length(unique(groupData$indID))
        ## Capture-recapture process
        for(i in 1:Nind) {
          tmp <- subset(groupData,indID==i) #groupData instead of simDat?
          for(k in 1:K) { 
            uall[i,k,1] <- tmp[k,"mu.x"]
            uall[i,k,2]  <- tmp[k,"mu.y"]
            d[i,,k] <- sqrt((uall[i,k,1] - X[,1])^2 +
                              (uall[i,k,2] - X[,2])^2)
            lambda[i,,k] <- lambda0*exp(-d[i,,k]^2 /
                                          (2*sigmaDet^2))
            yall[i,,k] <- rpois(J, lambda[i,,k])
          }
        }
        
        
        detected <- apply(yall>0, 1, any, na.rm=TRUE)
        y <- yall[detected,,,drop=FALSE]
      y.I.obs <- apply(y, c(1,2), sum)
      # To comment out: look at some summary stats
      # Observed group IDs
      packID <- as.numeric(apply(table(groupData$Group, groupData$indID), 2, function(x){which(x>0)}))[detected]
      packID <- as.numeric(as.factor(packID))
      pack_by_trap <- aggregate(y.I.obs, by=list(packID=packID), sum)[,-1]  
        
      # Individual and pack recaptures
      inds_captured <- y.I.obs>0
      ind_recaptures <- mean(apply(inds_captured, 1, sum))
      
      pack_recaptures <- mean(apply(pack_by_trap>0, 1, sum))
      
      # Data to save/create (include ind_recaptures and pack_recaptures ideally):
      # X, K, packID.obs, trapopprime (JxK I guess?), yprime (IxJx1 array), buff
      data <- list()
      data$X <- X
      data$K <- K
      data$packID.obs <- packID
      data$y.I.obs <- array(y.I.obs, dim=c(dim(y.I.obs), 1)) # maybe save yourself some effort below? Do I want to save any movement tracks?
      data$trapop <- NULL # assuming perfect trap operation
      data$buff <- buff
      data$Nind <- Nind
    save(data, file=paste("np", countone, "/lamp", counttwo, "/lami", countthree, "/attr", abs(countattract), "/rep", countfour, "/hmmdata.RData", sep=""))
    return(list(data=data))
    } # countfour
    stopCluster(cl.tmp.hmm)
    

    
    } # attraction parameter
  } # countthree
  } #
} # add progress bar

# Look at the basics for each scenario
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
      for(countattract in c(-1, -0.1)){
        scenfolder <- paste("D:/Workspace/emmetr/newest_october20_awd/Ben Analyses/sims/np", Np, "/lamp", lambda.P, "/lami", lambda.I, "/attr", abs(countattract), sep="")
        
      for(countfour in 1:100){
        
        #K=31 # why this?
        #spacing=2/3
        #X=expand.grid(1:11,1:12)*spacing*1
        #J=nrow(X)
        #buff=1*3
        
        if(!file.exists(paste(scenfolder,"/rep", countfour, "/hmmdata.RData", sep=""))){
          next 
        }
        #set.seed(18) # delete later
        #data=sim.Pack(Np=Np,lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I,K=K,X=X,buff=buff)
        load(file=paste(scenfolder, "/rep", countfour, "/hmmdata.RData",sep=""))
        # Possible correlates of model performance - nind, some recapture statistics, group overlap, individual-group recapture ratio, etc.
        nind <- dim(data$y.I.obs)[1]
        npack <- length(unique(data$packID.obs))
        
        # Figure out some measure of pack ID success
        packID.true <- data$packID.obs
        
        # Co-detection proportion and group overlap calculations
        total_det_sites <- sum(apply(data$y.I.obs,2,sum)>0)
        
        pack_by_trap <- aggregate(apply(data$y.I.obs,c(1,2),sum), by=list(packID=packID.true), sum)
        pack_by_trap <- pack_by_trap[,-1]
        
        # Pack recaptures and individual recaptures
        inds_captured <- apply(data$y.I.obs, c(1,2), sum)>0
        ind_recaptures <- mean(apply(inds_captured, 1, sum))
        
        pack_recaptures <- mean(apply(pack_by_trap>0, 1, sum))
        
        #pack_site_visits <- sum(data$y.P)/dim(data$y.P)[1]
        data_summaries <- rbind(data_summaries, c(countone, counttwo, countattract, countfour, nind, npack, ind_recaptures, pack_recaptures))
        
      } # countfour
      } # countattract
    } # countthree
  } #
} # add progress bar

# Brief analysis
names(data_summaries) <- c("Npack", "lambda.P", "attract", "rep", "nind", "npack", "ind_recaptures","pack_recaptures")


## Check which data sets didn't get simulated
for(countone in c(10, 20)){
  for(counttwo in c(5, 12.5)){
    for(countthree in c(0.3)){
      for(countattract in c(-1, -0.1)){
        for(countfour in 1:100){
          if(!file.exists(paste("np", countone, "/lamp", counttwo, "/lami", countthree, "/attr", abs(countattract), "/rep", countfour, "/hmmdata.RData", sep=""))){
            print(paste("np", countone, "/lamp", counttwo, "/lami", countthree, "/attr", abs(countattract), "/rep", countfour, "/hmmdata.RData", sep=""))
          }
        }
      }
    }
  }
}

# ^ Find some way to resimulate those data I guess?



library(snow)
library(doSNOW)
library(foreach)
cores=34 #try 20-30?
reps=34
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
  lam0.P=0.3 # initial value
  sigma.P=3 # initial value
  lambda.I=0.3 # Vary this by simulation # CHANGE BY SCENARIO
  countattract=-0.1
  scenfolder <- paste("D:/Workspace/emmetr/newest_october20_awd/Ben Analyses/sims/np", Np, "/lamp", lambda.P, "/lami", lambda.I, "/attr", abs(countattract), sep="")
  K=31 # why this?
  spacing=2/3
  X=expand.grid(1:11,1:12)*spacing*1
  J=nrow(X)
  #buff=1*3
  
  if(!file.exists(paste(scenfolder,"/rep", i, "/hmmdata.RData", sep=""))){
    return(NA)  
    }
  #set.seed(18) # delete later
  #data=sim.Pack(Np=Np,lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I,K=K,X=X,buff=buff)
  load(file=paste(scenfolder, "/rep", i, "/hmmdata.RData",sep="")) # if not tuning parms
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
  data$trapop <- matrix(K, nrow=J, ncol=1)
  niter=30000 #just running 250 to set proppars # try 30k for simulations
  nburn=0
  nthin=1
  M2=60 #Not sure how high this will need to be. Originally 40 # check if it needs to be higher
  if(Np==10 & lambda.P==5){
    M1=2200
    mincluster=30 # originally 30
    
    inits=list(lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I,psi1=0.5,psi2=0.5)
    proppars=list(lambda.P=0.9,lam0.P=0.02,sigma.P=0.2,lambda.I=0.04,s=2.5,ncluster=1,MoveOneProb=0.5)
  }
  
  if(Np==10 & lambda.P==12.5){
    M1=3600
    mincluster=50 # originally 30
    
    inits=list(lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I,psi1=0.5,psi2=0.5)
    proppars=list(lambda.P=1.3,lam0.P=0.02,sigma.P=0.25,lambda.I=0.04,s=2.5,ncluster=1,MoveOneProb=0.5)
  }
  
  if(Np==20 & lambda.P==5){
    M1=2200
    mincluster=30 # originally 30
    
    inits=list(lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I,psi1=0.5,psi2=0.5)
    proppars=list(lambda.P=0.9,lam0.P=0.02,sigma.P=0.25,lambda.I=0.04,s=2.5,ncluster=1,MoveOneProb=0.5)
  }
  
  if(Np==20 & lambda.P==12.5){
    M1=3600
    mincluster=50 # originally 30
    
    inits=list(lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I,psi1=0.5,psi2=0.5)
    proppars=list(lambda.P=1.3,lam0.P=0.02,sigma.P=0.25,lambda.I=0.04,s=2.5,ncluster=1,MoveOneProb=0.5)
  }
  
  
  
  #initTrue=TRUE #want to initTrue when setting proppars
  storeLatent=list(ID=FALSE,s=FALSE,z=FALSE)#set s to true to set proppars
  plotLatent=FALSE #turn off for simulations
  #J=nrow(X)
  #trapop=matrix(1,nrow=J,ncol=K)
  #yprime=apply(data$y.I.obs,c(1,2),sum)
  #yprime=array(yprime,dim=c(dim(yprime),1))
  #trapopprime=rowSums(trapop)
  #trapopprime=matrix(trapopprime,nrow=length(trapopprime),ncol=1)
  #data$trapop=trapopprime
  #data$y.I.obs=yprime
  #if you want to initTrue
  #yprimefull=apply(data$y.I,c(1,2),sum)
  #yprimefull=array(yprimefull,dim=c(dim(yprimefull),1))
  #yPprime=apply(data$y.P,c(1,2),sum)
  #yPprime=array(yPprime,dim=c(dim(yPprime),1))
  #data$y.I=yprimefull
  #data$y.P=yPprime
  
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
# Alternative for missing data files
bigout <- list()
for(i in 1:100){
  x <- out[[i]]
  if(all(is.na(x))){
    bigout[[i]] <- NA
    next
  }
  if(any(is.na(list(x$chain1, x$chain2, x$chain3)))){
    bigout[[i]] <- NA
  } else {
    bigout[[i]] <- mcmc.list(mcmc(x$chain1$out), mcmc(x$chain2$out), mcmc(x$chain3$out))
  }
}
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
      for(countattract in c(-0.1)){
      #dir.create(paste("np", countone, "/lamp", counttwo, "/lami", countthree, sep=""))
      for(countfour in 1:100){ # change as you add more reps
        tempout <- bigout[[countfour]] # right? -30? Think
        save(tempout,file=paste("D:/Workspace/emmetr/newest_october20_awd/Ben Analyses/sims/np", countone, "/lamp", counttwo, "/lami", countthree, "/attr", abs(countattract), "/rep", countfour, "/knownidhmmresults.RData", sep=""))
        #tempidout <- bigidout[[countfour]]
        #save(tempidout,file=paste("np", countone, "/lamp", counttwo, "/lami", countthree, "/rep", countfour, "/idsimresults.RData", sep=""))
        
      }
      }
    }
  }
}
