library(coda)
library(MCMCglmm)
# Keep only converged models and models with non-missing Rhat
# Thin and burnin for both models
# Use 2.5% and 97.5% quantiles to assess coverage
scenarios=data.frame(np=c(10,10,20,20),lamp=c(5,12.5,5,12.5),lami=rep(0.3,4))
setwd("D:/Workspace/emmetr/newest_october20_awd/Ben Analyses/")
#Scenario 1
scen=1
true=c(5,0.034,4.5,0.3,10/60,10,NA,NA)
modes=medians=means=covers=groupmodes=groupcovers=groupeffective=groupgelmans=effective=gelmans=matrix(NA,nrow=100,ncol=8)
HPDs=creds=groupcreds=array(NA,dim=c(100,8,2))
storeSV=storenp=storeNi=storeni=ind_recaptures=pack_recaptures=rep(NA,100)
for(i in 1:100){
  load(paste("sims/np", scenarios[scen,1], "/lamp", scenarios[scen,2], "/lami", scenarios[scen,3], "/rep", i, "/knownidsimresults.RData", sep=""))
  if(all(is.na(tempout)))next
  load(paste("sims/np", scenarios[scen,1], "/lamp", scenarios[scen,2], "/lami", scenarios[scen,3], "/rep", i, "/simdata.RData", sep=""))
  
  # Ind and group recaptures
  pack_by_trap <- aggregate(apply(data$y.I.obs,c(1,2),sum), by=list(packID=data$packID[1:dim(data$y.I.obs)[1]]), sum)
  pack_by_trap <- pack_by_trap[,-1]
  
  # Pack recaptures and individual recaptures
  inds_captured <- apply(data$y.I.obs, c(1,2), sum)>0
  ind_recaptures[i] <- mean(apply(inds_captured, 1, sum))
  
  pack_recaptures[i] <- mean(apply(pack_by_trap>0, 1, sum))
  
  gelmans[i,]=gelman.diag(tempout)$psrf[,1]
  if(max(gelmans[i,]>1.1)){next}
  effective[i,]=effectiveSize(tempout) # before burn-in and thinning
  #pool chains
  tempout=mcmc(rbind(tempout[[1]][seq(5001,30000,by=150),],tempout[[2]][seq(5001,30000,by=150),],tempout[[3]][seq(5001,30000,by=150),]))
  true[7]=storeNi[i]=data$Ni
  true[8]=storeSV[i]=sum(data$y.P)
  storeni[i]=nrow(data$y.I.obs)
  modes[i,]=posterior.mode(tempout)
  medians[i,]=summary(tempout)$quantiles[,c("50%")]
  means[i,]=summary(tempout)$statistics[, "Mean"]
  #HPDs[i,,]=HPDinterval(tempout)
  creds[i,,]=summary(tempout)$quantiles[,c("2.5%", "97.5%")]
  covers[i,]=1*(creds[i,,1]<=true&creds[i,,2]>=true)
  
  # Group SCR?
  load(paste("sims/np", scenarios[scen,1], "/lamp", scenarios[scen,2], "/lami", scenarios[scen,3], "/rep", i, "/basic_group_scr.RData", sep=""))
  if(all(is.na(tempout))){
    next
    }
  tempgelman=unlist(tempout$Rhat)#tryCatch(gelman.diag(tempout$samples)$psrf[,1],error=function(e){tempgelman<-NA})
  if(any(is.na(tempgelman))){
    groupgelmans[i,]=rep(NA,8)
    maxgelman <- 1.2 # just to make sure values are saved here?
    next
  } else {
    groupgelmans[i,]=c(NA, tempgelman["lambda0"], tempgelman["sigma"], NA, tempgelman["psi"], tempgelman["N"], NA, NA)
    maxgelman <- max(tempgelman)
    }
  if(maxgelman>1.1){next}
  tempeff=unlist(tempout$n.eff) #effectiveSize(tempout$samples)
  groupeffective[i,]=c(NA, tempeff["lambda0"], tempeff["sigma"], NA, tempeff["psi"], tempeff["N"], NA, NA) # before burn-in and thinning
  #pool chains
  tempout=mcmc(rbind(tempout$samples[[1]][seq(5001,30000,by=150),],tempout$samples[[2]][seq(5001,30000,by=150),],tempout$samples[[3]][seq(5001,30000,by=150),]))
  #true[7]=storeNi[i]=data$Ni
  #true[8]=storeSV[i]=sum(data$y.P)
  #storeni[i]=nrow(data$y.I.obs)
  tempmodes=posterior.mode(tempout)
  groupmodes[i,]=c(NA, tempmodes["lambda0"], tempmodes["sigma"], NA, tempmodes["psi"], tempmodes["N"], NA, NA)
  #HPDs[i,,]=HPDinterval(tempout)
  tempcreds=summary(tempout)$quantiles[,c("2.5%", "97.5%")]
  groupcreds[i,,]=rbind(rep(NA,2), tempcreds["lambda0",], tempcreds["sigma",], rep(NA,2), tempcreds["psi",], tempcreds["N",], rep(NA,2), rep(NA,2))
  groupcovers[i,]=1*(groupcreds[i,,1]<=true&groupcreds[i,,2]>=true)
}

# Convergence/model failure counts
which(is.na(rowSums(modes))) # 12 convergence or model failures
unique(which(gelmans>1.1, arr.ind=T)[,1]) # All convergence failures no model failures

which(rowSums(groupmodes, na.rm=T)==0) # 21 convergence or model failures
unique(which(groupgelmans>1.1, arr.ind=T)[,1]) # All model failures

ests1=rbind(true,colMeans(modes,na.rm=TRUE), colMeans(groupmodes, na.rm=TRUE), colMeans(covers,na.rm=TRUE), colMeans(groupcovers, na.rm=TRUE), colMeans(effective,na.rm=TRUE), colMeans(groupeffective, na.rm=TRUE))
ests1[1,7:8]=c(mean(storeNi,na.rm=TRUE),mean(storeSV,na.rm=TRUE))
colnames(ests1)=c("lambda.P","lam0.P","sigma.P","lambda.I","psi","Npack","Nind","SV")
rownames(ests1)=c("true","est", "groupest", "cover","groupcover", "eff", "groupeff")
round(ests1,3)
biasrow <- ((ests1["est",]-ests1["true",])/ests1["true",])*100
biasrow
ests1["cover",]

biasrowgroup <- ((ests1["groupest",]-ests1["true",])/ests1["true",])*100
biasrowgroup
ests1["groupcover",]

mean(ind_recaptures, na.rm=T)
mean(pack_recaptures, na.rm=T)


#Scenario 2
scen=2
true=c(12.5,0.034,4.5,0.3,10/60,10,NA,NA)
modes=covers=medians=means=groupmodes=groupcovers=groupeffective=groupgelmans=effective=gelmans=matrix(NA,nrow=100,ncol=8)
HPDs=creds=groupcreds=array(NA,dim=c(100,8,2))
storeSV=storenp=storeNi=storeni=ind_recaptures=pack_recaptures=rep(NA,100)
for(i in 1:100){
  load(paste("sims/np", scenarios[scen,1], "/lamp", scenarios[scen,2], "/lami", scenarios[scen,3], "/rep", i, "/knownidsimresults.RData", sep=""))
  if(all(is.na(tempout)))next
  load(paste("sims/np", scenarios[scen,1], "/lamp", scenarios[scen,2], "/lami", scenarios[scen,3], "/rep", i, "/simdata.RData", sep=""))
  
  # Ind and group recaptures
  pack_by_trap <- aggregate(apply(data$y.I.obs,c(1,2),sum), by=list(packID=data$packID[1:dim(data$y.I.obs)[1]]), sum)
  pack_by_trap <- pack_by_trap[,-1]
  
  # Pack recaptures and individual recaptures
  inds_captured <- apply(data$y.I.obs, c(1,2), sum)>0
  ind_recaptures[i] <- mean(apply(inds_captured, 1, sum))
  
  pack_recaptures[i] <- mean(apply(pack_by_trap>0, 1, sum))
  
  gelmans[i,]=gelman.diag(tempout)$psrf[,1]
  if(max(gelmans[i,]>1.1)){next}
  effective[i,]=effectiveSize(tempout) # before burn-in and thinning
  #pool chains
  tempout=mcmc(rbind(tempout[[1]][seq(5001,30000,by=150),],tempout[[2]][seq(5001,30000,by=150),],tempout[[3]][seq(5001,30000,by=150),]))
  true[7]=storeNi[i]=data$Ni
  true[8]=storeSV[i]=sum(data$y.P)
  storeni[i]=nrow(data$y.I.obs)
  modes[i,]=posterior.mode(tempout)
  medians[i,]=summary(tempout)$quantiles[,c("50%")]
  means[i,]=summary(tempout)$statistics[, "Mean"]
  #HPDs[i,,]=HPDinterval(tempout)
  creds[i,,]=summary(tempout)$quantiles[,c("2.5%", "97.5%")]
  covers[i,]=1*(creds[i,,1]<=true&creds[i,,2]>=true)
  
  # Group SCR?
  load(paste("sims/np", scenarios[scen,1], "/lamp", scenarios[scen,2], "/lami", scenarios[scen,3], "/rep", i, "/basic_group_scr.RData", sep=""))
  if(all(is.na(tempout))){
    next
  }
  tempgelman=unlist(tempout$Rhat)#tryCatch(gelman.diag(tempout$samples)$psrf[,1],error=function(e){tempgelman<-NA})
  if(any(is.na(tempgelman))){
    groupgelmans[i,]=rep(NA,8)
    maxgelman <- 1.2 # just to make sure values are saved here?
    next
  } else {
    groupgelmans[i,]=c(NA, tempgelman["lambda0"], tempgelman["sigma"], NA, tempgelman["psi"], tempgelman["N"], NA, NA)
    maxgelman <- max(tempgelman)
  }
  if(maxgelman>1.1){next}
  tempeff=unlist(tempout$n.eff) #effectiveSize(tempout$samples)
  groupeffective[i,]=c(NA, tempeff["lambda0"], tempeff["sigma"], NA, tempeff["psi"], tempeff["N"], NA, NA) # before burn-in and thinning
  #pool chains
  tempout=mcmc(rbind(tempout$samples[[1]][seq(5001,30000,by=150),],tempout$samples[[2]][seq(5001,30000,by=150),],tempout$samples[[3]][seq(5001,30000,by=150),]))
  #true[7]=storeNi[i]=data$Ni
  #true[8]=storeSV[i]=sum(data$y.P)
  #storeni[i]=nrow(data$y.I.obs)
  tempmodes=posterior.mode(tempout)
  groupmodes[i,]=c(NA, tempmodes["lambda0"], tempmodes["sigma"], NA, tempmodes["psi"], tempmodes["N"], NA, NA)
  #HPDs[i,,]=HPDinterval(tempout)
  tempcreds=summary(tempout)$quantiles[,c("2.5%", "97.5%")]
  groupcreds[i,,]=rbind(rep(NA,2), tempcreds["lambda0",], tempcreds["sigma",], rep(NA,2), tempcreds["psi",], tempcreds["N",], rep(NA,2), rep(NA,2))
  groupcovers[i,]=1*(groupcreds[i,,1]<=true&groupcreds[i,,2]>=true)
}

# Convergence/model failure counts
which(is.na(rowSums(modes))) # 2 convergence or model failures
unique(which(gelmans>1.1, arr.ind=T)[,1]) # All convergence failures no model failures

which(rowSums(groupmodes, na.rm=T)==0) # 17 convergence or model failures
unique(which(groupgelmans>1.1, arr.ind=T)[,1]) # 2 are non convergence

ests1=rbind(true,colMeans(modes,na.rm=TRUE), colMeans(groupmodes, na.rm=TRUE), colMeans(covers,na.rm=TRUE), colMeans(groupcovers, na.rm=TRUE), colMeans(effective,na.rm=TRUE), colMeans(groupeffective, na.rm=TRUE))
ests1[1,7:8]=c(mean(storeNi,na.rm=TRUE),mean(storeSV,na.rm=TRUE))
colnames(ests1)=c("lambda.P","lam0.P","sigma.P","lambda.I","psi","Npack","Nind","SV")
rownames(ests1)=c("true","est", "groupest", "cover","groupcover", "eff", "groupeff")
round(ests1,3)
biasrow <- ((ests1["est",]-ests1["true",])/ests1["true",])*100
biasrow
ests1["cover",]

biasrowgroup <- ((ests1["groupest",]-ests1["true",])/ests1["true",])*100
biasrowgroup
ests1["groupcover",]

mean(ind_recaptures, na.rm=T)
mean(pack_recaptures, na.rm=T)

#Scenario 3
scen=3
true=c(5,0.034,4.5,0.3,20/60,20,NA,NA)
modes=covers=medians=means=groupmodes=groupcovers=groupeffective=groupgelmans=effective=gelmans=matrix(NA,nrow=100,ncol=8)
HPDs=creds=groupcreds=array(NA,dim=c(100,8,2))
storeSV=storenp=storeNi=storeni=ind_recaptures=pack_recaptures=rep(NA,100)
for(i in 1:100){
  load(paste("sims/np", scenarios[scen,1], "/lamp", scenarios[scen,2], "/lami", scenarios[scen,3], "/rep", i, "/knownidsimresults.RData", sep=""))
  if(all(is.na(tempout)))next
  load(paste("sims/np", scenarios[scen,1], "/lamp", scenarios[scen,2], "/lami", scenarios[scen,3], "/rep", i, "/simdata.RData", sep=""))
  
  # Ind and group recaptures
  pack_by_trap <- aggregate(apply(data$y.I.obs,c(1,2),sum), by=list(packID=data$packID[1:dim(data$y.I.obs)[1]]), sum)
  pack_by_trap <- pack_by_trap[,-1]
  
  # Pack recaptures and individual recaptures
  inds_captured <- apply(data$y.I.obs, c(1,2), sum)>0
  ind_recaptures[i] <- mean(apply(inds_captured, 1, sum))
  
  pack_recaptures[i] <- mean(apply(pack_by_trap>0, 1, sum))
  
  gelmans[i,]=gelman.diag(tempout)$psrf[,1]
  if(max(gelmans[i,]>1.1)){next}
  effective[i,]=effectiveSize(tempout) # before burn-in and thinning
  #pool chains
  tempout=mcmc(rbind(tempout[[1]][seq(5001,30000,by=150),],tempout[[2]][seq(5001,30000,by=150),],tempout[[3]][seq(5001,30000,by=150),]))
  true[7]=storeNi[i]=data$Ni
  true[8]=storeSV[i]=sum(data$y.P)
  storeni[i]=nrow(data$y.I.obs)
  modes[i,]=posterior.mode(tempout)
  medians[i,]=summary(tempout)$quantiles[,c("50%")]
  means[i,]=summary(tempout)$statistics[, "Mean"]
  #HPDs[i,,]=HPDinterval(tempout)
  creds[i,,]=summary(tempout)$quantiles[,c("2.5%", "97.5%")]
  covers[i,]=1*(creds[i,,1]<=true&creds[i,,2]>=true)
  
  # Group SCR?
  load(paste("sims/np", scenarios[scen,1], "/lamp", scenarios[scen,2], "/lami", scenarios[scen,3], "/rep", i, "/basic_group_scr.RData", sep=""))
  if(all(is.na(tempout))){
    next
  }
  tempgelman=unlist(tempout$Rhat)#tryCatch(gelman.diag(tempout$samples)$psrf[,1],error=function(e){tempgelman<-NA})
  if(any(is.na(tempgelman))){
    groupgelmans[i,]=rep(NA,8)
    maxgelman <- 1.2 # just to make sure values are saved here?
    next
  } else {
    groupgelmans[i,]=c(NA, tempgelman["lambda0"], tempgelman["sigma"], NA, tempgelman["psi"], tempgelman["N"], NA, NA)
    maxgelman <- max(tempgelman)
  }
  if(maxgelman>1.1){next}
  tempeff=unlist(tempout$n.eff) #effectiveSize(tempout$samples)
  
  groupeffective[i,]=c(NA, tempeff["lambda0"], tempeff["sigma"], NA, tempeff["psi"], tempeff["N"], NA, NA) # before burn-in and thinning
  #pool chains
  tempout=mcmc(rbind(tempout$samples[[1]][seq(5001,30000,by=150),],tempout$samples[[2]][seq(5001,30000,by=150),],tempout$samples[[3]][seq(5001,30000,by=150),]))
  #true[7]=storeNi[i]=data$Ni
  #true[8]=storeSV[i]=sum(data$y.P)
  #storeni[i]=nrow(data$y.I.obs)
  tempmodes=posterior.mode(tempout)
  groupmodes[i,]=c(NA, tempmodes["lambda0"], tempmodes["sigma"], NA, tempmodes["psi"], tempmodes["N"], NA, NA)
  #HPDs[i,,]=HPDinterval(tempout)
  tempcreds=summary(tempout)$quantiles[,c("2.5%", "97.5%")]
  groupcreds[i,,]=rbind(rep(NA,2), tempcreds["lambda0",], tempcreds["sigma",], rep(NA,2), tempcreds["psi",], tempcreds["N",], rep(NA,2), rep(NA,2))
  groupcovers[i,]=1*(groupcreds[i,,1]<=true&groupcreds[i,,2]>=true)
}

# Convergence/model failure counts
which(is.na(rowSums(modes))) # 1 convergence or model failure
unique(which(gelmans>1.1, arr.ind=T)[,1]) # All convergence failures no model failures

which(rowSums(groupmodes, na.rm=T)==0) # 7 convergence or model failures
unique(which(groupgelmans>1.1, arr.ind=T)[,1]) # All model failures

ests1=rbind(true,colMeans(modes,na.rm=TRUE), colMeans(groupmodes, na.rm=TRUE), colMeans(covers,na.rm=TRUE), colMeans(groupcovers, na.rm=TRUE), colMeans(effective,na.rm=TRUE), colMeans(groupeffective, na.rm=TRUE))
ests1[1,7:8]=c(mean(storeNi,na.rm=TRUE),mean(storeSV,na.rm=TRUE))
colnames(ests1)=c("lambda.P","lam0.P","sigma.P","lambda.I","psi","Npack","Nind","SV")
rownames(ests1)=c("true","est", "groupest", "cover","groupcover", "eff", "groupeff")
round(ests1,3)
biasrow <- ((ests1["est",]-ests1["true",])/ests1["true",])*100
biasrow
ests1["cover",]

biasrowgroup <- ((ests1["groupest",]-ests1["true",])/ests1["true",])*100
biasrowgroup
ests1["groupcover",]

mean(ind_recaptures, na.rm=T)
mean(pack_recaptures, na.rm=T)

#Scenario 4
scen=4
true=c(12.5,0.034,4.5,0.3,20/60,20,NA,NA)
modes=covers=groupmodes=medians=means=groupcovers=groupeffective=groupgelmans=effective=gelmans=matrix(NA,nrow=100,ncol=8)
HPDs=creds=groupcreds=array(NA,dim=c(100,8,2))
storeSV=storenp=storeNi=storeni=ind_recaptures=pack_recaptures=rep(NA,100)
for(i in 1:100){
  load(paste("sims/np", scenarios[scen,1], "/lamp", scenarios[scen,2], "/lami", scenarios[scen,3], "/rep", i, "/knownidsimresults.RData", sep=""))
  if(all(is.na(tempout)))next
  load(paste("sims/np", scenarios[scen,1], "/lamp", scenarios[scen,2], "/lami", scenarios[scen,3], "/rep", i, "/simdata.RData", sep=""))
  
  # Ind and group recaptures
  pack_by_trap <- aggregate(apply(data$y.I.obs,c(1,2),sum), by=list(packID=data$packID[1:dim(data$y.I.obs)[1]]), sum)
  pack_by_trap <- pack_by_trap[,-1]

  # Pack recaptures and individual recaptures
  inds_captured <- apply(data$y.I.obs, c(1,2), sum)>0
  ind_recaptures[i] <- mean(apply(inds_captured, 1, sum))
  
  pack_recaptures[i] <- mean(apply(pack_by_trap>0, 1, sum))
  
  # Cluster SCR model
  gelmans[i,]=gelman.diag(tempout)$psrf[,1]
  if(max(gelmans[i,]>1.1)){next}
  effective[i,]=effectiveSize(tempout) # before burn-in and thinning
  #pool chains
  tempout=mcmc(rbind(tempout[[1]][seq(5001,30000,by=150),],tempout[[2]][seq(5001,30000,by=150),],tempout[[3]][seq(5001,30000,by=150),]))
  true[7]=storeNi[i]=data$Ni
  true[8]=storeSV[i]=sum(data$y.P)
  storeni[i]=nrow(data$y.I.obs)
  modes[i,]=posterior.mode(tempout)
  medians[i,]=summary(tempout)$quantiles[,c("50%")]
  means[i,]=summary(tempout)$statistics[, "Mean"]
  #HPDs[i,,]=HPDinterval(tempout)
  creds[i,,]=summary(tempout)$quantiles[,c("2.5%", "97.5%")]
  covers[i,]=1*(creds[i,,1]<=true&creds[i,,2]>=true)
  
  # Group SCR?
  load(paste("sims/np", scenarios[scen,1], "/lamp", scenarios[scen,2], "/lami", scenarios[scen,3], "/rep", i, "/basic_group_scr.RData", sep=""))
  if(all(is.na(tempout))){
    next
  }
  tempgelman=unlist(tempout$Rhat)#tryCatch(gelman.diag(tempout$samples)$psrf[,1],error=function(e){tempgelman<-NA})
  if(any(is.na(tempgelman))){
    groupgelmans[i,]=rep(NA,8)
    maxgelman <- 1.2 # just to make sure values are saved here?
    next
  } else {
    groupgelmans[i,]=c(NA, tempgelman["lambda0"], tempgelman["sigma"], NA, tempgelman["psi"], tempgelman["N"], NA, NA)
    maxgelman <- max(tempgelman)
  }
  if(maxgelman>1.1){next}
  tempeff=unlist(tempout$n.eff) #effectiveSize(tempout$samples)
  groupeffective[i,]=c(NA, tempeff["lambda0"], tempeff["sigma"], NA, tempeff["psi"], tempeff["N"], NA, NA) # before burn-in and thinning
  #pool chains
  tempout=mcmc(rbind(tempout$samples[[1]][seq(5001,30000,by=150),],tempout$samples[[2]][seq(5001,30000,by=150),],tempout$samples[[3]][seq(5001,30000,by=150),]))
  #true[7]=storeNi[i]=data$Ni
  #true[8]=storeSV[i]=sum(data$y.P)
  #storeni[i]=nrow(data$y.I.obs)
  tempmodes=posterior.mode(tempout)
  groupmodes[i,]=c(NA, tempmodes["lambda0"], tempmodes["sigma"], NA, tempmodes["psi"], tempmodes["N"], NA, NA)
  #HPDs[i,,]=HPDinterval(tempout)
  tempcreds=summary(tempout)$quantiles[,c("2.5%", "97.5%")]
  groupcreds[i,,]=rbind(rep(NA,2), tempcreds["lambda0",], tempcreds["sigma",], rep(NA,2), tempcreds["psi",], tempcreds["N",], rep(NA,2), rep(NA,2))
  groupcovers[i,]=1*(groupcreds[i,,1]<=true&groupcreds[i,,2]>=true)
}

# Convergence/model failure counts
which(is.na(rowSums(modes))) # 0 convergence or model failures
unique(which(gelmans>1.1, arr.ind=T)[,1]) # All convergence failures no model failures

which(rowSums(groupmodes, na.rm=T)==0) # 16 convergence or model failures
unique(which(groupgelmans>1.1, arr.ind=T)[,1]) # All model failures

ests1=rbind(true,colMeans(modes,na.rm=TRUE), colMeans(groupmodes, na.rm=TRUE), colMeans(covers,na.rm=TRUE), colMeans(groupcovers, na.rm=TRUE), colMeans(effective,na.rm=TRUE), colMeans(groupeffective, na.rm=TRUE))
ests1[1,7:8]=c(mean(storeNi,na.rm=TRUE),mean(storeSV,na.rm=TRUE))
colnames(ests1)=c("lambda.P","lam0.P","sigma.P","lambda.I","psi","Npack","Nind","SV")
rownames(ests1)=c("true","est", "groupest", "cover","groupcover", "eff", "groupeff")
round(ests1,3)
biasrow <- ((ests1["est",]-ests1["true",])/ests1["true",])*100
biasrow
ests1["cover",]

biasrowgroup <- ((ests1["groupest",]-ests1["true",])/ests1["true",])*100
biasrowgroup
ests1["groupcover",]

mean(ind_recaptures, na.rm=T)
mean(pack_recaptures, na.rm=T)
