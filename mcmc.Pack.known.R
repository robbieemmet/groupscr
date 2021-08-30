dztpois=function(k,lambda){
  return(lambda^k/((exp(lambda)-1)*factorial(k)))
}
e2dist=function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

mcmc.Pack.known<-
  function(data,niter=5000,nburn=0, nthin=1, Z = NA,W=NA, inits=NA, proppars=NA,mincluster=NA,
           storeLatent=list(ID=TRUE, s=TRUE, z=TRUE),plotLatent=plotLatent){
    if(length(storeLatent)!=3)stop("storeLatent should be a list of length 3, with elements 'ID', 's', and 'z'")
    
    library(abind)
    library(VGAM)
    buff<- data$buff
    X<-as.matrix(data$X)
    trapop <- data$trapop
    xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
    ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
    zeta=inits$zeta
    lam0=inits$lam0
    sigma=inits$sigma
    theta=inits$theta
    psi1=inits$psi1
    psi=inits$psi
    y.I=data$y.I.obs
    ni=nrow(y.I)
    J=dim(y.I)[2]
    K=dim(y.I)[3]
    if(is.null(trapop)){
      trapop <- array(1,dim=c(Z,J,K))
      warning("No trap file in use. Assuming perfect trap operation.")
    }else{
      if(!all(dim(trapop)==c(J,K)))stop("Trap operation must be a matrix of dimension J x K")
      trapop=array(trapop,dim=c(J,K,Z))
      trapop=aperm(trapop,c(3,1,2))
    }
    
    known.vector=c(rep(1,ni),rep(0,W-ni))
    #Augment y.I
    y.I<- abind(y.I,array(0, dim=c( W-dim(y.I)[1],J, K)), along=1)
    w=1*(apply(y.I,1,sum)>0)
    w[w==0]=rbinom(sum(w==0),1,psi1)
    
    #initialize pack activity centers
    s=cbind(runif(Z,xlim[1],xlim[2]),runif(Z,ylim[1],ylim[2]))
    
    #assign captured individuals to known packs
    gammai=rep(NA,W)
    gammai[1:ni]=data$packID[which(rowSums(data$y.I)>0)]
    cappacks=unique(gammai)
    cappacks=cappacks[-which(is.na(cappacks))]
    npack=length(cappacks)

    #assign uncaptured individuals to packs
    idx=unique(gammai[1:ni])
    idx2=ni+1
    for(i in idx){
      assigned=sum(gammai==i,na.rm=TRUE)
      add=mincluster-assigned
      if(add>0){
        gammai[idx2:(idx2+add-1)]=i
        idx2=idx2+add
      }
    }
    empty=setdiff(1:Z,unique(gammai[1:(idx2-1)]))
    for(i in empty){
      gammai[idx2:(idx2+mincluster-1)]=i
      idx2=idx2+mincluster
    }
    maxZ=sum(!is.na(gammai))
    if(maxZ>=W)stop(paste("must set W to at least",maxZ+1,"given Z and mincluster"))
    #Randomly assign remaining individuals in a pack
    gammai[idx2:W]=sample(1:Z,length(idx2:W),replace=TRUE)
    
    #make sure packs with captured individuals turned on
    z=rep(0,Z)
    z[unique(gammai)]=1
    #turn off empty packs
    rem=setdiff( which(z==1),unique(gammai[w==1]))
    z[rem]=0
    
    #improve initial s for captured packs
    y.I2D=apply(y.I,c(1,2),sum)
    for(i in 1:length(cappacks)){
      trapcaps=which(y.I2D[cappacks[i],]>0)
      if(length(trapcaps)>1){
        meanloc=colMeans(X[trapcaps,])
      }else{
        meanloc=X[trapcaps,]
      }
      s[cappacks[i],]=meanloc
    }
    
    
    D<- array(rep(e2dist(s, X),K),dim=c(Z,J,K))
    lamd=lamd.cand=lam0*exp(-D*D/(2*sigma*sigma))
    v=array(0,dim=c(Z,J,K))
    for(i in 1:Z){
      for(j in 1:J){
        for(k in 1:K){
          v[i,j,k]=rpois(1,lamd[i,j,k]*trapop[i,j,k])
        }
      }
    }
    #make sure v>0 if sample assigned to it
    for(i in 1:ni){
      idx=which(y.I[i,,]>0,arr.ind=TRUE)
      if(K>1){
        for(j in 1:nrow(idx)){
          if(v[gammai[i],idx[j,1],idx[j,2]]==0){
            v[gammai[i],idx[j,1],idx[j,2]]=1
          }
          
        }
      }else{
        for(j in 1:length(idx)){
          if(v[gammai[i],idx[j],1]==0){
            v[gammai[i],idx[j],1]=1
          }
        }
      }
    }
    
    #pack site use likelihood
    D<- array(rep(e2dist(s, X),K),dim=c(Z,J,K)) # UPDATE D
    lamd=lamd.cand=lam0*exp(-D*D/(2*sigma*sigma))
    ll.v=dpois(v,lamd,log=TRUE)
    ll.v.cand=ll.v
    ll.v.sum=sum(ll.v)
    if(!is.finite(ll.v.sum))stop("Starting pack site visitation likelihood is not finite. Try raising lam0 or sigma")
    
    #individual likelihood given pack site use
    lambda.tmp=w*v[gammai,,]*theta #reorder by individual gammai
    ll.yi=dpois(y.I,lambda.tmp,log=TRUE)
    ll.yi.cand=ll.yi
    ll.yi.sum=sum(ll.yi)
    
    ll.pack=n_g=rep(0,Z)
    for(i in 1:Z){
      if(z[i]==1){
        n_g[i]=sum(gammai==i&w==1)
        ll.pack[i]=log(dztpois(n_g[i],zeta))
      }else{
        n_g[i]=rzapois(1,zeta)
        ll.pack[i]=log(dztpois(n_g[i],zeta))
      }
    }
    ll.pack.cand=ll.pack
    ll.pack.sum=sum(ll.pack)
    n_g.cand=n_g
    maxout=0
    
    # some objects to hold the MCMC simulation output
    nstore=(niter-nburn)/nthin
    if(nburn%%nthin!=0){
      nstore=nstore+1
    }
    out<-matrix(NA,nrow=nstore,ncol=8)
    dimnames(out)<-list(NULL,c("zeta","lam0","sigma","theta","psi","Npack","Nind","SiteVisits"))
    if(storeLatent$ID){
      IDout<- matrix(NA, nrow=nstore, ncol=W)
    }else{
      IDout=NA
    }
    if(storeLatent$s){
      sout<- array(NA, dim=c(nstore,Z,2))
    }else{
      sout=NA
    }
    if(storeLatent$z){
      wout<- matrix(NA,nrow=nstore,ncol=W)
      zout<- matrix(NA,nrow=nstore,ncol=Z)
    }else{
      wout=zout=NA
    }
    iteridx=1 #for storing output not recorded every iteration
    
    for(iter in 1:niter){
      #Update lam0
      ll.v.sum=sum(ll.v)
      lam0.cand<- rnorm(1,lam0,proppars$lam0)
      if(lam0.cand > 0){
        lamd.cand<- lam0.cand*exp(-D*D/(2*sigma*sigma))
        ll.v.cand=dpois(v,lamd.cand*trapop,log=TRUE)
        ll.v.cand.sum=sum(ll.v.cand)
        if(runif(1) < exp(ll.v.cand.sum - ll.v.sum)){
          lam0=lam0.cand
          lamd=lamd.cand
          ll.v=ll.v.cand
          ll.v.sum=ll.v.cand.sum
        }
      }
      #Update sigma
      sigma.cand<- rnorm(1,sigma,proppars$sigma)
      if(lam0.cand > 0){
        lamd.cand<- lam0*exp(-D*D/(2*sigma.cand*sigma.cand))
        ll.v.cand=dpois(v,lamd.cand*trapop,log=TRUE)
        ll.v.cand.sum=sum(ll.v.cand)
        if(runif(1) < exp(ll.v.cand.sum - ll.v.sum)){
          sigma<- sigma.cand
          lamd=lamd.cand
          ll.v=ll.v.cand
        }
      }
      #update theta, individual use conditional on site use
      theta.cand<- rnorm(1,theta,proppars$theta)
      ll.yi.cand=ll.yi
      if(theta.cand > 0){
        lambda.tmp=w*v[gammai,,]*theta.cand
        ll.yi.cand=dpois(y.I,lambda.tmp,log=TRUE)
        ll.yi.cand.sum=sum(ll.yi.cand)
        ll.yi.sum=sum(ll.yi)
        if(runif(1) < exp(ll.yi.cand.sum-ll.yi.sum)){
          theta<- theta.cand
          ll.yi=ll.yi.cand
        }
      }
      #Update zeta for mean pack size
      zeta.cand<- rnorm(1,zeta,proppars$zeta)
      if(zeta.cand > 0){
        ll.pack.cand=log(dztpois(n_g,zeta.cand))
        ll.pack.cand.sum=sum(ll.pack.cand)
        ll.pack.sum=sum(ll.pack)
        if(runif(1) < exp(ll.pack.cand.sum-ll.pack.sum)){
          zeta<- zeta.cand
          ll.pack=ll.pack.cand
        }
      }
      #update n_g for z[i]=0
      zeros=which(z==0)
      if(length(zeros)>0){
        n_g[zeros]=rzapois(length(zeros),zeta)
        ll.pack[zeros]=log(dztpois(n_g[zeros],zeta))
      }
      
      #update latent pack use, v, using Rcpp below
      #  # First for z=1
      #  v.cand=v
      #  for(i in 1:Z){
      #    if(z[i]==0)next
      #    guys=which(gammai==i)
      #    sampmat=matrix(sample(c(-1,1),J*K,replace=TRUE),nrow=J,ncol=K)#faster to draw at once
      #    for(j in 1:J){
      #      for(k in 1:K){
      #        v.cand[i,j,k]=v[i,j,k]+sampmat[j,k]#Randomly subtract or add 1
      #        if(v.cand[i,j,k]<0)next
      #        obs=y.I[guys,j,k]
      #        if(any(obs>0)&v.cand[i,j,k]==0)next#test for v.cand=0 and obs>0
      #        ll.yi.cand[guys,j,k]=dpois(y.I[guys,j,k],w[guys]*v.cand[i,j,k]*theta,log=TRUE)
      #        ll.v.cand[i,j,k]=dpois(v.cand[i,j,k],lamd[i,j],log=TRUE)
      #        if(runif(1) < exp((sum(ll.yi.cand[guys,j,k])+ll.v.cand[i,j,k])-
      #                          (sum(ll.yi[guys,j,k])+ll.v[i,j,k]))){
      #          ll.yi[guys,j,k]=ll.yi.cand[guys,j,k]
      #          v[i,j,k]=v.cand[i,j,k]
      #          ll.v[i,j,k]=ll.v.cand[i,j,k]
      #        }
      #      }
      #    }
      #  }
      # 
      # #update latent pack use for v, z=0. simulate from full conditional
      # for(i in 1:Z){
      #   if(z[i]==1)next
      #   for(j in 1:J){
      #     v[i,j,]=rpois(K,lamd[i,j])
      #   }
      #   ll.v[i,,]=dpois(v[i,,],lamd[i,],log=TRUE)
      # }
      
      # update latent pack use in Rcpp
      sampmat=array(sample(c(-1,1),J*K*Z,replace=TRUE),dim=c(Z,J,K))
      a=siteup3D(w,z,y.I,v,sampmat,theta,gammai-1,ll.yi,ll.v,lamd,trapop)
      
      
      #update z, which inds in population
      #z update
      for(i in 1:Z){
        z.cand=z
        z.cand[i]=1-z[i]
        #priors
        prior=dbinom(z[i], 1, psi, log = TRUE)
        prior.cand=dbinom(z.cand[i], 1, psi, log = TRUE)
        if(z[i]==1){#turn off if all augmented individuals
          onpack=which(gammai==i&w==1)
          #if all on pack members are augmented, turn off z and all w
          if(all(known.vector[onpack]==0)){#test for all pack members augmented
            w.cand=w
            w.cand[onpack]=0
            ll.pack.cand[i]=log(dztpois(n_g.cand[i],zeta))
            ll.yi.cand[onpack,,]=0
            ll.yi.cand.sum=0#Don't need to sum a bunch of 0s
            ll.yi.sum=sum(ll.yi[onpack,,])
            if(runif(1) < exp((ll.yi.cand.sum+prior.cand) -
                              (ll.yi.sum+prior))){
              w[onpack]=w.cand[onpack]
              z[i]=z.cand[i]
              ll.yi[onpack,,]=0
            }
          }#otherwise cannot turn off
        }else{#turn on
          offpack=which(gammai==i&w==0)
          Nonu=n_g[i] #latent structure count already there
          n_g.cand[i]=n_g[i]
          if(Nonu>length(offpack)){
            Nonu=length(offpack)
            n_g.cand[i]=Nonu
            warning(paste("Not enough ws to turn z on with correct n_g[i] on iteration",iter,". Probably should augment more."))
            maxout=maxout+1
          }
          onpack=sample(offpack,Nonu)
          w.cand=w
          w.cand[onpack]=1
          #cluster LL
          ll.pack.cand[i]=log(dztpois(Nonu,zeta))
          #ll.yi
          lambda.tmp=v[gammai[onpack],,]*theta #pull out site visits for this pack, multiply by theta
          ll.yi.cand[onpack,,]=dpois(y.I[onpack,,],lambda.tmp,log=TRUE)
          
          ll.yi.cand.sum=sum(ll.yi.cand[onpack,,])
          ll.yi.sum=0 #Don't need to sum a bunch of 0s
          #count likelihood does not change unless you get the warning message above
          if (runif(1) < exp((ll.yi.cand.sum+prior.cand) -
                             (ll.yi.sum+prior))) {
            z[i] <- z.cand[i]
            w[onpack]=w.cand[onpack]
            ll.yi[onpack,,]=ll.yi.cand[onpack,,]
            n_g[i]=n_g.cand[i]
            ll.pack[i]=ll.pack.cand[i]
          }
        }
      }
      
      psi <- rbeta(1, 1 + sum(z), 1 + Z - sum(z))
      
      #update w, which inds in population
      for (i in 1:W) {
        w.cand=w
        if(known.vector[i]==1)next
        w.cand[i]=1-w[i]
        thispack=gammai[i]
        if(z[thispack]==0)next #only turn on or off w if z already on
        curr.members=which(gammai[w==1]==thispack)
        new.members=which(gammai[w.cand==1]==thispack)
        n_g.cand[thispack]=length(new.members)
        if(n_g.cand[thispack]==0)next #can't have 0 pack members.
        #cluster likelihood
        ll.pack.cand[thispack]=log(dztpois(n_g.cand[thispack],zeta))
        #individual likelihood
        if(w.cand[i]==1){
          lambda.tmp=v[thispack,,]*theta #pull out site visits for this pack, multiply by theta
          ll.yi.cand[i,,]=dpois(y.I[i,,],lambda.tmp,log=TRUE)
          ll.yi.cand.sum=sum(ll.yi.cand[i,,])
          ll.yi.sum=0 #Don't need to sum 0s
        }else{
          ll.yi.cand[i,,]=0
          ll.yi.cand.sum=0 #Don't need to sum 0s
          ll.yi.sum=sum(ll.yi[i,,])
        }
        idx2=which(gammai==thispack)
        #proposal probabilities for counts
        pfor=sum(w[idx2]==w[i])/length(idx2) #prob propose w.cand[i]
        pback=sum(w.cand[idx2]==w.cand[i])/length(idx2) #prob get back to w[i]
        
        if(runif(1) < exp((ll.yi.cand.sum+ll.pack.cand[thispack]) -
                          (ll.yi.sum+ll.pack[thispack]))*(pback/pfor)){
          w[i]=w.cand[i]
          n_g[thispack]=n_g.cand[thispack]
          ll.pack[thispack]=ll.pack.cand[thispack]
          ll.yi[i,,]=ll.yi.cand[i,,]
        }
      }
      
      ##update the activity centers
      for (i in 1:Z) {
        Scand <- c(rnorm(1, s[i, 1], proppars$s), rnorm(1, s[i, 2], proppars$s))
        inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
        if (inbox) {
          dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
          dtmp <- matrix(dtmp, nrow=J,ncol=K,byrow=T)
          lamd.cand[i,,]=lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
          ll.v.cand[i,,]=dpois(v[i,,],lamd.cand[i,,]*trapop[i,,],log=TRUE)
          if (runif(1) < exp(sum(ll.v.cand[i,,])- sum(ll.v[i,,]))) {
            s[i,]=Scand
            D[i,,]=dtmp
            lamd[i,,]=lamd.cand[i,,]
            ll.v[i,,]=ll.v.cand[i,,]
          }
        }
      }
      #Do we record output on this iteration?
      if(iter>nburn&iter%%nthin==0){
        if(storeLatent$ID){
          IDout[iteridx,]=gammai
        }
        if(storeLatent$s){
          sout[iteridx,,]=s
        }
        if(storeLatent$z){
          wout[iteridx,]=w
          zout[iteridx,]=z
        }
        out[iteridx,]<- c(zeta,lam0,sigma,theta,psi,sum(z),sum(w),sum(v[z==1,,]))
        iteridx=iteridx+1
      }
      #Do we plot?
      if(plotLatent){
        if(iter%%100==0){
          par(mfrow=c(1,1),ask=FALSE)
          par(mar=c(8.1, 4.1, 4.1, 2.1), xpd=TRUE)
          plot(NA,pch=16,xlim=xlim,ylim=ylim,xlab="X",ylab="Y",
               main=paste("N Packs = ", sum(z), ", n Packs = ", npack, ", N Ind = ", sum(w)))
          visiti=which(rowSums(v)>0)
          v2D=apply(v,c(1,2),sum)
          for(i in visiti){
            if(z[i]==0)next
            trapcap=which(v2D[i,]>0)
            for(j in trapcap){
              lines(x=c(s[i,1],X[j,1]),y=c(s[i,2],X[j,2]),lty=3)
            }
          }
          capi=which(rowSums(y.I)>0)
          y.I2D=apply(y.I,c(1,2),sum)
          for(i in capi){
            trapcap=which(y.I2D[i,]>0)
            for(j in trapcap){
              lines(x=c(s[gammai[i],1],X[j,1]),y=c(s[gammai[i],2],X[j,2]),lty=1,col="#999999",lwd=2)
            }
          }
          points(s[z==1,],cex=2.5,col="#56B4E9",pch=16)
          points(s[z==0,],cex=2.5,col="#D55E00",pch=16)
          text(s[,1],s[,2],label=as.character(n_g),cex=0.75)
          points(X,pch=4)
          legend("bottomleft", inset=c(0,-0.25), legend=c("True","Augmented"), pch=c(16,16),
                 col=c("#56B4E9","#D55E00"),title="Packs",cex=0.85)
          legend("bottomright", inset=c(-0.05,-0.25), legend=c("Detector","Site Visit","Pack Member Detection"), pch=c(4,NA,NA),
                 col=c("black","black","#999999"),lty=c(NA,3,1),lwd=c(NA,1,2),cex=0.85)
          
          
          
        }
      }
    }
    
    if(any(c(storeLatent$ID, storeLatent$s, storeLatent$z))){
      out2=list(out=out,sout=sout,wout=wout,zout=zout)
    }else{ 
      out2=list(out=out)
    }
    return(out2)
  }
