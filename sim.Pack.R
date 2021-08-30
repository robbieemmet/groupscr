
sim.Pack=function(Np=10,zeta=6,lam0=0.2,sigma=0.50,theta=0.2,K=10,X=X,buff=3){
  library(VGAM)
  # In here, packID is the observed/realized "gamma" ("gammai" in the model)
  # # simulate a population of group activity centers
  s=cbind(runif(Np, min(X[,1])-buff,max(X[,1])+buff), runif(Np,min(X[,2])-buff,max(X[,2])+buff))
  D=e2dist(s,X)
  lamd<- lam0*exp(-D*D/(2*sigma*sigma))
  J<- nrow(X)
  # Simulate latent group encounter history
  v=array(0,dim=c(Np,J,K))
  for(i in 1:Np){
    for(j in 1:J){
      for(k in 1:K){
        v[i,j,k]=rpois(1,lamd[i,j])
      }
    }
  }
  #Pack membership stuff
  NperPack=rzapois(Np,zeta)
  Ni=sum(NperPack)
  packID=rep(1:Np,times=NperPack)
  
  #Simulate observed individual history
  y.I=array(0,dim=c(Ni,J,K))
  for(i in 1:Ni){
    for(j in 1:J){
      for(k in 1:K){
        y.I[i,j,k]=rpois(1,theta*v[packID[i],j,k])
      }
    }
  }
  si=cbind(rep(s[,1],times=NperPack),rep(s[,2],times=NperPack))
  keep=which(rowSums(y.I)>0)
  y.I.obs=y.I[keep,,]
  if(K==1){
    y.I.obs=array(y.I.obs,dim=c(dim(y.I.obs),K))
  }
  si.obs=si[keep,]
  
  
  out<-list(v=v,s=s,X=X, K=K,y.I=y.I,si=si,y.I.obs=y.I.obs,si.obs=si.obs,NperPack=NperPack,
            packID=packID,buff=buff,Ni=Ni)
  return(out)
}