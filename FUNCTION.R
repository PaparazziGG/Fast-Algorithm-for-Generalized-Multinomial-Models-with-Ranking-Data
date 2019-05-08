Transfer<-function(Data,Response)
{
  W_Indicator<-Data==1
  L_Indicator<-Data==-1
  P_Indicator<-Data!=0
  
  single_select<-apply(W_Indicator,1,sum)==1
  multi_select<-apply(W_Indicator,1,sum)!=1
  part_candidator<-apply(P_Indicator,1,prod)==0
  a<-rep(0,ncol(Data))
  b<-0
  Delta<-matrix(-1,nrow = 1,ncol=ncol(Data))
  if(sum(single_select)!=0)
  {a<-(t(W_Indicator[single_select,]*Response[single_select])%*%rep(1,sum(single_select)))[,1]}
  if(sum(multi_select)+sum(part_candidator)!=0)
    {b<-c(Response[multi_select],-Response[part_candidator])
  Delta<-rbind(W_Indicator[multi_select,],P_Indicator[part_candidator,])}
  
  return(list(Win=W_Indicator,Loss=L_Indicator,Participant=P_Indicator,a=a,b=b,Delta=Delta))
}

weaver<-function(a,b,Delta,tol=1e-6)
{
  if(prod(a)==0)
  {
    return("Invalid")
  }
  else
  {
    p_0_Weaver<-0
    p_1_Weaver<-rep(1/length(a),length(a))
    path_Weaver<-NULL
    path_Weaver<-rbind(path_Weaver,p_1_Weaver)
    k_Weaver<-0
    s=sum(a)+sum(b)
    
    while(sqrt(sum((p_1_Weaver-p_0_Weaver)^2))>tol)
    {
      p_0_Weaver<-p_1_Weaver
      tau<-b/((Delta%*%p_0_Weaver)[,1])
      p_1_Weaver<-a/(s-(t(Delta)%*%tau)[,1])
      p_1_Weaver<-p_1_Weaver/sum(p_1_Weaver)
      k_Weaver<-k_Weaver+1
      path_Weaver<-rbind(path_Weaver,p_1_Weaver)
    }
    return(list(Estimator=p_1_Weaver,Iteration=k_Weaver,Path=path_Weaver))
  }
}

Bayesian_weaver<-function(a,b,Delta,tol=1e-6)
{
  gamma<-sum(abs(b))
  p_0_Weaver<-0
  p_1_Weaver<-rep(1/length(a),length(a))
  path_Weaver<-NULL
  path_Weaver<-rbind(path_Weaver,p_1_Weaver)
  k_Weaver<-0
  s=sum(a)+sum(b)
  
  while(sqrt(sum((p_1_Weaver-p_0_Weaver)^2))>tol)
  {
    p_0_Weaver<-p_1_Weaver
    tau<-b/((Delta%*%p_0_Weaver)[,1])
    p_1_Weaver<-(a+p_0_Weaver*gamma)/(s+gamma-(t(Delta)%*%tau)[,1])
    p_1_Weaver<-p_1_Weaver/sum(p_1_Weaver)
    k_Weaver<-k_Weaver+1
    path_Weaver<-rbind(path_Weaver,p_1_Weaver)
    p_1_Weaver
  }
  return(list(Estimator=p_1_Weaver,Iteration=k_Weaver,Path=path_Weaver))
}

Bayesian_weaver2<-function(a,b,Delta,tol=1e-6)
{
  gamma<-2*sum(abs(b))
  p_0_Weaver<-0
  p_1_Weaver<-rep(1/length(a),length(a))
  path_Weaver<-NULL
  path_Weaver<-rbind(path_Weaver,p_1_Weaver)
  k_Weaver<-0
  s=sum(a)+sum(b)
  
  while(sqrt(sum((p_1_Weaver-p_0_Weaver)^2))>tol)
  {
    p_0_Weaver<-p_1_Weaver
    tau<-b/((Delta%*%p_0_Weaver)[,1])
    p_1_Weaver<-(a+p_0_Weaver*gamma)/(s+gamma-(t(Delta)%*%tau)[,1])
    p_1_Weaver<-p_1_Weaver/sum(p_1_Weaver)
    k_Weaver<-k_Weaver+1
    path_Weaver<-rbind(path_Weaver,p_1_Weaver)
    p_1_Weaver
  }
  return(list(Estimator=p_1_Weaver,Iteration=k_Weaver,Path=path_Weaver))
}

MC<-function(W,L,P,Response,tol=1e-6)
{
  library(DTMCPack)
  p_0_MC<-0
  p_1_MC<-rep(1/ncol(P),ncol(P))
  path_MC<-NULL
  path_MC<-rbind(path_MC,p_1_MC)
  k_MC<-0
  while(sqrt(sum((p_1_MC-p_0_MC)^2))>tol)
  {
    #Sigma<-matrix(0,nrow = 6,ncol = 6)
    p_0_MC<-p_1_MC
    Sigma<-(p_0_MC*t(W))%*%((Response/(W%*%p_0_MC)[,1]/(P%*%p_0_MC)[,1])*L)
    
    sum_Sigma<-(t(rep(1,ncol(P)))%*%Sigma)[1,]
    diag(Sigma)<-diag(Sigma)-sum_Sigma
    diag(Sigma)<-diag(Sigma)-min(diag(Sigma))
    Sigma<-Sigma/sum(Sigma)*ncol(Sigma)
    
    p_1_MC<-statdistr(t(Sigma))[1,]
    
    path_MC<-rbind(path_MC,p_1_MC)
    k_MC<-k_MC+1
  }
  return(list(Estimator=p_1_MC,Iteration=k_MC,Path=path_MC))
}

MC_approx<-function(W,L,P,Response)
{
  library(DTMCPack)
  p_0_MC<-0
  p_1_MC<-rep(1/ncol(P),ncol(P))
  
  p_0_MC<-p_1_MC
  Sigma<-(p_0_MC*t(W))%*%((Response/(W%*%p_0_MC)[,1]/(P%*%p_0_MC)[,1])*L)
  
  sum_Sigma<-(t(rep(1,ncol(P)))%*%Sigma)[1,]
  diag(Sigma)<-mean(sum_Sigma)-sum_Sigma
  Sigma<-Sigma/sum(Sigma)*ncol(Sigma)
  
  p_1_MC<-statdistr(t(Sigma))[1,]
  return(list(Estimator=p_1_MC))
}

EM<-function(a,b,Delta,tol=1e-6)
{
  p_0_Weaver<-0
  p_1_Weaver<-rep(1/length(a),length(a))
  path_Weaver<-NULL
  path_Weaver<-rbind(path_Weaver,p_1_Weaver)
  k_Weaver<-0
  s=sum(a)+sum(b)
  
  while(sqrt(sum((p_0_Weaver-p_1_Weaver)^2))>tol)
  {
    p_0_Weaver<-p_1_Weaver
    tau<-b/((Delta%*%p_0_Weaver)[,1])
    p_1_Weaver<-(a+((p_0_Weaver*t(Delta))%*%tau)[,1])/s
    p_1_Weaver<-p_1_Weaver/sum(p_1_Weaver)
    k_Weaver<-k_Weaver+1
    path_Weaver<-rbind(path_Weaver,p_1_Weaver)
  }
  return(list(Estimator=p_1_Weaver,Iteration=k_Weaver,Path=path_Weaver))
}

MM<-function(a,b,Delta,tol=1e-6)
{
  p_0_Weaver<-0
  p_1_Weaver<-rep(1/length(a),length(a))
  path_Weaver<-NULL
  path_Weaver<-rbind(path_Weaver,p_1_Weaver)
  k_Weaver<-0
  while(sqrt(sum((p_0_Weaver-p_1_Weaver)^2))>tol)
  {
    p_0_Weaver<-p_1_Weaver
    tau<-b/((Delta%*%p_0_Weaver)[,1])
    p_1_Weaver<-a/((t(Delta)%*%tau)[,1])
    p_1_Weaver<-p_1_Weaver/sum(p_1_Weaver)
    k_Weaver<-k_Weaver+1
    path_Weaver<-rbind(path_Weaver,p_1_Weaver)
  }
  return(list(Estimator=p_1_Weaver,Iteration=k_Weaver,Path=path_Weaver))
}

MMEM<-function(a,b,Delta,tol=1e-6)
{
  p_0_Weaver<-0
  p_1_Weaver<-rep(1/length(a),length(a))
  path_Weaver<-NULL
  path_Weaver<-rbind(path_Weaver,p_1_Weaver)
  k_Weaver<-0
  b1<-b[b>=0]
  b2<-b[b<0]
  Delta1<-Delta[b>=0,]
  Delta2<-Delta[b<0,]
  s<-sum(a)+sum(b)
  
  while(sqrt(sum((p_0_Weaver-p_1_Weaver)^2))>tol)
  {
    p_0_Weaver<-p_1_Weaver
    tau1<-b1/((Delta1%*%p_0_Weaver)[,1])
    tau2<-b2/((Delta2%*%p_0_Weaver)[,1])
    p_1_Weaver<-(a+((p_0_Weaver*t(Delta1))%*%tau1)[,1])/(s-((p_0_Weaver*t(Delta2))%*%tau2)[,1])
    p_1_Weaver<-p_1_Weaver/sum(p_1_Weaver)
    k_Weaver<-k_Weaver+1
    path_Weaver<-rbind(path_Weaver,p_1_Weaver)
  }
  return(list(Estimator=p_1_Weaver,Iteration=k_Weaver,Path=path_Weaver))
}



KL<-function(pfit,pbase)
{
  sum(pfit*(log(pfit/pbase)))
}

Mdist<-function(pfit,pbase,W,P,Response)
{
  WP<-W/((W%*%pfit)[,1])
  PP<-P/((P%*%pfit)[,1])
  I<-t(WP)%*%diag(Response)%*%WP-t(PP)%*%diag(Response)%*%PP
  (t(pfit-pbase)%*%I%*%(pfit-pbase))[1]
}