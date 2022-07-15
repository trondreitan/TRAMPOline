
# Requires data and betabin definitions to have been read already.

hyper.this=list(logsd.mu=0, logsd.sd=4, log.betabin.s.mu=0, log.betabin.s.sd=4,
                logit.mu=0, logit.sd=4, log.lambda.mu=0, log.lambda.sd=4)

init.params.flat=function(hyp)
{
  N.sp=hyp$N.sp
  N=hyp$N
  N.forms=hyp$N.forms
  Nf=hyp$Nf

  ret=list(
   log.betabin.s=rnorm(N.sp, hyper$log.betabin.s.mu, hyper$log.betabin.s.sd/4),
   log.const.lambda=rnorm(N.sp,hyper$log.lambda.mu,hyper$log.lambda.sd),
   logit.const.occu=rnorm(N.sp-1,hyper$logit.mu,hyper$logit.sd),
   log.sd.f.occu=rnorm(1,hyper$logsd.mu, hyper$logsd.sd/4),
   log.sd.f.lambda=rnorm(1,hyper$logsd.mu, hyper$logsd.sd/4),
   log.sd.sf.occu=rnorm(N.sp-1,hyper$logsd.mu-2, hyper$logsd.sd/4),
   log.sd.sf.lambda=rnorm(N.sp,hyper$logsd.mu-2, hyper$logsd.sd/4),
   f.occu=rep(0,Nf),
   f.lambda=rep(0,Nf),
   sf.occu=matrix(0, N.forms, N.sp-1),
   sf.lambda=matrix(0, N.forms, N.sp)
   )
 
 ret$f.occu=rnorm(Nf, 0, exp(ret$log.sd.f.occu)/4)
 ret$f.lambda=rnorm(Nf, 0, exp(ret$log.sd.f.lambda)/4)
 
 for(i in 1:(N.sp-1))
  ret$sf.occu[,i]=rnorm(Nf, 0, exp(ret$log.sd.sf.occu[i])/4)
    
 for(i in 1:N.sp)
  ret$sf.lambda[,i]=rnorm(Nf, 0, exp(ret$log.sd.sf.lambda[i])/4)
    
 return(unlist(ret))
}




init.params.flat.prevmodel=function(hyp,prevmodel.mu,prevmodel.sd)
{
  N.sp=hyp$N.sp
  N=hyp$N
  N.forms=hyp$N.forms
  Nf=hyp$Nf

  ret=list(
   log.betabin.s=rnorm(N.sp, prevmodel.mu[1:N.sp], prevmodel.sd[1:N.sp]), 
   log.const.lambda=rnorm(N.sp, prevmodel.mu[N.sp+1:N.sp], prevmodel.sd[N.sp+1:N.sp]),
   logit.const.occu=rnorm(N.sp-1,prevmodel.mu[2*N.sp+1:(N.sp-1)],prevmodel.sd[2*N.sp+1:(N.sp-1)]),
   log.sd.f.occu=rnorm(1,prevmodel.mu[3*N.sp], prevmodel.sd[3*N.sp]),
   log.sd.f.lambda=rnorm(1,prevmodel.mu[3*N.sp+1], prevmodel.sd[3*N.sp+1]),
   
   log.sd.sf.occu=rnorm(N.sp-1,hyper$logsd.mu-2, hyper$logsd.sd/4),
   log.sd.sf.lambda=rnorm(N.sp,hyper$logsd.mu-2, hyper$logsd.sd/4),
   f.occu=rnorm(Nf,prevmodel.mu[3*N.sp+1+1:Nf],prevmodel.sd[3*N.sp+1+1:Nf]),
   f.lambda=rnorm(Nf,prevmodel.mu[Nf+3*N.sp+1+1:Nf],prevmodel.sd[Nf+3*N.sp+1+1:Nf]),
   sf.occu=matrix(0, N.forms, N.sp-1),
   sf.lambda=matrix(0, N.forms, N.sp)
   )
 
 for(i in 1:(N.sp-1))
  ret$sf.occu[,i]=rnorm(Nf, 0, exp(ret$log.sd.sf.occu[i])/4)
    
 for(i in 1:N.sp)
  ret$sf.lambda[,i]=rnorm(Nf, 0, exp(ret$log.sd.sf.lambda[i])/4)
 
 return(unlist(ret))
}


init.params.flat.prevmodel2=function(hyp,prevmodel.mu,prevmodel.sd)
{
  N.sp=hyp$N.sp
  N=hyp$N
  N.forms=hyp$N.forms
  Nf=hyp$Nf

  ret=list(
   log.betabin.s=prevmodel.mu[1:N.sp] ,
   log.const.lambda=prevmodel.mu[N.sp+1:N.sp],
   logit.const.occu=prevmodel.mu[2*N.sp+1:(N.sp-1)],
   log.sd.f.occu=prevmodel.mu[3*N.sp],
   log.sd.f.lambda=prevmodel.mu[3*N.sp+1],
   log.sd.sf.occu=prevmodel.mu[(3*N.sp+2):(4*N.sp)],
   log.sd.sf.lambda=rnorm(N.sp,hyper$logsd.mu, hyper$logsd.sd/4),
   f.occu=prevmodel.mu[4*N.sp+1:Nf],
   f.lambda=prevmodel.mu[Nf+4*N.sp+1:Nf],
   sf.occu=matrix(prevmodel.mu[2*Nf+4*N.sp+1:(Nf*(N.sp-1))], N.forms, N.sp-1),
   sf.lambda=matrix(0, N.forms, N.sp)
   )
 
 for(i in 1:N.sp)
 { #show(i)
  ret$sf.lambda[,i]=rnorm(Nf, 0, 0.001)
 }
 
 return(unlist(ret))
}


logprior.flat=function(par, hyp)
{
  N.sp=hyp$N.sp
  N=hyp$N
  N.forms=hyp$N.forms
  Nf=hyp$Nf

  i=0
  
  # Get top parameters:
  log.betabin.s=par[(i+1):(i+N.sp)]; i=i+N.sp
  log.const.lambda=par[(i+1):(i+N.sp)]; i=i+N.sp
  logit.const.occu=par[(i+1):(i+N.sp-1)]; i=i+N.sp-1
  log.sd.f.occu=par[i+1]; i=i+1
  log.sd.f.lambda=par[i+1]; i=i+1
  log.sd.sf.occu=par[(i+1):(i+(N.sp-1))]; i=i+N.sp-1
  log.sd.sf.lambda=par[(i+1):(i+N.sp)]; i=i+N.sp

  if(sum(log.betabin.s< -20)>0)
    return(-1e+200)

  sum(dnorm(log.betabin.s, hyper$log.betabin.s.mu, hyper$log.betabin.s.sd, log=T)) +
  sum(dnorm(log.const.lambda, hyper$log.lambda.mu, hyper$log.lambda.sd, log=T)) +
  sum(dnorm(logit.const.occu, hyper$logit.mu, hyper$logit.sd, log=T)) +
  sum(dnorm(log.sd.f.occu, hyper$logsd.mu, hyper$logsd.sd, log=T)) +
  sum(dnorm(log.sd.f.lambda, hyper$logsd.mu, hyper$logsd.sd, log=T)) +
  sum(dnorm(log.sd.sf.occu, hyper$logsd.mu, hyper$logsd.sd, log=T))  +
  sum(dnorm(log.sd.sf.lambda, hyper$logsd.mu, hyper$logsd.sd, log=T))
}

loglik.flat=function(Y,X,par)
# PS: we're going to take data Y and X right from new.data anyhow, no
# need to send it, really...
{
  N.sp=X$N.sp
  N=X$N
  N.forms=X$N.forms
  Nf=X$Nf
  counts=Y[,(which(names(Y)=="num.with.colonies")+1):which(names(Y)=="Superspecies")]

  form=X$form
  N.rep=X$N.rep

  i=0
  
  # Get top parameters:
  log.betabin.s=par[(i+1):(i+N.sp)]; i=i+N.sp
  betabin.s=exp(log.betabin.s)
  log.const.lambda=par[(i+1):(i+N.sp)]; i=i+N.sp
  logit.const.occu=par[(i+1):(i+N.sp-1)]; i=i+N.sp-1
  log.sd.f.occu=par[i+1]; i=i+1
  sd.f.occu=exp(log.sd.f.occu)
  log.sd.f.lambda=par[i+1]; i=i+1
  sd.f.lambda=exp(log.sd.f.lambda)
  log.sd.sf.occu=par[(i+1):(i+N.sp-1)]; i=i+N.sp-1
  sd.sf.occu=exp(log.sd.sf.occu)
  log.sd.sf.lambda=par[(i+1):(i+N.sp)]; i=i+N.sp
  sd.sf.lambda=exp(log.sd.sf.lambda)
  
  # Get random variables ...
  f.occu=as.numeric(par[(i+1):(i+N.forms)]); i=i+N.forms
  f.lambda=as.numeric(par[(i+1):(i+N.forms)]); i=i+N.forms

  sf.occu=matrix(par[(i+1):(i+Nf*(N.sp-1))], ncol=N.sp-1); i=i+Nf*(N.sp-1)
  sf.lambda=matrix(par[(i+1):(i+Nf*N.sp)], ncol=N.sp); i=i+Nf*N.sp


  # Start likelihood calculations:
  ll=0
  betabin.s.use=rep(1,N)%*%t(betabin.s)
  
  ll=ll+sum(dnorm(f.occu, 0, sd.f.occu, log=T))
  f.occu.use=f.occu%*%t(rep(1,N.sp-1))
  
  ll=ll+sum(dnorm(f.lambda, 0, sd.f.lambda, log=T))
  f.lambda.use=f.lambda%*%t(rep(1,N.sp))
  
  sd.sf.occu.use=rep(1,Nf)%*%t(sd.sf.occu)
  ll=ll+sum(dnorm(sf.occu, 0, sd.sf.occu.use,log=T))
  
  sd.sf.lambda.use=rep(1,Nf)%*%t(sd.sf.lambda)
  ll=ll+sum(dnorm(sf.lambda, 0, sd.sf.lambda.use,log=T))
  
  const.occu.mat=rep(1,N)%*%t(logit.const.occu)
  occu.mat=array(1, c(N, N.sp))
  occu.mat[,(1:(N.sp-1))]=ilogit(const.occu.mat+sf.occu[Y$form.nr,] +
                                 f.occu.use[Y$form.nr,] )

  const.lambda.mat=rep(1,N)%*%t(log.const.lambda)
  lambda.mat=exp(const.lambda.mat + f.lambda.use[Y$form.nr,] + sf.lambda[Y$form.nr,] )
    
  if(sum(lambda.mat==0 & counts>0)>0)
    return(-1e+200)
  
  nrep=Y$Total%*%t(rep(1,N.sp))
  ll=ll+sum(dlbetabin.zero(counts, nrep, 1-exp(-lambda.mat), betabin.s.use, occu.mat))
  
  return(ll)
}

