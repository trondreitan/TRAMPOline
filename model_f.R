
# Requires data and betabin definitions to have been read already.

hyper.this=list(logsd.mu=0, logsd.sd=4, log.betabin.s.mu=0, log.betabin.s.sd=4,
                logit.mu=0, logit.sd=4)

init.params.flat=function(hyp)
{
  N.sp=hyp$N.sp
  N=hyp$N
  N.forms=hyp$N.forms
  Nf=hyp$Nf

  ret=list(
   log.betabin.s=rnorm(N.sp, hyper$log.betabin.s.mu, hyper$log.betabin.s.sd/4),
   logit.const.bin=rnorm(N.sp,hyper$logit.mu,hyper$logit.sd),
   logit.const.occu=rnorm(N.sp-1,hyper$logit.mu,hyper$logit.sd),
   log.sd.f.occu=rnorm(1,hyper$logsd.mu, hyper$logsd.sd/4),
   log.sd.f.bin=rnorm(1,hyper$logsd.mu, hyper$logsd.sd/4),
   f.occu=rep(0,Nf),
   f.bin=rep(0,Nf)
   )
 
 ret$f.occu=rnorm(Nf, 0, exp(ret$log.sd.f.occu))
 ret$f.bin=rnorm(Nf, 0, exp(ret$log.sd.f.bin))
 
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
   logit.const.bin=rnorm(N.sp, prevmodel.mu[N.sp+1:N.sp], prevmodel.sd[N.sp+1:N.sp]),
   logit.const.occu=rnorm(N.sp-1,prevmodel.mu[2*N.sp+1:(N.sp-1)],prevmodel.sd[2*N.sp+1:(N.sp-1)]),
   log.sd.f.occu=rnorm(1,prevmodel.mu[3*N.sp], prevmodel.sd[3*N.sp]),
   log.sd.f.bin=rnorm(1,hyper$logsd.mu, hyper$logsd.sd/4),
   f.occu=rnorm(Nf,prevmodel.mu[3*N.sp+1:Nf],prevmodel.sd[3*N.sp+1:Nf]),
   f.bin=rep(0,Nf)
   )
 
 ret$f.bin=rnorm(Nf, 0, exp(ret$log.sd.f.bin))

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
  logit.const.bin=par[(i+1):(i+N.sp)]; i=i+N.sp
  logit.const.occu=par[(i+1):(i+N.sp-1)]; i=i+N.sp-1
  log.sd.f.occu=par[i+1]; i=i+1
  log.sd.f.bin=par[i+1]; i=i+1
  
  sum(dnorm(log.betabin.s, hyper$log.betabin.s.mu, hyper$log.betabin.s.sd, log=T)) +
  sum(dnorm(logit.const.bin, hyper$logit.mu, hyper$logit.sd, log=T)) +
  sum(dnorm(logit.const.occu, hyper$logit.mu, hyper$logit.sd, log=T)) +
  sum(dnorm(log.sd.f.occu, hyper$logsd.mu, hyper$logsd.sd, log=T)) +
  sum(dnorm(log.sd.f.bin, hyper$logsd.mu, hyper$logsd.sd, log=T))
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


  i=0
  
  # Get top parameters:
  log.betabin.s=par[(i+1):(i+N.sp)]; i=i+N.sp
  betabin.s=exp(log.betabin.s)
  logit.const.bin=par[(i+1):(i+N.sp)]; i=i+N.sp
  logit.const.occu=par[(i+1):(i+N.sp-1)]; i=i+N.sp-1
  log.sd.f.occu=par[i+1]; i=i+1
  sd.f.occu=exp(log.sd.f.occu)
  log.sd.f.bin=par[i+1]; i=i+1
  sd.f.bin=exp(log.sd.f.bin)

  # Get random variables ...
  f.occu=as.numeric(par[(i+1):(i+N.forms)]); i=i+N.forms
  f.bin=as.numeric(par[(i+1):(i+N.forms)]); i=i+N.forms

  
  ll=0
  betabin.s.use=rep(1,N)%*%t(betabin.s)
  
  ll=ll+sum(dnorm(f.occu, 0, sd.f.occu, log=T))
  f.occu.use=f.occu%*%t(rep(1,N.sp-1))
  
  ll=ll+sum(dnorm(f.bin, 0, sd.f.bin, log=T))
  f.bin.use=f.bin%*%t(rep(1,N.sp))
  
  const.occu.mat=rep(1,N)%*%t(logit.const.occu)
  occu.mat=array(1, c(N, N.sp))
  occu.mat[,(1:(N.sp-1))]=ilogit(const.occu.mat+f.occu.use[Y$form.nr,] )

  const.bin.mat=rep(1,N)%*%t(logit.const.bin)
  bin.mat=ilogit(const.bin.mat + f.bin.use[Y$form.nr,] )
    
  if(sum(bin.mat==0 & counts>0)>0)
    return(-1e+200)
  
  nrep=Y$Total%*%t(rep(1,4))
  ll=ll+sum(dlbetabin.zero(counts, nrep, bin.mat, betabin.s.use, occu.mat))

  return(ll)
}

