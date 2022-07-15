
# Requires data and betabin definitions to have been read already.

hyper.this=list(log.betabin.s.mu=0, log.betabin.s.sd=5,
                logit.mu=0, logit.sd=5)

init.params.flat=function(hyp)
{
  N.sp=hyp$N.sp
  N=hyp$N
  N.forms=hyp$N.forms
  Nf=hyp$Nf

  ret=list(
   log.betabin.s=rnorm(N.sp, hyper$log.betabin.s.mu, hyper$log.betabin.s.sd/4),
   logit.const.bin=rnorm(N.sp,hyper$logit.mu,hyper$logit.sd),
   logit.const.occu=rnorm(N.sp-1,hyper$logit.mu,hyper$logit.sd)
   )

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
  
  sum(dnorm(log.betabin.s, hyper$log.betabin.s.mu, hyper$log.betabin.s.sd, log=T)) +
  sum(dnorm(logit.const.bin, hyper$logit.mu, hyper$logit.sd, log=T)) +
  sum(dnorm(logit.const.occu, hyper$logit.mu, hyper$logit.sd, log=T))
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

  # Get random variables ... none
  
  ll=0
  betabin.s.use=rep(1,N)%*%t(betabin.s)

  const.occu.mat=rep(1,N)%*%t(logit.const.occu)
  occu.mat=array(1, c(N, N.sp))
  occu.mat[,(1:(N.sp-1))]=ilogit(const.occu.mat)
  

  const.bin.mat=rep(1,N)%*%t(logit.const.bin)
  bin.mat=ilogit(const.bin.mat )
    
  if(sum(bin.mat==0 & counts>0)>0)
    return(-1e+200)
  
  nrep=Y$Total%*%t(rep(1,N.sp))
  ll=ll+sum(dlbetabin.zero(counts, nrep, bin.mat, betabin.s.use, occu.mat))

  return(ll)
}

