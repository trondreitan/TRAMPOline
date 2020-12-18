
# Requires data and betabin definitions to have been read already.

hyper.this=list(logsd.mu=0, logsd.sd=4,
                log.betabin.s.mu=0, log.betabin.s.sd=4,
                logit.mu=0, logit.sd=4,
		o18.optim.mu=4.0, o18.optim.sd=2,
		o18.lwidth.mu=0, o18.lwidth.sd=1)

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
   log.sd.sf.occu=rnorm(N.sp-1,hyper$logsd.mu-2, hyper$logsd.sd/4),
   log.sd.sf.bin=rnorm(N.sp,hyper$logsd.mu-2, hyper$logsd.sd/4),
   o18.optim.bin=rnorm(N.sp,hyper$o18.optim.mu,hyper$o18.optim.sd/4),
   o18.lwidth.bin=rnorm(N.sp,hyper$o18.lwidth.mu+7,hyper$o18.lwidth.sd/4),
   f.occu=rep(0,Nf),
   f.bin=rep(0,Nf),
   sf.occu=matrix(0, N.forms, N.sp-1),
   sf.bin=matrix(0, N.forms, N.sp)
   )
 
 ret$f.occu=rnorm(Nf, 0, exp(ret$log.sd.f.occu)/4)
 ret$f.bin=rnorm(Nf, 0, exp(ret$log.sd.f.bin)/4)
 
 for(i in 1:(N.sp-1))
  ret$sf.occu[,i]=rnorm(Nf, 0, exp(ret$log.sd.sf.occu[i])/4)
    
 for(i in 1:N.sp)
  ret$sf.bin[,i]=rnorm(Nf, 0, exp(ret$log.sd.sf.bin[i])/4)
    
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
   log.sd.f.bin=rnorm(1,prevmodel.mu[3*N.sp+1], prevmodel.sd[3*N.sp+1]),
   log.sd.sf.occu=prevmodel.mu[(3*N.sp+2):(4*N.sp)],
   log.sd.sf.bin=rnorm(N.sp,hyper$logsd.mu, hyper$logsd.sd/4),
   o18.optim.bin=rnorm(N.sp,hyper$o18.optim.mu,hyper$o18.optim.sd/4),
   o18.lwidth.bin=rnorm(N.sp,hyper$o18.lwidth.mu+7,hyper$o18.lwidth.sd/4),
   f.occu=rnorm(Nf,prevmodel.mu[3*N.sp+1+1:Nf],prevmodel.sd[3*N.sp+1+1:Nf]),
   f.bin=rnorm(Nf,prevmodel.mu[Nf+3*N.sp+1+1:Nf],prevmodel.sd[Nf+3*N.sp+1+1:Nf]),
   sf.occu=matrix(0, N.forms, N.sp-1),
   sf.bin=matrix(0, N.forms, N.sp)
   )
 
 for(i in 1:(N.sp-1))
  ret$sf.occu[,i]=rnorm(Nf, 0, exp(ret$log.sd.sf.occu[i])/4)
    
 for(i in 1:N.sp)
  ret$sf.bin[,i]=rnorm(Nf, 0, exp(ret$log.sd.sf.bin[i])/4)
 
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
   logit.const.bin=prevmodel.mu[N.sp+1:N.sp],
   logit.const.occu=prevmodel.mu[2*N.sp+1:(N.sp-1)],
   log.sd.f.occu=prevmodel.mu[3*N.sp],
   log.sd.f.bin=prevmodel.mu[3*N.sp+1],
   log.sd.sf.occu=prevmodel.mu[(3*N.sp+2):(4*N.sp)],
   log.sd.sf.bin=rnorm(N.sp,hyper$logsd.mu, hyper$logsd.sd/4),
   o18.optim.bin=rnorm(N.sp,hyper$o18.optim.mu,hyper$o18.optim.sd/4),
   o18.lwidth.bin=rnorm(N.sp,hyper$o18.lwidth.mu+7,hyper$o18.lwidth.sd/4),
   f.occu=prevmodel.mu[4*N.sp+1:Nf],
   f.bin=prevmodel.mu[Nf+4*N.sp+1:Nf],
   sf.occu=matrix(prevmodel.mu[2*Nf+4*N.sp+1:(Nf*(N.sp-1))], N.forms, N.sp-1),
   sf.bin=matrix(0, N.forms, N.sp)
   )
 
 for(i in 1:N.sp)
 { #show(i)
  ret$sf.bin[,i]=rnorm(Nf, 0, 0.001)
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
  logit.const.bin=par[(i+1):(i+N.sp)]; i=i+N.sp
  logit.const.occu=par[(i+1):(i+N.sp-1)]; i=i+N.sp-1
  log.sd.f.occu=par[i+1]; i=i+1
  log.sd.f.bin=par[i+1]; i=i+1
  log.sd.sf.occu=par[(i+1):(i+(N.sp-1))]; i=i+N.sp-1
  log.sd.sf.bin=par[(i+1):(i+N.sp)]; i=i+N.sp
  o18.optim.bin=par[(i+1):(i+N.sp)]; i=i+N.sp
  o18.lwidth.bin=par[(i+1):(i+N.sp)]; i=i+N.sp

  if(sum(log.betabin.s< -20)>0)
    return(-1e+200)

  sum(dnorm(log.betabin.s, hyper$log.betabin.s.mu, hyper$log.betabin.s.sd, log=T)) +
  sum(dnorm(logit.const.bin, hyper$logit.mu, hyper$logit.sd, log=T)) +
  sum(dnorm(logit.const.occu, hyper$logit.mu, hyper$logit.sd, log=T)) +
  sum(dnorm(log.sd.f.occu, hyper$logsd.mu, hyper$logsd.sd, log=T)) +
  sum(dnorm(log.sd.f.bin, hyper$logsd.mu, hyper$logsd.sd, log=T)) +
  sum(dnorm(log.sd.sf.occu, hyper$logsd.mu, hyper$logsd.sd, log=T))  +
  sum(dnorm(log.sd.sf.bin, hyper$logsd.mu, hyper$logsd.sd, log=T))  +
  sum(dnorm(o18.optim.bin, hyper$o18.optim.mu, hyper$o18.optim.sd, log=T)) +
  sum(dnorm(o18.lwidth.bin, hyper$o18.lwidth.mu, hyper$o18.lwidth.sd, log=T))
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
  logit.const.bin=par[(i+1):(i+N.sp)]; i=i+N.sp
  logit.const.occu=par[(i+1):(i+N.sp-1)]; i=i+N.sp-1
  log.sd.f.occu=par[i+1]; i=i+1
  sd.f.occu=exp(log.sd.f.occu)
  log.sd.f.bin=par[i+1]; i=i+1
  sd.f.bin=exp(log.sd.f.bin)
  log.sd.sf.occu=par[(i+1):(i+N.sp-1)]; i=i+N.sp-1
  sd.sf.occu=exp(log.sd.sf.occu)
  log.sd.sf.bin=par[(i+1):(i+N.sp)]; i=i+N.sp
  sd.sf.bin=exp(log.sd.sf.bin)
  o18.optim.bin=par[(i+1):(i+N.sp)]; i=i+N.sp
  o18.lwidth.bin=par[(i+1):(i+N.sp)]; i=i+N.sp
  o18.width.bin=exp(o18.lwidth.bin)
  
  # Get random variables ...
  f.occu=as.numeric(par[(i+1):(i+N.forms)]); i=i+N.forms
  f.bin=as.numeric(par[(i+1):(i+N.forms)]); i=i+N.forms

  sf.occu=matrix(par[(i+1):(i+Nf*(N.sp-1))], ncol=N.sp-1); i=i+Nf*(N.sp-1)
  sf.bin=matrix(par[(i+1):(i+Nf*N.sp)], ncol=N.sp); i=i+Nf*N.sp


  # Start likelihood calculations:
  ll=0
  betabin.s.use=rep(1,N)%*%t(betabin.s)
  
  ll=ll+sum(dnorm(f.occu, 0, sd.f.occu, log=T))
  f.occu.use=f.occu%*%t(rep(1,N.sp-1))
  
  ll=ll+sum(dnorm(f.bin, 0, sd.f.bin, log=T))
  f.bin.use=f.bin%*%t(rep(1,N.sp))
  
  sd.sf.occu.use=rep(1,Nf)%*%t(sd.sf.occu)
  ll=ll+sum(dnorm(sf.occu, 0, sd.sf.occu.use,log=T))
  
  sd.sf.bin.use=rep(1,Nf)%*%t(sd.sf.bin)
  ll=ll+sum(dnorm(sf.bin, 0, sd.sf.bin.use,log=T))
  
  const.occu.mat=rep(1,N)%*%t(logit.const.occu)
  occu.mat=array(1, c(N, N.sp))
  occu.mat[,(1:(N.sp-1))]=ilogit(const.occu.mat+sf.occu[Y$form.nr,] +
                                 f.occu.use[Y$form.nr,] )

  const.bin.mat=rep(1,N)%*%t(logit.const.bin)
  bin.mat=ilogit(const.bin.mat + f.bin.use[Y$form.nr,] + sf.bin[Y$form.nr,] -
    (Y$O18.mean%*%t(rep(1,N.sp))-rep(1,N)%*%t(o18.optim.bin))^2/
	(rep(1,N)%*%t(o18.width.bin))^2 )
	
  if(sum((occu.mat==0 | bin.mat==0) & counts>0)>0)
    return(-1e+200)
  
  nrep=Y$Total%*%t(rep(1,4))
  ll=ll+sum(dlbetabin.zero(counts, nrep, bin.mat, betabin.s.use, occu.mat))

  return(ll)
}

