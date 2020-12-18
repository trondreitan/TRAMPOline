
library(coda)

find.best.par.hpc=function(mcmc.fit,burnin=100,doshow=F)
{
  best.para=NA
  best.lp=best.p=best.l=best.i=best.j=-1e+200
    
  for(j in 1:length(mcmc.fit))
  {
    fit=mcmc.fit[[j]]
    res=mcmc(fit$Monitor[burnin:dim(fit$Monitor)[1],])
    lp=rep(NA,dim(res)[1])
    l=rep(NA,dim(res)[1])
    p=rep(NA,dim(res)[1])
  
    for(i in 1:dim(res)[1])
    {
      l[i]=loglik.flat(y,data,res[i,])
      p[i]=logprior.flat(res[i,],data$hyper)
      lp[i]=l[i]+p[i]
      if(lp[i]>best.lp)
      {
        best.para=res[i,]
        best.lp=lp[i]
        best.p=p[i]
        best.l=l[i]
        best.i=i
	best.j=j
      }
      if(doshow)
        print.srcref(sprintf("%d %d %7.2f %7.2f %7.2f",j, i,
          l[i],p[i],lp[i]))
    }
  }

  return(list(par=best.para, best.i=best.i, best.lp=best.lp))
}

find.best.par=function(mcmc.fit,data,burnin=100,doshow=F)
{
  if(class(mcmc.fit)=="demonoid.hpc")
    return(find.best.par.hpc(mcmc.fit,burnin,doshow))
  
  res=mcmc(mcmc.fit$Monitor[burnin:dim(mcmc.fit$Monitor)[1],])
  lp=rep(NA,dim(res)[1])
  l=rep(NA,dim(res)[1])
  p=rep(NA,dim(res)[1])
  
  best.para=NA
  best.lp=best.p=best.l=best.i=-1e+200
  for(i in 1:dim(res)[1])
  {
    l[i]=loglik.flat(y,data,res[i,])
    p[i]=logprior.flat(res[i,],data$hyper)
    lp[i]=l[i]+p[i]
    if(lp[i]>best.lp)
    {
      best.para=res[i,]
      best.lp=lp[i]
      best.p=p[i]
      best.l=l[i]
      best.i=i
    }
    if(doshow)
      print.srcref(sprintf("%d %7.2f %7.2f %7.2f",i,
        l[i],p[i],lp[i]))
  }

  return(list(par=best.para, best.i=best.i, best.lp=best.lp, lp=lp, l=l, p=p))
}
