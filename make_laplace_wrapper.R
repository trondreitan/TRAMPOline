
library("LaplacesDemon")


Model<-function(parm, Data)
{
  hyper=Data$hyper
  beta.prior=logprior.flat(parm,hyper)
  LL=loglik.flat(Data$y,Data,parm)
  LP=LL+beta.prior
  Model.out=list(LP=LP, Dev=-2*LL, Monitor=parm,
    yhat=rep(0,length(Data$y)),
    parm=parm)
  return(Model.out)
}

make.data=function(new.par0,hyp,y)
{
  mon.names=names(new.par0)
  parm.names=names(new.par0)
  PGF=function(Data)
  {
    beta=init.params.flat2(10000,Data$hyper)
    return(beta)
  }

  data=list(PGF=PGF,x=NA,y=y,mon.names=mon.names,parm.names=parm.names,
     hyper=hyp,N=hyp$N,Nf=hyp$Nf,N.forms=hyp$N.forms,N.sp=hyp$N.sp,
     forms=hyp$forms, K=hyp$K, K1=hyp$K1, K2=hyp$K2)
  
  return(data)
}

