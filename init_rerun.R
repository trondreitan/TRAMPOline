
init.params.flat2=function(iter=100000,hyp)
{
  par0=init.params.flat(hyp)
  lp0=logprior.flat(par0,hyper)
  ll0=loglik.flat(y,hyp,par0)
  while(!is.numeric(lp0) | !is.numeric(ll0) | is.na(lp0) | is.na(ll0))
  {
    par0=init.params.flat(hyp)
    lp0=logprior.flat(par0,hyper)
    ll0=loglik.flat(y,hyp,par0)
  }
  for(i in 1:iter)
  {
    par1=init.params.flat(hyp)
    lp1=logprior.flat(par1,hyper)
    ll1=loglik.flat(y,hyp,par1)
    if(is.numeric(lp1) & is.numeric(ll1) & !is.na(lp1) & !is.na(ll1))
     if(lp1+ll1>lp0+ll0)
    {
      par0=par1
      ll0=ll1
      lp0=lp1
    }
  }
  return(par0)
}


init.params.flat2.prevmodel=function(iter=100000,hyp,prevmodel.mu,prevmodel.sd)
{
  par0=init.params.flat.prevmodel(hyp,prevmodel.mu,prevmodel.sd)
  lp0=logprior.flat(par0,hyper)
  ll0=loglik.flat(y,hyp,par0)
  while(!is.numeric(lp0) | !is.numeric(ll0) | is.na(lp0) | is.na(ll0))
  {
    par0=init.params.flat.prevmodel(hyp,prevmodel.mu,prevmodel.sd)
    lp0=logprior.flat(par0,hyper)
    ll0=loglik.flat(y,hyp,par0)
  }
  for(i in 1:iter)
  {
    par1=init.params.flat.prevmodel(hyp,prevmodel.mu,prevmodel.sd)
    lp1=logprior.flat(par1,hyper)
    ll1=loglik.flat(y,hyp,par1)
    if(is.numeric(lp1) & is.numeric(ll1) & !is.na(lp1) & !is.na(ll1))
     if(lp1+ll1>lp0+ll0)
    {
      par0=par1
      ll0=ll1
      lp0=lp1
    }
  }
  return(par0)
}

