# Beta binmial and zero-inflated beta-binomial distribution function,
# original scale and log-transformed.

# Zero-inflated binomial:

dbinom.zero=function(x,n,p,p.nonzero)
  p.nonzero*dbinom(x,n,p,p.nonzero)+(1-p.nonzero)*(x==0)

dlbinom.zero=function(x,n,p,p.nonzero)
  log(p.nonzero*dbinom(x,n,p,p.nonzero)+(1-p.nonzero)*(x==0))


# Beta-binomial:

dbetabin=function(x,n,p,s)
  exp(lgamma(n+1)-lgamma(x+1)-lgamma(n-x+1) +
      lgamma(x+p/s)+lgamma(n-x+(1-p)/s)-lgamma(n+1/s) +
      lgamma(1/s)-lgamma(p/s)-lgamma((1-p)/s))
    
dlbetabin=function(x,n,p,s)
  lgamma(n+1)-lgamma(x+1)-lgamma(n-x+1) +
  lgamma(x+p/s)+lgamma(n-x+(1-p)/s)-lgamma(n+1/s) +
  lgamma(1/s)-lgamma(p/s)-lgamma((1-p)/s)


# Zero-inflated betea-binomial:

dbetabin.zero=function(x,n,p,s,p.nonzero)
  p.nonzero*dbetabin(x,n,p,s)+(1-p.nonzero)*(x==0)
      
dlbetabin.zero=function(x,n,p,s,p.nonzero)
  log(p.nonzero*dbetabin(x,n,p,s)+(1-p.nonzero)*(x==0))

