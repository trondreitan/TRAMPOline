
prevmodel="foccu"
model="f"

N.mcmc=20000
run.nr=4
thin=20

source("3stage.R")



# Collate (if multiple runs)
N.runs=4

ll=rep(0,8*N.runs)
# all
for(i in 1:N.runs)
{
  load(sprintf("%s_%02d.Rdata",model,i))
  for(j in 1:8)
   ll[j+(i-1)*8]=Fit3[[j]]$LML
  if(i==1)
    res=mcmc.list(mcmc(Fit3[[1]]$Monitor),
              mcmc(Fit3[[2]]$Monitor),
              mcmc(Fit3[[3]]$Monitor),
              mcmc(Fit3[[4]]$Monitor),
              mcmc(Fit3[[5]]$Monitor),
              mcmc(Fit3[[6]]$Monitor),
              mcmc(Fit3[[7]]$Monitor),
              mcmc(Fit3[[8]]$Monitor))
  k=length(res)
  if(i>1)
   for(j in 1:8)
    res[[j+k]]=mcmc(Fit3[[j]]$Monitor)
}
ll=ll[!is.na(ll)]
ll=ll[ll!=0]
# Nothing! :(


s=summary(res)
par.mean=s$statistics[,1]
save(par.mean,file=sprintf("par_mean_%s.RData",model))
par.sd=s$statistics[,2]
save(par.sd,file=sprintf("par_sd_%s.RData",model))


res2=as.matrix(res)
par.median=par0
for(i in 1:length(par0))
  par.median[i]=median(res2[,i])
save(par.median,file=sprintf("par_median_%s.RData",model))



gelman.diag(res)
# 1.7
# :/

plot(res)






source("http://folk.uio.no/trondr/R/mcmc.R")

ll2=rep(0,32)
for(i in 1:32)
{
  show(i)
  ll2[i]=model.loglik(loglik.flat, logprior.flat, y, hyper, hyper, par0, t(res2), 1000)
}
t.test(ll2)
# -1001.833 -1001.612
#-1001.723 

sd(ll2)
# 0.3064665
 






