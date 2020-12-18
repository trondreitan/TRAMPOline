
prevmodel=""
model="global"

N.mcmc=20000
run.nr=4
thin=20

source("3stage.R")



# Collate (if multiple runs)
N.runs=4
runs=c(1,2,3,4)

ll=rep(0,8*N.runs)
# all
for(i in 1:N.runs)
{
  load(sprintf("%s_%02d.Rdata",model,runs[i]))
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

t.test(ll)
# -1105.354 -1104.463
# -1104.909 

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
# Multivariate psrf 1.01

plot(res)

s
#                     2.5%      25%      50%     75%  97.5%
#log.betabin.s1    -1.7817 -1.49923 -1.35398 -1.2058 -0.924
#log.betabin.s2    -2.7449 -2.33095 -2.09955 -1.8598 -1.294
#log.betabin.s3    -2.8504 -2.28467 -1.96120 -1.6557 -1.210
#log.betabin.s4    -2.0142 -1.83255 -1.73577 -1.6378 -1.449
#logit.const.bin1  -2.6366 -2.42796 -2.31261 -2.1914 -1.933
#logit.const.bin2  -2.2812 -1.86946 -1.74506 -1.6400 -1.455
#logit.const.bin3  -3.0823 -2.69135 -2.38792 -2.1703 -1.900
#logit.const.bin4  -0.4055 -0.30829 -0.25568 -0.2042 -0.106
#logit.const.occu1  1.1342  2.67336  4.25350  6.4289 11.661
#logit.const.occu2 -0.3445 -0.06474  0.09565  0.2787  1.275
#logit.const.occu3 -0.2595  0.16549  0.58087  1.9499  8.835






# subset:
index=9:11
# all
for(i in 1:N.runs)
{
  load(sprintf("%s_%02d.Rdata",model,i))
  if(i==1)
    res=mcmc.list(mcmc(Fit3[[1]]$Monitor[,index]),
              mcmc(Fit3[[2]]$Monitor[,index]),
              mcmc(Fit3[[3]]$Monitor[,index]),
              mcmc(Fit3[[4]]$Monitor[,index]),
              mcmc(Fit3[[5]]$Monitor[,index]),
              mcmc(Fit3[[6]]$Monitor[,index]),
              mcmc(Fit3[[7]]$Monitor[,index]),
              mcmc(Fit3[[8]]$Monitor[,index]))
  k=length(res)
  if(i>1)
   for(j in 1:8)
    res[[j+k]]=mcmc(Fit3[[j]]$Monitor[,index])
}
plot(res)





source("http://folk.uio.no/trondr/R/mcmc.R")

ll2=rep(0,32)
for(i in 1:32)
{
  show(i)
  ll2[i]=model.loglik(loglik.flat, logprior.flat, y, hyper, hyper, par0, t(res2), 1000)
}
t.test(ll2)
# -1137.724 -1137.622
# -1137.673 

