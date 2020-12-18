# Look at a single model:

model="f_sf"

source("read_data.R")
source("betabin.R")
source(sprintf("model_%s.R",model))
hyper=hyper.this
hyper$N=N
hyper$N.sp=N.sp
hyper$Nf=Nf
hyper$N.forms=N.forms
source("make_laplace_wrapper.R")
library(coda)

par0=init.params.flat(hyper)

# Collate (if multiple runs)
N.start=1
N.runs=16

ll=rep(NA,100*N.runs)
# all
f=list()
for(i in 1:N.runs)
{
  load(sprintf("%s_%02d.Rdata",model,i+N.start-1))
  start=1
  if(i==1)
  {
    res=mcmc.list(mcmc(Fit3[[1]]$Monitor))
    start=2
    f=Fit3
  }
  if(i>1)
  {
    for(j in 1:length(Fit3))
      {
        k=length(f)
	show(c(i,j,k))
	f[[k+1]]=Fit3[[j]]
      }
  }
  cc=length(Fit3)
  for(j in 1:cc)
   ll[j+(i-1)*cc]=Fit3[[j]]$LML
  k=length(res)
  for(j in start:cc)
    res[[j+k-start+1]]=mcmc(Fit3[[j]]$Monitor)
  
}
ll=ll[!is.na(ll)]
ll=ll[ll!=0]; ll

t.test(ll)
# f_sf 1-16:
# -722.0639 -678.4025 -700.2332
# fdet_sf 1-8:
# -758.8911 -721.3820 -740.1365 
# f_sf_o18bin 1-9:
# -697.0253 -627.0155 -662.0204
# f_sf_o18bin 11-20:
# -701.9372 -657.3354 -679.6363 
# f_sf_o18bin 21-40:
# -715.5151 -680.4033 -697.9592 
# f_sf_ou 1-40:
# -714.6724 -672.5636 -693.618 


source("make_laplace_wrapper.R")
data=make.data(par0,hyper,y)
d=Combine(f,Data=data)
d$LML
# global 1-4:
# -1107
# foccu 1-4:
# -1028.8
# f 1-4:
# -942.5
# f_sf 1-16:
# -678.2954
# fdet_sfdet 1-8:
# -807.6239
# fdet_sf 1-8:
# -754.8282
# f_sf_o18bin 1-9:
# -700.5181
# f_sf_o18bin 1-10:
# -638.866
# f_sf_o18bin 11-20:
# -677.5055
# f_sf_o18bin 21-30:
# -726.2049
# f_sf_o18bin 31-40:
# -674.57
# f_sf_o18bin 11-14:
# -692.0427
# f_sf_ou: 1-20:
# -578.1121
# f_sf_ou: 21-40:
# -663.8038
# f_sf_ou: 51-56
# -677.3947

par.mean=as.numeric(substr(summary(d$Monitor)[4,],9,nchar(summary(d$Monitor)[4,])))
names(par.mean)=names(par0)

par.median=as.numeric(substr(summary(d$Monitor)[4,],9,nchar(summary(d$Monitor)[3,])))
names(par.median)=names(par0)

par.sd=(as.numeric(substr(summary(d$Monitor)[4,],9,nchar(summary(d$Monitor)[3,])))-
 as.numeric(substr(summary(d$Monitor)[5,],9,nchar(summary(d$Monitor)[2,]))))/
 (qnorm(0.75)-qnorm(0.25))
names(par.sd)=names(par0)

save(par.mean,file=sprintf("par_mean_%s.RData",model))
save(par.sd,file=sprintf("par_sd_%s.RData",model))
save(par.median,file=sprintf("par_median_%s.RData",model))




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
# f_sf_o18bin 1-9:
# 1.18
# f_sf_o18bin 1-10:
# 1.42

plot(res)







# Lookat subset of parameters:

# 4 worst parameters:
g=gelman.diag(res)
index=which(g$psrf[,1]>=sort(g$psrf[,1],decreasing=TRUE)[4])
#index=which(attributes(d$Monitor)$dimnames[[2]]=="log.sd.sf.bin1")
for(i in 1:N.runs)
{
  load(sprintf("%s_%02d.Rdata",model,i+N.start-1))
  start=1
  if(i==1)
  {
    res.index=mcmc.list(mcmc(Fit3[[1]]$Monitor[,index]))
    start=2
  }

  k=length(res.index)
  cc=length(Fit3)
  for(j in start:cc)
    res.index[[j+k-start+1]]=mcmc(Fit3[[j]]$Monitor[,index])
}
plot(res.index)







