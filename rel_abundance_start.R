
model="f_sf"

is.o18bin=FALSE
if(length(grep("o18bin",model))>0)
  is.o18bin=TRUE

is.cabin=FALSE
if(length(grep("cabin",model))>0)
  is.cabin=TRUE

source("read_data.R")
source("betabin.R")
source(sprintf("model_%s.R",model))
hyper=hyper.this
hyper$N=N
hyper$N.sp=N.sp
hyper$Nf=Nf
hyper$N.forms=N.forms

library(coda)


# Collate (if multiple runs)
N.start=1
N.runs=146

for(i in 1:N.runs)
{
  show(i)
  load(sprintf("%s_%02d.Rdata",model,i+N.start-1))
  start=1
  if(i==1)
  {
    res=mcmc.list(mcmc(Fit3[[1]]$Monitor))
    start=2
  }
  cc=length(Fit3)
  k=length(res)
  for(j in start:cc)
    res[[j+k-start+1]]=mcmc(Fit3[[j]]$Monitor)
}


n1=dim(res[[1]])[2]
n2=0
for(i in 1:length(res))
{
  n2=n2+dim(res[[i]])[1]
}
res2=array(NA,c(n2,n1))
j=1
for(i in 1:length(res))
{
  nn=dim(res[[i]])[1]
  res2[j:(j+nn-1),]=res[[i]]
  j=j+nn
}

NN=dim(res2)[1]

par0=init.params.flat(hyper)
pname=names(par0)

const.occu.index=which(substr(pname,1,16)=="logit.const.occu")
const.logit.occu=res2[,const.occu.index]

const.bin.index=which(substr(pname,1,15)=="logit.const.bin")
const.logit.bin=res2[,const.bin.index]

f.occu.index=which(substr(pname,1,6)=="f.occu")
f.occu=res2[,f.occu.index]

f.bin.index=which(substr(pname,1,5)=="f.bin")
f.bin=res2[,f.bin.index]

sf.occu.index=which(substr(pname,1,7)=="sf.occu")
sf.occu=array(NA,c(NN,Nf,N.sp-1))
for(i in 1:NN)
  sf.occu[i,,]=matrix(res2[i,sf.occu.index], ncol=N.sp-1)

sf.bin.index=which(substr(pname,1,6)=="sf.bin")
sf.bin=array(NA,c(NN,Nf,N.sp))
for(i in 1:NN)
  sf.bin[i,,]=matrix(res2[i,sf.bin.index], ncol=N.sp)

o18.formation=rep(0,Nf)
ca.formation=rep(0,Nf)
for(f in 1:Nf)
{
  o18.formation[f]=y$O18.mean[y$form.nr==f][1]
  ca.formation[f]=y$ca.mean[y$form.nr==f][1]	
}


logit.occu.formation.species=array(NA,c(NN,Nf,N.sp-1))
for(i in 1:Nf)
 logit.occu.formation.species[,i,]=const.logit.occu+f.occu[,i]%*%t(rep(1,N.sp-1))+sf.occu[,i,]
occu.formation.species=ilogit(logit.occu.formation.species)

logit.psi.star=array(NA,c(NN,Nf,N.sp-1))
for(i in 1:Nf)
 logit.psi.star[,i,]=const.logit.occu+sf.occu[,i,]+f.occu[,i]%*%t(rep(1,N.sp-1))
psi.star=ilogit(logit.psi.star)



d=dim(occu.formation.species)
occu.formation.species2=array(NA,c(d[1],d[2],d[3]+1))
occu.formation.species2[,,1:d[3]]=occu.formation.species
occu.formation.species2[,,d[3]+1]=1
psi.star2=array(NA,c(d[1],d[2],d[3]+1))
psi.star2[,,1:d[3]]=psi.star
psi.star2[,,d[3]+1]=1
rm(psi.star)
rm(logit.psi.star)

cabin=array(0,c(NN,Nf,N.sp))
if(is.cabin)
{
  for(f in 1:Nf)
   for(s in 1:N.sp)
   {
     ca.o.ind=which(names(par0)==sprintf("ca.optim.bin%d",s))
     ca.w.ind=which(names(par0)==sprintf("ca.lwidth.bin%d",s))

     cabin[,f,s]=-(ca.formation[f]-res2[,ca.o.ind])^2/exp(res2[,ca.w.ind])^2
   }
}

o18bin=array(0,c(NN,Nf,N.sp))
if(is.o18bin)
{
  for(f in 1:Nf)
   for(s in 1:N.sp)
   {
     o18.o.ind=which(names(par0)==sprintf("o18.optim.bin%d",s))
     o18.w.ind=which(names(par0)==sprintf("o18.lwidth.bin%d",s))

     o18bin[,f,s]=-(o18.formation[f]-res2[,o18.o.ind])^2/exp(res2[,o18.w.ind])^2
   }
}


logit.bin.formation.species=array(NA,c(NN,Nf,N.sp))
for(i in 1:Nf)
 logit.bin.formation.species[,i,]=const.logit.bin+f.bin[,i]%*%t(rep(1,N.sp))+sf.bin[,i,]+o18bin[,i,]+cabin[,i,]
bin.formation.species=ilogit(logit.bin.formation.species)

logit.p.star=array(NA,c(NN,Nf,N.sp))
for(i in 1:Nf)
 logit.p.star[,i,]=const.logit.bin+sf.bin[,i,]+o18bin[,i,]+cabin[,i,]
p.star=ilogit(logit.p.star)

rm(cabin)
rm(o18bin)


lower=function(x) quantile(x, 0.025)
upper=function(x) quantile(x, 0.975)

forms=sort(unique(y$form.nr))
Nf2=length(forms) # check that this is the same as Nf
Nf3=max(forms) # check that this is also the same as Nf

K1=rep(NA,Nf)
K2=rep(NA,Nf)
for(i in 1:Nf)
{
  K1[i]=y$time.start[y$form.nr==i][1]
  K2[i]=y$time.end[y$form.nr==i][1]
}



K.mid=(K1+K2)/2

K=K.mid
Ks=K1
Ke=K2

 
sp=names(y)[22:25]



#####################
# Summary plots:
#####################


# Raw detection rate data:

mean.raw=array(NA,c(Nf,N.sp))
for(s in 1:N.sp)
 for(f in 1:Nf)
   mean.raw[f,s]=mean(y[y$form.nr==f,names(y)==sp[s]]/y$Total[y$form.nr==f])
  

# Formation-wise adjustments to overall occupancy:
f.occu.mean=apply(f.occu,2,mean)
f.occu.expmean=apply(exp(f.occu),2,mean)
f.occu.median=apply(f.occu,2,median)
f.occu.lower=apply(f.occu,2,lower)
f.occu.upper=apply(f.occu,2,upper)

# Formation-wise adjustments to overall detection:
f.bin.mean=apply(f.bin,2,mean)
f.bin.expmean=apply(exp(f.bin),2,mean)
f.bin.median=apply(f.bin,2,median)
f.bin.lower=apply(f.bin,2,lower)
f.bin.upper=apply(f.bin,2,upper)



sf.occu.mean=array(NA,c(N.sp-1,Nf))
sf.occu.expmean=array(NA,c(N.sp-1,Nf))
sf.occu.median=array(NA,c(N.sp-1,Nf))
sf.occu.lower=array(NA,c(N.sp-1,Nf))
sf.occu.upper=array(NA,c(N.sp-1,Nf))

for(s in 1:N.sp-1)
{
  sf.occu.mean[s,]=apply(sf.occu[,,s],2,mean)
  sf.occu.expmean[s,]=apply(exp(sf.occu[,,s]),2,mean)
  sf.occu.median[s,]=apply(sf.occu[,,s],2,median)
  sf.occu.lower[s,]=apply(sf.occu[,,s],2,lower)
  sf.occu.upper[s,]=apply(sf.occu[,,s],2,upper)
}  

sf.bin.mean=array(NA,c(N.sp,Nf))
sf.bin.expmean=array(NA,c(N.sp,Nf))
sf.bin.median=array(NA,c(N.sp,Nf))
sf.bin.lower=array(NA,c(N.sp,Nf))
sf.bin.upper=array(NA,c(N.sp,Nf))

for(s in 1:N.sp)
{
  sf.bin.mean[s,]=apply(sf.bin[,,s],2,mean)
  sf.bin.expmean[s,]=apply(exp(sf.bin[,,s]),2,mean)
  sf.bin.median[s,]=apply(sf.bin[,,s],2,median)
  sf.bin.lower[s,]=apply(sf.bin[,,s],2,lower)
  sf.bin.upper[s,]=apply(sf.bin[,,s],2,upper)
}  




# Per formation per species occupancy

# Make raw data illustration data:
occu.raw=array(NA,c(Nf,N.sp-1))
for(s in 1:(N.sp-1))
 for(f in 1:Nf)
 {
   ind=which(names(y)==sp[s])
   occu.raw[f,s]=mean(y[y$form.nr==f,ind]>0)
 }

occu.sf.mean=array(NA,c(N.sp-1,Nf))
occu.sf.expmean=array(NA,c(N.sp-1,Nf))
occu.sf.median=array(NA,c(N.sp-1,Nf))
occu.sf.lower=array(NA,c(N.sp-1,Nf))
occu.sf.upper=array(NA,c(N.sp-1,Nf))

for(s in 1:N.sp-1)
{
  occu.sf.mean[s,]=apply(occu.formation.species[,,s],2,mean)
  occu.sf.expmean[s,]=apply(exp(occu.formation.species[,,s]),2,mean)
  occu.sf.median[s,]=apply(occu.formation.species[,,s],2,median)
  occu.sf.lower[s,]=apply(occu.formation.species[,,s],2,lower)
  occu.sf.upper[s,]=apply(occu.formation.species[,,s],2,upper)
}  

rm(occu.formation.species)

bin.sf.mean=array(NA,c(N.sp,Nf))
bin.sf.expmean=array(NA,c(N.sp,Nf))
bin.sf.median=array(NA,c(N.sp,Nf))
bin.sf.lower=array(NA,c(N.sp,Nf))
bin.sf.upper=array(NA,c(N.sp,Nf))

for(s in 1:N.sp)
{
  bin.sf.mean[s,]=apply(bin.formation.species[,,s],2,mean)
  bin.sf.expmean[s,]=apply(exp(bin.formation.species[,,s]),2,mean)
  bin.sf.median[s,]=apply(bin.formation.species[,,s],2,median)
  bin.sf.lower[s,]=apply(bin.formation.species[,,s],2,lower)
  bin.sf.upper[s,]=apply(bin.formation.species[,,s],2,upper)
}  



pstar.mean=array(NA,c(N.sp,Nf))
pstar.median=array(NA,c(N.sp,Nf))
pstar.lower=array(NA,c(N.sp,Nf))
pstar.upper=array(NA,c(N.sp,Nf))

for(s in 1:N.sp)
{
  pstar.mean[s,]=apply(p.star[,,s],2,mean)  
  pstar.median[s,]=apply(p.star[,,s],2,median)
  pstar.lower[s,]=apply(p.star[,,s],2,lower)
  pstar.upper[s,]=apply(p.star[,,s],2,upper)
}  





# Now make relative abundance, formation by formation:
# abundance.given.occupancy=bin.formation.species/(1-bin.formation.species)
#abundance.given.occupancy=-log(1-bin.formation.species)
#abundance.given.occupancy=log(1+exp(logit.bin.formation.species))
#abundance=abundance.given.occupancy*occu.formation.species2

abundance=log(1+exp(logit.p.star))*psi.star2



rel.abundance=abundance
for(i in 1:Nf)
  rel.abundance[,i,]=abundance[,i,]/apply(abundance[,i,],1,sum)

mean.rel.abundance=array(NA,c(Nf,N.sp))
median.rel.abundance=array(NA,c(Nf,N.sp))
lower.rel.abundance=array(NA,c(Nf,N.sp))
upper.rel.abundance=array(NA,c(Nf,N.sp))
for(f in 1:Nf)
  for(s in 1:N.sp)
  {
    mean.rel.abundance[f,s]=mean(rel.abundance[,f,s])
    median.rel.abundance[f,s]=median(rel.abundance[,f,s])
    lower.rel.abundance[f,s]=lower(rel.abundance[,f,s])
    upper.rel.abundance[f,s]=upper(rel.abundance[,f,s])
  }
  




A=array(NA,c(N.sp,Nf))
for(i in 1:Nf)
  A[,i]=apply(rel.abundance[,i,],2,mean)
species=c("Ant.Ton.","Esc.Ex.","Ara.Uni.","Super")


rm(rel.abundance)


Q=abundance
for(s in 1:S)
  Q[,,s]=abundance[,,s]/apply(abundance[,,s],1,mean)
mean.Q=array(NA,c(Nf,N.sp))
median.Q=array(NA,c(Nf,N.sp))
lower.Q=array(NA,c(Nf,N.sp))
upper.Q=array(NA,c(Nf,N.sp))
for(f in 1:Nf)
  for(s in 1:N.sp)
  {
    mean.Q[f,s]=mean(Q[!is.infinite(Q[,f,s]),f,s],na.rm=T)
    median.Q[f,s]=median(Q[,f,s],na.rm=T)
    lower.Q[f,s]=lower(Q[,f,s])
    upper.Q[f,s]=upper(Q[,f,s])
  }

rm(Q)






rm(abundance)


# Raw detection rate data with ocpcuancy*detection

detection.given.occupancy=bin.formation.species
detection=detection.given.occupancy*occu.formation.species2

mean.raw=array(NA,c(Nf,N.sp))
for(s in 1:N.sp)
 for(f in 1:Nf)
   mean.raw[f,s]=mean(y[y$form.nr==f,names(y)==sp[s]]/y$Total[y$form.nr==f])



detection.sf.mean=array(NA,c(N.sp,Nf))
detection.sf.expmean=array(NA,c(N.sp,Nf))
detection.sf.median=array(NA,c(N.sp,Nf))
detection.sf.lower=array(NA,c(N.sp,Nf))
detection.sf.upper=array(NA,c(N.sp,Nf))

for(s in 1:N.sp)
{
  detection.sf.mean[s,]=apply(detection[,,s],2,mean)
  detection.sf.expmean[s,]=apply(exp(detection[,,s]),2,mean)
  detection.sf.median[s,]=apply(detection[,,s],2,median)
  detection.sf.lower[s,]=apply(detection[,,s],2,lower)
  detection.sf.upper[s,]=apply(detection[,,s],2,upper)
}  


# Remove all MCMC stuff (the big chunks):

rm(bin.formation.species)
rm(d)
rm(const.logit.bin)
rm(const.logit.occu)
rm(detection)
rm(detection.given.occupancy)
rm(f.bin)
rm(f.occu)
rm(Fit3)
rm(logit.bin.formation.species)
rm(logit.occu.formation.species)
rm(occu.formation.species2)
rm(res)
rm(res2)
rm(sf.bin)
rm(sf.occu)
rm(psi.star2)
rm(p.star)
rm(logit.p.star)


save.image(file="plotting_variables.Rdata")
