
model="f_sf_cabin"

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
N.runs=40

for(i in 1:N.runs)
{
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
logit.psi.star=array(NA,c(NN,Nf,N.sp-1))
for(i in 1:Nf)
{
  logit.occu.formation.species[,i,]=
    const.logit.occu+
    f.occu[,i]%*%t(rep(1,N.sp-1))+sf.occu[,i,]
  logit.psi.star[,i,]=
    const.logit.occu+
    f.occu[,i]%*%t(rep(1,N.sp-1))+sf.occu[,i,]

}
occu.formation.species=ilogit(logit.occu.formation.species)
psi.star=ilogit(logit.occu.formation.species)
rm(logit.occu.formation.species)


d=dim(occu.formation.species)
occu.formation.species2=array(NA,c(d[1],d[2],d[3]+1))
occu.formation.species2[,,1:d[3]]=occu.formation.species
occu.formation.species2[,,d[3]+1]=1
rm(occu.formation.species)

psi.star2=array(NA,c(d[1],d[2],d[3]+1))
psi.star2[,,1:d[3]]=psi.star
psi.star2[,,d[3]+1]=1
rm(psi.star)




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

logit.bin.formation.species=array(NA,c(NN,Nf,N.sp))
logit.p.star=array(NA,c(NN,Nf,N.sp))
for(i in 1:Nf)
{
  logit.bin.formation.species[,i,]=
    const.logit.bin+
    f.bin[,i]%*%t(rep(1,N.sp))+
    sf.bin[,i,]+o18bin[,i,]+cabin[,i,]
  logit.p.star[,i,]=
    const.logit.bin+
    sf.bin[,i,]+o18bin[,i,]+cabin[,i,]
}
bin.formation.species=ilogit(logit.bin.formation.species)
p.star=ilogit(logit.p.star)

rm(logit.bin.formation.species)
rm(logit.p.star)



lower=function(x) quantile(x, 0.025, na.rm=T)
upper=function(x) quantile(x, 0.975, na.rm=T)

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


# Plot raw detection rate data:

mean.raw=array(NA,c(Nf,N.sp))
for(s in 1:N.sp)
 for(f in 1:Nf)
   mean.raw[f,s]=mean(y[y$form.nr==f,names(y)==sp[s]]/y$Total[y$form.nr==f])
  
png("raw_data.png",height=1200,width=1600)
par(cex=3) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
par(mfrow=c(2,2)) 
plot(K[y$form.nr], y[,names(y)==sp[1]]/y$Total,xlab="Time",ylab=sprintf("Raw ratio: %s",sp[1]),ylim=c(0,1))
lines(K, mean.raw[,1])
plot(K[y$form.nr], y[,names(y)==sp[2]]/y$Total,xlab="Time",ylab=sprintf("Raw ratio: %s",sp[2]),ylim=c(0,1))
lines(K, mean.raw[,2])
plot(K[y$form.nr], y[,names(y)==sp[3]]/y$Total,xlab="Time",ylab=sprintf("Raw ratio: %s",sp[3]),ylim=c(0,1))
lines(K, mean.raw[,3])
plot(K[y$form.nr], y[,names(y)==sp[4]]/y$Total,xlab="Time",ylab=sprintf("Raw ratio: %s",sp[4]),ylim=c(0,1))
lines(K, mean.raw[,4])
dev.off()




# Formation-wise adjustments to overall occupancy:
png("log_odds_occupancy_formation.png",height=1200,width=1600)
par(cex=3) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
plot(K, apply(f.occu,2,mean), type="b", ylim=c(-5,7),
  xlab="time", ylab="Log odds ratio, occupancy",lwd=5)
lines(K, apply(f.occu,2,median),type="b",lty=2,lwd=2)
lines(K, apply(f.occu,2,lower),type="b",lty=2,lwd=2)
lines(K, apply(f.occu,2,upper),type="b",lty=2,lwd=2)
dev.off()

png("odds_occupancy_formation.png",height=1200,width=1600)
par(cex=3) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
plot(K, exp(apply(f.occu,2,mean)), type="b", ylim=c(0.01,500),
  xlab="time", ylab="Odds ratio, occupancy",log="y",lwd=5)
lines(K, apply(exp(f.occu),2,median),type="b",lty=2,lwd=2)
lines(K, apply(exp(f.occu),2,lower),type="b",lty=2,lwd=2)
lines(K, apply(exp(f.occu),2,upper),type="b",lty=2,lwd=2)
dev.off()

# Formation-wise adjustments to detection:
png("log_odds_detection_formation.png",height=1200,width=1600)
par(cex=3) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
plot(K, apply(f.bin,2,mean), type="b", ylim=c(-2.4,2),
  xlab="time", ylab="Log-ddds ratio, detection",lwd=5)
lines(K, apply(f.bin,2,median),type="b",lty=2,lwd=2)
lines(K, apply(f.bin,2,lower),type="b",lty=2,lwd=2)
lines(K, apply(f.bin,2,upper),type="b",lty=2,lwd=2)
dev.off()

png("odds_detection_formation.png",height=1200,width=1600)
par(cex=3) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
plot(K, exp(apply(f.bin,2,mean)), type="b", ylim=c(0.1,10),
  xlab="time", ylab="Odds ratio, detection",log="y",lwd=5)
lines(K, exp(apply(f.bin,2,median)),type="b",lty=2,lwd=2)
lines(K, exp(apply(f.bin,2,lower)),type="b",lty=2,lwd=2)
lines(K, exp(apply(f.bin,2,upper)),type="b",lty=2,lwd=2)
dev.off()



# Formation-species-wise adjustments to overall occupancy:
png("log_odds_occupancy_species_formation.png",height=1200,width=1600)
par(cex=3, mfrow=c(2,2)) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
plot(K, apply(sf.occu[,,1],2,mean), type="b", ylim=c(-4,4),
  xlab="time", ylab="Log odds ratio, occupancy",lwd=5,main=sp[1])
lines(K, apply(sf.occu[,,1],2,median),type="b",lty=2,lwd=2)
lines(K, apply(sf.occu[,,1],2,lower),type="b",lty=2,lwd=2)
lines(K, apply(sf.occu[,,1],2,upper),type="b",lty=2,lwd=2)
plot(K, apply(sf.occu[,,2],2,mean), type="b", ylim=c(-31,31),
  xlab="time", ylab="Log odds ratio, occupancy",lwd=5,main=sp[2])
lines(K, apply(sf.occu[,,2],2,median),type="b",lty=2,lwd=2)
lines(K, apply(sf.occu[,,2],2,lower),type="b",lty=2,lwd=2)
lines(K, apply(sf.occu[,,2],2,upper),type="b",lty=2,lwd=2)
plot(K, apply(sf.occu[,,1],2,mean), type="b", ylim=c(-25,25),
  xlab="time", ylab="Log odds ratio, occupancy",lwd=5,main=sp[3])
lines(K, apply(sf.occu[,,3],2,median),type="b",lty=2,lwd=2)
lines(K, apply(sf.occu[,,3],2,lower),type="b",lty=2,lwd=2)
lines(K, apply(sf.occu[,,3],2,upper),type="b",lty=2,lwd=2)
dev.off()

png("odds_occupancy_species_formation.png",height=1200,width=1600)
par(cex=3, mfrow=c(2,2)) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
plot(K, exp(apply(sf.occu[,,1],2,mean)), type="b", ylim=c(0.01,100),
  xlab="time", ylab="Odds ratio, occupancy",log="y",lwd=5,main=sp[1])
lines(K, exp(apply(sf.occu[,,1],2,median)),type="b",lty=2,lwd=2)
lines(K, exp(apply(sf.occu[,,1],2,lower)),type="b",lty=2,lwd=2)
lines(K, exp(apply(sf.occu[,,1],2,upper)),type="b",lty=2,lwd=2)
plot(K, exp(apply(sf.occu[,,2],2,mean)), type="b", ylim=c(0.000000000001,5000000000000),
  xlab="time", ylab="Odds ratio, occupancy",log="y",lwd=5,main=sp[2])
lines(K, exp(apply(sf.occu[,,2],2,median)),type="b",lty=2,lwd=2)
lines(K, exp(apply(sf.occu[,,2],2,lower)),type="b",lty=2,lwd=2)
lines(K, exp(apply(sf.occu[,,2],2,upper)),type="b",lty=2,lwd=2)
plot(K, exp(apply(sf.occu[,,3],2,mean)), type="b", ylim=c(0.000000001,5000000000),
  xlab="time", ylab="Odds ratio, occupancy",log="y",lwd=5,main=sp[3])
lines(K, exp(apply(sf.occu[,,3],2,median)),type="b",lty=2,lwd=2)
lines(K, exp(apply(sf.occu[,,3],2,lower)),type="b",lty=2,lwd=2)
lines(K, exp(apply(sf.occu[,,3],2,upper)),type="b",lty=2,lwd=2)
dev.off()

# Formation-species-wise adjustments to detection:
png("log_odds_detection_species_formation.png",height=1200,width=1600)
par(cex=3, mfrow=c(2,2)) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
plot(K, apply(sf.bin[,,1],2,mean), type="b", ylim=c(-2.4,2),
  xlab="time", ylab="Log-ddds ratio, detection",lwd=5,main=sp[1])
lines(K, apply(sf.bin[,,1],2,median),type="b",lty=2,lwd=2)
lines(K, apply(sf.bin[,,1],2,lower),type="b",lty=2,lwd=2)
lines(K, apply(sf.bin[,,1],2,upper),type="b",lty=2,lwd=2)
plot(K, apply(sf.bin[,,2],2,mean), type="b", ylim=c(-8,5),
  xlab="time", ylab="Log-ddds ratio, detection",lwd=5,main=sp[2])
lines(K, apply(sf.bin[,,2],2,median),type="b",lty=2,lwd=2)
lines(K, apply(sf.bin[,,2],2,lower),type="b",lty=2,lwd=2)
lines(K, apply(sf.bin[,,2],2,upper),type="b",lty=2,lwd=2)
plot(K, apply(sf.bin[,,3],2,mean), type="b", ylim=c(-4,4),
  xlab="time", ylab="Log-ddds ratio, detection",lwd=5,main=sp[3])
lines(K, apply(sf.bin[,,3],2,median),type="b",lty=2,lwd=2)
lines(K, apply(sf.bin[,,3],2,lower),type="b",lty=2,lwd=2)
lines(K, apply(sf.bin[,,3],2,upper),type="b",lty=2,lwd=2)
plot(K, apply(sf.bin[,,4],2,mean), type="b", ylim=c(-1.5,1.5),
  xlab="time", ylab="Log-ddds ratio, detection",lwd=5,main=sp[4])
lines(K, apply(sf.bin[,,4],2,median),type="b",lty=2,lwd=2)
lines(K, apply(sf.bin[,,4],2,lower),type="b",lty=2,lwd=2)
lines(K, apply(sf.bin[,,4],2,upper),type="b",lty=2,lwd=2)
dev.off()

png("odds_detection_species_formation.png",height=1200,width=1600)
par(cex=3, mfrow=c(2,2)) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
plot(K, exp(apply(sf.bin[,,1],2,mean)), type="b", ylim=c(0.1,10),
  xlab="time", ylab="Odds ratio, detection",log="y",lwd=5,main=sp[1])
lines(K, exp(apply(sf.bin[,,1],2,median)),type="b",lty=2,lwd=2)
lines(K, exp(apply(sf.bin[,,1],2,lower)),type="b",lty=2,lwd=2)
lines(K, exp(apply(sf.bin[,,1],2,upper)),type="b",lty=2,lwd=2)
plot(K, exp(apply(sf.bin[,,2],2,mean)), type="b", ylim=c(0.001,100),
  xlab="time", ylab="Odds ratio, detection",log="y",lwd=5,main=sp[2])
lines(K, exp(apply(sf.bin[,,2],2,median)),type="b",lty=2,lwd=2)
lines(K, exp(apply(sf.bin[,,2],2,lower)),type="b",lty=2,lwd=2)
lines(K, exp(apply(sf.bin[,,2],2,upper)),type="b",lty=2,lwd=2)
plot(K, exp(apply(sf.bin[,,3],2,mean)), type="b", ylim=c(0.01,100),
  xlab="time", ylab="Odds ratio, detection",log="y",lwd=5,main=sp[3])
lines(K, exp(apply(sf.bin[,,3],2,median)),type="b",lty=2,lwd=2)
lines(K, exp(apply(sf.bin[,,3],2,lower)),type="b",lty=2,lwd=2)
lines(K, exp(apply(sf.bin[,,3],2,upper)),type="b",lty=2,lwd=2)
plot(K, exp(apply(sf.bin[,,4],2,mean)), type="b", ylim=c(0.1,10),
  xlab="time", ylab="Odds ratio, detection",log="y",lwd=5,main=sp[4])
lines(K, exp(apply(sf.bin[,,4],2,median)),type="b",lty=2,lwd=2)
lines(K, exp(apply(sf.bin[,,4],2,lower)),type="b",lty=2,lwd=2)
lines(K, exp(apply(sf.bin[,,4],2,upper)),type="b",lty=2,lwd=2)
dev.off()




# Per formation per species occupancy

# Make raw data illustration data:
occu.raw=array(NA,c(Nf,N.sp-1))
for(s in 1:(N.sp-1))
 for(f in 1:Nf)
 {
   ind=which(names(y)==sp[s])
   occu.raw[f,s]=mean(y[y$form.nr==f,ind]>0)
 }

# All species:
png("occupancy_species_formation.png",height=1200,width=1600)
par(cex=3,cex.lab=2.9,cex.sub=1.9,cex.axis=1.9,cex.main=1.9,mar=c(10,5,6,1))
cols=c("red","green","blue", "black")
tstr=sprintf("Red=%s,Green=%s,Blue=%s",sp[1],sp[2],sp[3])
par(mfrow=c(1,1))
plot(K, apply(occu.formation.species[,,1],2,mean), type="b", ylim=c(0.4,1),
  xlab="time", ylab="Occupancy",main=tstr,col=cols[1],lwd=4,log="y")
lines(K, apply(occu.formation.species[,,2],2,mean),col=cols[2],lwd=4)
lines(K, apply(occu.formation.species[,,3],2,mean), type="b",col=cols[3],lwd=4)
dev.off()

png("occupancy_species_formation_uncert.png",height=1200,width=1600)
par(cex=3,cex.lab=2.9,cex.sub=1.9,cex.axis=1.9,mar=c(10,5,1,1))
par(mfrow=c(2,2))
plot(K, apply(occu.formation.species[,,1],2,mean), type="b", ylim=c(0,1),
  xlab="time", ylab=sprintf("Occupancy: %s",sp[1]),lwd=4)
lines(K, occu.raw[,1],lty=2,lwd=2,type="b",pch=2)
lines(K, apply(occu.formation.species[,,1],2,lower),lty=3,lwd=3,type="b",pch=4)
lines(K, apply(occu.formation.species[,,1],2,upper),lty=3,lwd=3,type="b",pch=4)

plot(K, apply(occu.formation.species[,,2],2,mean), type="b", ylim=c(0,1),
  xlab="time", ylab=sprintf("Occupancy: %s",sp[2]),lwd=4)
lines(K, occu.raw[,2],lty=2,lwd=2,type="b",pch=2)
lines(K, apply(occu.formation.species[,,2],2,lower),lty=3,lwd=3,type="b",pch=4)
lines(K, apply(occu.formation.species[,,2],2,upper),lty=3,lwd=3,type="b",pch=4)

plot(K, apply(occu.formation.species[,,3],2,mean), type="b", ylim=c(0,1),
  xlab="time", ylab=sprintf("Occupancy: %s",sp[3]),lwd=4)
lines(K, occu.raw[,3],lty=2,lwd=2,type="b",pch=2)
lines(K, apply(occu.formation.species[,,3],2,lower),lty=3,lwd=3,type="b",pch=4)
lines(K, apply(occu.formation.species[,,3],2,upper),lty=3,lwd=3,type="b",pch=4)
dev.off()




# Per formation per species detection
png("detection_species_formation.png",height=1200,width=1600)
par(cex=3,cex.lab=2.9,cex.sub=1.9,cex.axis=1.9,cex.main=1.9,mar=c(10,5,6,1))
cols=c("red","green","blue", "black")
tstr=sprintf("red=%s,green=%s\nblue=%s,black=superspecies",sp[1],sp[2],sp[3])
par(mfrow=c(1,1))
plot(K, apply(bin.formation.species[,,1],2,mean), type="b", ylim=c(0.005,1),
  xlab="time", ylab="Detection",main=tstr,col=cols[1],lwd=4,log="y")
lines(K, apply(bin.formation.species[,,2],2,mean),col=cols[2],lwd=4)
lines(K, apply(bin.formation.species[,,3],2,mean), type="b",col=cols[3],lwd=4)
lines(K, apply(bin.formation.species[,,4],2,mean), type="b",col=cols[4],lwd=4)
dev.off()

png("detection_species_formation_uncert.png",height=1200,width=1600)
par(cex=3,cex.lab=2.9,cex.sub=1.9,cex.axis=1.9,mar=c(10,5,1,1))
par(mfrow=c(2,2))
plot(K, apply(bin.formation.species[,,1],2,mean), type="b", ylim=c(0.0001,1),
  xlab="time", ylab=sprintf("Detection: %s",sp[1]),lwd=4,log="y")
lines(K, apply(bin.formation.species[,,1],2,median),lty=2)
lines(K, apply(bin.formation.species[,,1],2,lower),lty=2)
lines(K, apply(bin.formation.species[,,1],2,upper),lty=2)

plot(K, apply(bin.formation.species[,,2],2,mean), type="b", ylim=c(0.0001,1),
  xlab="time", ylab=sprintf("Detection: %s",sp[2]),lwd=4,log="y")
lines(K, apply(bin.formation.species[,,2],2,median),lty=2)
lines(K, apply(bin.formation.species[,,2],2,lower),lty=2)
lines(K, apply(bin.formation.species[,,2],2,upper),lty=2)

plot(K, apply(bin.formation.species[,,3],2,mean), type="b", ylim=c(0.0001,1),
  xlab="time", ylab=sprintf("Detection: %s",sp[3]),lwd=4,log="y")
lines(K, apply(bin.formation.species[,,3],2,median),lty=2)
lines(K, apply(bin.formation.species[,,3],2,lower),lty=2)
lines(K, apply(bin.formation.species[,,3],2,upper),lty=2)

plot(K, apply(bin.formation.species[,,4],2,mean), type="b", ylim=c(0.0001,1),
  xlab="time", ylab=sprintf("Detection: %s","Superspecies"),lwd=4,log="y")
lines(K, apply(bin.formation.species[,,4],2,median),lty=2)
lines(K, apply(bin.formation.species[,,4],2,lower),lty=2)
lines(K, apply(bin.formation.species[,,4],2,upper),lty=2)
dev.off()

rm(bin.formation.species)
rm(occu.formation.species2)






# Now make relative abundance, formation by formation:
# abundance.given.occupancy=bin.formation.species/(1-bin.formation.species)
#abundance.given.occupancy=-log(1-bin.formation.species)
#abundance.given.occupancy=log(1+exp(logit.bin.formation.species))
#abundance=abundance.given.occupancy*occu.formation.species2


abundance=-log(1-p.star)*psi.star2




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
save(mean.rel.abundance, file=sprintf("mean_ra_%s.Rdata",model))
save(median.rel.abundance, file=sprintf("median_ra_%s.Rdata",model))
save(lower.rel.abundance, file=sprintf("lower_ra_%s.Rdata",model))
save(upper.rel.abundance, file=sprintf("upper_ra_%s.Rdata",model))




png("rel_abundance.png",height=1200,width=1600)
par(cex=3,cex.lab=2.9,cex.sub=1.9,cex.axis=2.2,cex.main=1.9,mar=c(10,5,6,1))
cols=c("red","green","blue", "black")
tstr=sprintf("red=%s,green=%s\nblue=%s,black=superspecies",sp[1],sp[2],sp[3])
par(mfrow=c(1,1))
plot(K, mean.rel.abundance[,1], type="b", ylim=c(0.001,1),
  xlab="time", ylab="Rel. Abundance",main=tstr,col=cols[1],lwd=4,log="y")
lines(K, mean.rel.abundance[,2], type="b",col=cols[2],lwd=4)
lines(K, mean.rel.abundance[,3], type="b",col=cols[3],lwd=4)
lines(K, mean.rel.abundance[,4], type="b",col=cols[4],lwd=4)
dev.off()





png("rel_abundance_uncert.png",height=1200,width=1600)
par(cex=3,cex.lab=2.9,cex.sub=1.9,cex.axis=1.9,mar=c(10,5,1,1))
par(mfrow=c(2,2))
plot(K, mean.rel.abundance[,1], type="b", ylim=c(0.0001,1),
  xlab="time", ylab=sprintf("Relative abundance: %s",sp[1]),lwd=4,log="y")
lines(K, median.rel.abundance[,1],lty=2)
lines(K, lower.rel.abundance[,1],lty=2)
lines(K, upper.rel.abundance[,1],lty=2)

plot(K, mean.rel.abundance[,2], type="b", ylim=c(0.0001,1),
  xlab="time", ylab=sprintf("Relative abundance: %s",sp[2]),lwd=4,log="y")
lines(K, median.rel.abundance[,2],lty=2)
lines(K, lower.rel.abundance[,2],lty=2)
lines(K, upper.rel.abundance[,2],lty=2)

plot(K, mean.rel.abundance[,3], type="b", ylim=c(0.0001,1),
  xlab="time", ylab=sprintf("Relative abundance: %s",sp[3]),lwd=4,log="y")
lines(K, median.rel.abundance[,3],lty=2)
lines(K, lower.rel.abundance[,3],lty=2)
lines(K, upper.rel.abundance[,3],lty=2)

plot(K, mean.rel.abundance[,4], type="b", ylim=c(0.0001,1),
  xlab="time", ylab=sprintf("Relative abundance: %s","Superspecies"),lwd=4,log="y")
lines(K, median.rel.abundance[,4],lty=2)
lines(K, lower.rel.abundance[,4],lty=2)
lines(K, upper.rel.abundance[,4],lty=2)
dev.off()



A=array(NA,c(N.sp,Nf))
for(i in 1:Nf)
  A[,i]=apply(rel.abundance[,i,],2,mean)
species=c("Ant.Ton.","Esc.Ex.","Ara.Uni.","Super")

png("rel_abundance_formation.png",height=1200,width=1600)
par(cex=3,cex.lab=2, cex.sub=2,cex.main=2,cex.axis=2)
par(mfrow=c(3,3))
for(i in 1:Nf)
   barplot(A[1:(N.sp-1),i],names.arg=species[1:(N.sp-1)],
     main=sprintf("Formation age %2.2f-%2.2f",-Ke[i], -Ks[i]))
dev.off()






# Plot raw detection rate data with ocpcuancy*detection

detection.given.occupancy=bin.formation.species
detection=detection.given.occupancy*occu.formation.species2

mean.raw=array(NA,c(Nf,N.sp))
for(s in 1:N.sp)
 for(f in 1:Nf)
   mean.raw[f,s]=mean(y[y$form.nr==f,names(y)==sp[s]]/y$Total[y$form.nr==f])
  
png("raw_data_with_occu_times_det.png",height=1200,width=1600)
par(cex=4,cex.lab=1.5,cex.sub=1.5,cex.axis=1.5,mfrow=c(2,2),
  mar=c(5,8,4,2)) 
plot(K[y$form.nr], y[,names(y)==sp[1]]/y$Total,xlab="Time",ylab=sprintf("Ratio: %s",sp[1]),ylim=c(0,0.5),lwd=2)
lines(K, mean.raw[,1],lwd=2)
lines(K, apply(detection[,,1],2,mean), lty=3, lwd=3)
lines(K, apply(detection[,,1],2,lower), lty=2, lwd=1)
lines(K, apply(detection[,,1],2,upper), lty=2, lwd=1)
plot(K[y$form.nr], y[,names(y)==sp[2]]/y$Total,xlab="Time",ylab=sprintf("Ratio: %s",sp[2]),ylim=c(0,0.5),lwd=2)
lines(K, mean.raw[,2],lwd=2)
lines(K, apply(detection[,,2],2,mean), lty=3, lwd=3)
lines(K, apply(detection[,,2],2,lower), lty=2, lwd=1)
lines(K, apply(detection[,,2],2,upper), lty=2, lwd=1)
plot(K[y$form.nr], y[,names(y)==sp[3]]/y$Total,xlab="Time",ylab=sprintf("Ratio: %s",sp[3]),ylim=c(0,0.5),lwd=2)
lines(K, mean.raw[,3],lwd=2)
lines(K, apply(detection[,,3],2,mean), lty=3, lwd=3)
lines(K, apply(detection[,,3],2,lower), lty=2, lwd=1)
lines(K, apply(detection[,,3],2,upper), lty=2, lwd=1)
plot(K[y$form.nr], y[,names(y)==sp[4]]/y$Total,xlab="Time",ylab=sprintf("Ratio: %s",sp[4]),ylim=c(0,1),lwd=2)
lines(K, mean.raw[,4],lwd=2)
lines(K, apply(detection[,,4],2,mean), lty=3, lwd=3)
lines(K, apply(detection[,,4],2,lower), lty=2, lwd=1)
lines(K, apply(detection[,,4],2,upper), lty=2, lwd=1)
dev.off()

