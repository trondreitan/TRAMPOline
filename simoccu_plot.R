# Read mean parameters for full model:
load("../par_mean_f_sf.RData")

# Some extra info:
load("plotting_variables.Rdata")

load("par_sim.RData")

n.shell=60
n.site=10
F=9

N.start=1
N.sim=100

library("coda")
library("LaplacesDemon")


par.mean=list(p1=par.sim)
par.median=list(p1=par.sim)

for(i in 1:N.sim)
{
  show(i)
  load(sprintf("f_sf_r%02d.Rdata",N.start+i-1))

  n=length(Fit3)
  m=dim(Fit3[[1]]$Monitor)[1]
  k=dim(Fit3[[1]]$Monitor)[2]

  res=array(0,c(n,m,k))
  res2=list(res2=mcmc(Fit3[[1]]$Monitor))
  res3=array(0,c(m*n,k))
  for(j in 1:n)
  {
    res[j,,]=mcmc(Fit3[[j]]$Monitor)
    if(j>1)
      res2[[j]]=mcmc(Fit3[[j]]$Monitor)
    res3[(j-1)*m+1:m,]=mcmc(Fit3[[j]]$Monitor)
  }

  show(gelman.diag(res2)$mpsrf)

  par.mean[[i]]=par.sim
  par.median[[i]]=par.sim
  for(j in 1:k)
  {
    par.mean[[i]][j]=mean(res3[,j])
    par.median[[i]][j]=median(res3[,j])
  }
}


p.par=0*par.sim
bias.par=0*par.sim
sd.par=0*par.sim
for(j in 1:length(par.sim))
{
  x=rep(0,N.sim)
  for(i in 1:N.sim)
    x[i]=par.mean[[i]][j]
  p.par[j]=as.numeric(t.test(x,mu=par.sim[j])$p.value)
  bias.par[j]=par.sim[j]-mean(x)
  sd.par[j]=sd(x)
}
p.non.sd=p.par[substr(names(p.par),1,6)!="log.sd"]
bias.non.sd=bias.par[substr(names(p.par),1,6)!="log.sd"]
sd.non.sd=sd.par[substr(names(p.par),1,6)!="log.sd"]

sum(p.non.sd<0.05/length(p.non.sd))
# 25



p.non.sd[which(p.non.sd<0.05/length(p.non.sd))]


N.start1=1
N.end1=100


png("par_sim_01.png",height=3000,width=3000)
par(mfrow=c(5,1),cex=2)
par(cex=2.9, cex.main=2.9, 
    cex.lab=2.9, cex.sub=2.0, cex.axis=1.6,
    #cex.lab=2.9,cex.sub=2.4,cex.axis=2.4,
    mar=c(4,7,2,1), 
    oma=c(2,3,2,1),font.main = 3)
for(g in 1:5)
{
  jj=(g-1)*4+1:4
  show(jj)
  plot(as.factor(sprintf("%02d-%s",jj,
    names(par.sim)[jj])),
    as.numeric(par.sim)[jj],type="h",ylim=c(-7,3),
    ylab="",xlab="", axes=FALSE, lwd=3)
  axis(1, at=1:4, labels=names(par.sim[jj]))
  axis(2)
  box()
  for(j in jj)
   for(i in N.start1:N.end1)
    points(j-min(jj)+1+rnorm(1,0,0.05),par.mean[[i]][j])
}
mtext("Parameters", outer=T, padj=0, side=1, cex=9.5)
mtext("Parameter value", outer=T, padj=0.5, side=2, cex=9.5)
dev.off()

png("par_sim_02.png",height=2000,width=4000)
par(mfrow=c(2,1))
par(cex=2.9, cex.main=2.0, 
    cex.lab=1.9, cex.sub=1.0, cex.axis=1.7,
    #cex.lab=2.9,cex.sub=2.4,cex.axis=2.4,
    mar=c(4,5,2,1), 
    oma=c(2,3,2,1),font.main = 3)
plot(as.factor(names(par.sim)[21:29]),as.numeric(par.sim)[21:29],type="h",ylim=c(-1,1),axes=FALSE,lwd=3,
  ylab="",xlab="",main="Occupancy")
  axis(1, at=1:9, labels=1:9)
  axis(2)
  box()
for(j in 1:9)
 for(i in N.start1:N.end1)
  points(j+rnorm(1,0,0.05),par.mean[[i]][j+20])
plot(as.factor(names(par.sim)[30:38]),as.numeric(par.sim)[30:38],type="h",ylim=c(-1,1),axes=FALSE,lwd=3,
  ylab="",xlab="",main="Detection")
  axis(1, at=1:9, labels=1:9)
  axis(2)
  box()
for(j in 1:9)
 for(i in N.start1:N.end1)
  points(j+rnorm(1,0,0.05),par.mean[[i]][j+29])
mtext("Formation number", outer=T, padj=0, side=1, cex=6.5)
mtext("Formation-dependent random effect value", outer=T, padj=0.5, side=2, cex=6.5)
dev.off()

png("par_sim_03.png",height=2500,width=4000)
par(mfrow=c(3,1))
par(cex=2.9, cex.main=2.0, 
    cex.lab=1.9, cex.sub=1.0, cex.axis=1.66,
    #cex.lab=2.9,cex.sub=2.4,cex.axis=2.4,
    mar=c(4,5,2,1), 
    oma=c(2,3,2,1),font.main = 3)
plot(as.factor(names(par.sim)[39:47]),as.numeric(par.sim)[39:47],type="h",ylim=c(-3,3),axes=FALSE,lwd=3,
  ylab="",xlab="", main="Species 1")
  axis(1, at=1:9, labels=1:9)
  axis(2, at=seq(-3,3,2), labels=seq(-3,3,2) )
  box()
for(j in 1:9)
 for(i in N.start1:N.end1)
  points(j,par.mean[[i]][j+38])
plot(as.factor(names(par.sim)[48:56]),as.numeric(par.sim)[48:56],type="h",ylim=c(-3,3),axes=FALSE,lwd=3,
  ylab="",xlab="", main="Species 2")
  axis(1, at=1:9, labels=1:9)
  axis(2, at=seq(-3,3,2), labels=seq(-3,3,2) )
  box()
for(j in 1:9)
 for(i in N.start1:N.end1)
  points(j,par.mean[[i]][j+47])
plot(as.factor(names(par.sim)[57:65]),as.numeric(par.sim)[57:65],type="h",ylim=c(-3,3),axes=FALSE,lwd=3,
  ylab="",xlab="", main="Species 3")
  axis(1, at=1:9, labels=1:9)
  axis(2, at=seq(-3,3,2), labels=seq(-3,3,2) )
  box()
for(j in 1:9)
 for(i in N.start1:N.end1)
  points(j,par.mean[[i]][j+56])
mtext("Formation number", outer=T, padj=0, side=1, cex=6.5)
mtext("Species- and formation-dependent occupancy random effect value", outer=T, padj=0.5, side=2, cex=6.5)
dev.off()



png("par_sim_04.png",height=3000,width=4000)
par(mfrow=c(4,1))
par(cex=2.9, cex.main=2.0, 
    cex.lab=1.9, cex.sub=1.0, cex.axis=1.76,
    #cex.lab=2.9,cex.sub=2.4,cex.axis=2.4,
    mar=c(4,5,2,1), 
    oma=c(2,3,2,1),font.main = 3)
plot(as.factor(names(par.sim)[66:74]),as.numeric(par.sim)[66:74],type="h",ylim=c(-2,2),axes=FALSE,lwd=3,
  ylab="",xlab="", main="Species 1" )
  axis(1, at=1:9, labels=1:9)
  axis(2)
  box()

for(j in 1:9)
 for(i in N.start1:N.end1)
  points(j,par.mean[[i]][j+56+9])
plot(as.factor(names(par.sim)[75:83]),as.numeric(par.sim)[75:83],type="h",ylim=c(-2,2),axes=FALSE,lwd=3,
  ylab="",xlab="", main="Species 2" )
  axis(1, at=1:9, labels=1:9)
  axis(2)
  box()

for(j in 1:9)
 for(i in N.start1:N.end1)
  points(j,par.mean[[i]][j+56+2*9])
plot(as.factor(names(par.sim)[84:92]),as.numeric(par.sim)[84:92],type="h",ylim=c(-2,2),axes=FALSE,lwd=3,
  ylab="",xlab="", main="Species 3" )
  axis(1, at=1:9, labels=1:9)
  axis(2)
  box()

for(j in 1:9)
 for(i in N.start1:N.end1)
  points(j,par.mean[[i]][j+56+3*9])
plot(as.factor(names(par.sim)[93:101]),as.numeric(par.sim)[93:101],type="h",ylim=c(-2,2),axes=FALSE,lwd=3,
  ylab="",xlab="", main="Superspecies")
  axis(1, at=1:9, labels=1:9)
  axis(2)
  box()
for(j in 1:9)
 for(i in N.start1:N.end1)
  points(j,par.mean[[i]][j+56+4*9])
mtext("Formation number", outer=T, padj=0, side=1, cex=7.5)
mtext("Species- and formation-dependent detection random effect value", outer=T, padj=0.5, side=2, cex=7.5)
dev.off()





# Occupancy
const.occu=rep(1,F)%*%t(par.sim[const.occu.fi+1:(S-1)-1])
f.occu=matrix(rep(par.sim[f.occu.fi+1:F-1],N.sp-1),ncol=N.sp-1)
occu.sf.mean=t(matrix(c(ilogit(matrix(par.sim[sf.occu.fi+1:((S-1)*F)-1],ncol=N.sp-1)+const.occu+f.occu),rep(1,F)), ncol=N.sp))

occu.sf.sim=array(NA,c(N.sim,N.sp,Nf))
psistar.sim=array(NA,c(N.sim,N.sp,Nf))
p.occu=array(NA,c(N.sp-1,Nf))
bias.occu=array(NA,c(N.sp-1,Nf))
sd.occu=array(NA,c(N.sp-1,Nf))

#Fig. 3 style
png("occupancy_sim.png",height=1600,width=2400)
par(cex=9.9, cex.main=5.9, 
    cex.lab=9.9,cex.sub=6.9,cex.axis=4.9,
    mar=c(10,10,9,1), 
    oma=c(3,3,3,1),font.main = 3)
par(mfrow=c(1,3))
species=c("Species 1\n", "Species 2\n" , "Species 3\n" ,"Superspecies\n")
for(s in 1:(N.sp-1))
{
 plot(1:F, occu.sf.mean[s,], type="b", ylim=c(0,1.0),
   xlab="", ylab="", main=paste(species[s]),lwd=8, 
   axes=FALSE, pch=20)
  axis(1, at=c(3,6,9), labels=c("3","6","9"),padj=1)
  if(s==1)
    axis(2, at=c(0,0.5,1),labels=c(0,0.5,1))
  box()

  for(i in 1:N.sim)
  {
    const.occu.sim=rep(1,F)%*%t(par.mean[[i]][const.occu.fi+
                                              1:(S-1)-1])
    f.occu.sim=matrix(rep(par.mean[[i]][f.occu.fi+
                          1:F-1],N.sp-1),ncol=N.sp-1)
    occu.sf.sim[i,,]=
      t(matrix(c(ilogit(matrix(par.mean[[i]][sf.occu.fi+
      1:((S-1)*F)-1],ncol=N.sp-1)+const.occu+f.occu),
      rep(1,F)),ncol=N.sp))

    psistar.sim[i,,]=
      t(matrix(c(ilogit(matrix(par.mean[[i]][sf.occu.fi+
      1:((S-1)*F)-1],ncol=N.sp-1)+const.occu+f.occu),
      rep(1,F)),ncol=N.sp))


   if(i>=N.start1 & i<=N.end1)
     points(jitter(1:F,factor=0.5), occu.sf.sim[i,s,],lwd=2,col="grey")
  }

  for(f in 1:Nf)
  {
    p.occu[s,f]=as.numeric(t.test(occu.sf.sim[,s,f],
                           mu=occu.sf.mean[s,f])$p.value)
    bias.occu[s,f]=occu.sf.mean[s,f]-mean(occu.sf.sim[,s,f])
    sd.occu[s,f]=sd(occu.sf.sim[,s,f])
  }
}
mtext("Formation number", outer=T, padj=0, side=1, cex=4.5)
mtext("Occupancy probability", outer=T, padj=0.5, side=2, cex=4.5)
dev.off()

sum(p.occu<0.05/length(p.occu))
# 0

sum(sd.occu<abs(bias.occu))
# 0


mean(abs(bias.occu)/sd.occu)




# Detection:
const.bin=rep(1,F)%*%t(par.sim[const.bin.fi+1:S-1])
f.bin=matrix(rep(par.sim[f.bin.fi+1:F-1],N.sp),ncol=N.sp)
bin.sf.mean=t(ilogit(const.bin+f.bin+matrix(par.sim[sf.bin.fi+1:(S*F)-1],ncol=N.sp)))

bin.sf.sim=array(NA,c(N.sim,N.sp,Nf))
pstar.sim=array(NA,c(N.sim,N.sp,Nf))
p.bin=array(NA,c(N.sp,Nf))
bias.bin=array(NA,c(N.sp,Nf))
sd.bin=array(NA,c(N.sp,Nf))

png("detection_sim.png",height=1600,width=2400)
par(cex=5.9, cex.main=5.9, 
    cex.lab=9.9,cex.sub=6.9,cex.axis=4.9,
    mar=c(10,10,9,1), 
    oma=c(3,3,3,1),font.main = 3)
par(mfrow=c(2,2))
species=c("Species 1", "Species 2" , "Species 3" ,"Superspecies")
for(s in 1:N.sp)
{
 plot(1:F, bin.sf.mean[s,], type="b", ylim=c(0,1.1),
   xlab="", ylab="", main=paste(species[s]),lwd=8, 
   axes=FALSE, pch=20)
  if(s==3 | s==4)
    axis(1, at=c(3,6,9), labels=c(3,6,9), padj=1)
  if(s==1 | s==3)
    axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
  box()

  for(i in 1:N.sim)
  {
    const.bin.sim=rep(1,F)%*%t(par.mean[[i]][const.bin.fi+
                                             1:S-1])
    f.bin.sim=matrix(rep(par.mean[[i]][f.bin.fi+
                                       1:F-1],N.sp),ncol=N.sp)
    bin.sf.sim[i,,]=
      t(ilogit(matrix(par.mean[[i]][sf.bin.fi+
                                    1:(S*F)-1],
        ncol=N.sp)+const.bin+f.bin))

    pstar.sim[i,,]=
      t(ilogit(matrix(par.mean[[i]][sf.bin.fi+
                                    1:(S*F)-1],
        ncol=N.sp)+const.bin))


   if(i>=N.start1 & i<=N.end1)
     points(jitter(1:F,factor=0.5), bin.sf.sim[i,s,],lwd=2,col="grey")
  }

  for(f in 1:Nf)
  {
    p.bin[s,f]=as.numeric(t.test(bin.sf.sim[,s,f],
                           mu=bin.sf.mean[s,f])$p.value)
    bias.bin[s,f]=bin.sf.mean[s,f]-mean(bin.sf.sim[,s,f])
    sd.bin[s,f]=sd(bin.sf.sim[,s,f])
  }
}
mtext("Formation number", outer=T, padj=0, side=1, cex=4.5)
mtext("Detection probability", outer=T, padj=0.5, side=2, cex=4.5)
dev.off()


sum(p.bin<0.05/length(p.bin))
# 7

sum(sd.bin<abs(bias.bin))
# 2




# Relative abundance:
N.start1=1
N.end1=95
N.start2=96
N.end2=100

abundance.given.occupancy=-log(1-bin.formation.species)
abundance=abundance.given.occupancy*occu.formation.species2
rel.abundance=abundance
for(f in 1:Nf)
  rel.abundance[f,]=abundance[f,]/sum(abundance[f,])
rel.abundance=t(rel.abundance)

ra.sf.sim=array(NA,c(N.sim,N.sp,Nf))
p.ra=array(NA,c(N.sp,Nf))
bias.ra=array(NA,c(N.sp,Nf))
sd.ra=array(NA,c(N.sp,Nf))


png("R_sim.png",height=1600,width=2400)
par(cex=5.9, cex.main=5.9, 
    cex.lab=9.9,cex.sub=6.9,cex.axis=4.9,
    mar=c(10,10,9,1), 
    oma=c(3,3,3,1),font.main = 3)
par(mfrow=c(2,2))
species=c("Species 1", "Species 2" , "Species 3" ,"Superspecies")
use.max=c(0.2,0.2,0.1,0.9)
use.min=c(0.0,0.0,0.0,0.7)
for(s in 1:N.sp)
{
 plot(1:F, rel.abundance[s,], type="b", 
   ylim=c(use.min[s],use.max[s]),
   xlab="", ylab="", main=paste(species[s]),lwd=8, 
   axes=FALSE, pch=20, col="grey")
  if(s==3 | s==4)
    axis(1, at=c(3,6,9), labels=c(3,6,9),padj=1)
  if(s==4)
    axis(2, at=c(0.6,0.7,0.8,0.9), labels=c(0.6,0.7,0.8,0.9))
  if(s==1)
    axis(2, at=c(0,0.1,0.2, 0.3), labels=c(0,0.1,0.2,0.3))
  if(s==3)
    axis(2, at=c(0,0.05,0.1), labels=c(0,0.05,0.1))
  box()

  for(i in N.start1:N.end1)
  {
    abundance.given.occu.sim=-log(1-pstar.sim[i,,])
    abundance.sim=abundance.given.occu.sim*psistar.sim[i,,]
    for(f in 1:Nf)
      ra.sf.sim[i,,f]=abundance.sim[,f]/sum(abundance.sim[,f])

    points(jitter(1:F,factor=0.25), ra.sf.sim[i,s,], lty=3,
      lwd=2,col="black")
  }

  for(i in N.start2:N.end2)
  {
    abundance.given.occu.sim=-log(1-bin.sf.sim[i,,])
    abundance.sim=abundance.given.occu.sim*occu.sf.sim[i,,]
    for(f in 1:Nf)
      ra.sf.sim[i,,f]=abundance.sim[,f]/sum(abundance.sim[,f])

    gr=(i+2-N.start2)/7
    show(gr)
    lines(jitter(1:F,factor=0.5), ra.sf.sim[i,s,], 
          lwd=2,col=grey(gr),lty=i-N.start2+2)
  }

  for(f in 1:Nf)
  {
    p.ra[s,f]=as.numeric(t.test(ra.sf.sim[,s,f],
                           mu=rel.abundance[s,f])$p.value)
    bias.ra[s,f]=rel.abundance[s,f]-mean(ra.sf.sim[,s,f])
    sd.ra[s,f]=sd(ra.sf.sim[,s,f])
  }
}
mtext("Formation number", outer=T, padj=0, side=1, cex=4.5)
mtext("Relative species abundance, R", outer=T, padj=0.5, side=2, cex=4.5)
dev.off()


sum(p.ra<0.05/length(p.ra))
# 16

sum(sd.ra<abs(bias.ra))
# 2



