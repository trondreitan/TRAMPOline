# Read mean parameters for full model:
load("../par_mean_f_sf.RData")

# Some extra info:
load("plotting_variables.Rdata")

# betabinomial distribution:
dbetabin=function(x,n,p,s)
  exp(lgamma(n+1)-lgamma(x+1)-lgamma(n-x+1) +
      lgamma(x+p/s)+lgamma(n-x+(1-p)/s)-lgamma(n+1/s) +
      lgamma(1/s)-lgamma(p/s)-lgamma((1-p)/s))
    
rbetabin=function(N,n,p,s)
  rbinom(N,size=n,prob=rbeta(N,p/s,(1-p)/s))

# Read data:
y=read.csv("../allsamples_with_counts_and_metainfo.csv",sep=";",header=T)

# Useful definitions:
N.sp=dim(y)[2]-which(names(y)=="num.with.colonies")
N.sites=dim(y)[1]
N.forms=length(unique(y$Formation_name))

# Useful functions:
logit=function(x) log(x/(1-x))
ilogit=function(x) 1/(1+exp(-x))




# Simulate just as the data to start with:
sim.data=function(nshells, nsites, par, data)
{
  y.new=data[1:min(dim(data)[1],nsites*Nf),]

  for(f in 1:Nf)
  {
    #show(f)
    j=which(data$form.nr==f)
    y.new[(f-1)*nsites+1:nsites,]=data[sample(j,nsites,replace=TRUE),]
    y.new$SAMPLE_ID[(f-1)*nsites+1:nsites]=sprintf("%d",(f-1)*nsites+1:nsites)
  }

  S=4
  s.fi=which(names(data)=="Superspecies")-S+1
  const.bin.fi=which(names(par)=="logit.const.bin1")
  const.occu.fi=which(names(par)=="logit.const.occu1")
  overdisp.fi=which(names(par)=="log.betabin.s1")
  f.bin.fi=which(names(par)=="f.bin1")
  f.occu.fi=which(names(par)=="f.occu1")
  sf.bin.fi=which(names(par)=="sf.bin1")
  sf.occu.fi=which(names(par)=="sf.occu1")

  for(i in 1:dim(y.new)[1])
   for(s in 1:S)
    { 
     f=y.new$form.nr[i]
     p.bin=ilogit(par[const.bin.fi+s-1]+par[f.bin.fi+f-1]+
                  par[sf.bin.fi+N.forms*(s-1)+f-1])
     if(s<S)
       p.occu=ilogit(par[const.occu.fi+s-1]+par[f.occu.fi+f-1]+
                    par[sf.occu.fi+N.forms*(s-1)+f-1])
     #p .bin=bin.sf.mean[s,f]
     #if(s<S)
     #  p.occu=occu.sf.mean[s,f]
     if(s==S)
       p.occu=1

     overdisp=exp(par[overdisp.fi+s-1])

     rr=runif(1)
     y.new[i,s.fi+s-1]=(rr<p.occu)*rbetabin(1,nshells,p.bin,overdisp)
    }

  y.new$Total=nshells
  y.new
}


n.shell=60
n.site=10
N.sim=200
F=9

par=par.sim=par.mean

S=4
s.fi=which(names(y)=="Superspecies")-S+1
const.bin.fi=which(names(par)=="logit.const.bin1")
const.occu.fi=which(names(par)=="logit.const.occu1")
overdisp.fi=which(names(par)=="log.betabin.s1")
f.bin.fi=which(names(par)=="f.bin1")
f.occu.fi=which(names(par)=="f.occu1")
sf.bin.fi=which(names(par)=="sf.bin1")
sf.occu.fi=which(names(par)=="sf.occu1")

Psi=array(0,c(4,9))
Psi[1,]=ilogit(logit(0.5)+seq(-2,2,length.out=F))
Psi[2,]=rep(0.9,9)
Psi[3,]=rep(0.5,9)
Psi[4,]=rep(1,9)

A=array(0,c(4,9))
A[1,]=0.2*Psi[1,]
A[2,]=0.2*seq(1.0,0.1,length.out=9)
A[3,]=0.2*rep(0.5,9)
A[4,]=rep(2,9)

lambda=A/Psi
p=1-exp(-lambda)


par.sim[overdisp.fi+1:S-1]=log(rep(0.05,S))

par.sim[f.bin.fi+1:F-1]=0
par.sim[f.occu.fi+1:F-1]=0

par.sim[const.bin.fi+1:S-1]=apply(logit(p),1,median)
par.sim[const.occu.fi+1:(S-1)-1]=apply(logit(Psi[1:3,]),1,median)

# First species will have variation in occupancy:
par.sim[sf.occu.fi+1:((S-1)*F)-1]=c(logit(Psi[1,])-par.sim[const.occu.fi],
     rep(0,2*F))

# Second species will have variation in detection:
par.sim[sf.bin.fi+1:(S*F)-1]=c(rep(0,F),logit(p[2,])-par.sim[const.bin.fi+1],
  rep(0,2*F))

par.sim[12:13]=log(0.001)
par.sim[14]=sd(par.sim[sf.occu.fi+1:F-1])
par.sim[15:16]=log(0.001)
par.sim[17]=0.001
par.sim[18]=sd(par.sim[sf.bin.fi+F+1:F-1])
par.sim[19:20]=log(0.001)


save.image(file="par_sim.RData")

const.bin=rep(1,F)%*%t(par.sim[const.bin.fi+1:S-1])
f.bin=matrix(rep(par.sim[f.bin.fi+1:F-1],N.sp),ncol=N.sp)
bin.formation.species=ilogit(const.bin+f.bin+matrix(par.sim[sf.bin.fi+1:(S*F)-1],ncol=N.sp))
abundance.given.occupancy=-log(1-bin.formation.species)

const.occu=rep(1,F)%*%t(par.sim[const.occu.fi+1:(S-1)-1])
f.occu=matrix(rep(par.sim[f.occu.fi+1:F-1],N.sp-1),ncol=N.sp-1)
occu.formation.species2=matrix(c(ilogit(matrix(par.sim[sf.occu.fi+1:((S-1)*F)-1],ncol=N.sp-1)+const.occu+f.occu),rep(1,F)), ncol=N.sp)
abundance=abundance.given.occupancy*occu.formation.species2

rel.abundance=abundance
for(i in 1:Nf)
  rel.abundance[i,]=abundance[i,]/sum(abundance[i,])



for(i in 1:N.sim)
{
  y.new=sim.data(n.shell, n.site, par.sim, y)

  write.table(y.new , file=sprintf("sim_s%03d.csv",
      i), sep=";", row.names=FALSE, col.names=TRUE)
}
