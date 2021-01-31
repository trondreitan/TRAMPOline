
##########################
# Fetch and process data:
##########################

y=read.csv(sprintf("sim_s%03d.csv",sim.nr),sep=";",header=T)

N.sp=dim(y)[2]-which(names(y)=="num.with.colonies")
N.sites=dim(y)[1]
N.forms=length(unique(y$Formation_name))

logit=function(x) log(x/(1-x))
ilogit=function(x) 1/(1+exp(-x))


#####################################################

#Shortform for essential data sizes:
# Number of species:
S=N.sp
# Number of sites:
N=N.sites
# Number of formations:
Nf=N.forms



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
