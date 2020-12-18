
model="f_sf_modelselect"

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
N.runs=50

# all
for(i in 1:N.runs)
{
  show(sprintf("%s_%02d.Rdata",model,i+N.start-1))
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
  {
    res[[j+k-start+1]]=mcmc(Fit3[[j]]$Monitor)
  }
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

for(i in 1:length(par0))
{
  nn=names(par0)[i]
  samples=res2[,i]

  if(substr(nn,1,6)=="logit.")
  {
    nn=substr(nn,7,nchar(nn))
    samples=ilogit(samples)
  }
  if(substr(nn,1,4)=="log.")
  {
    nn=substr(nn,5,nchar(nn))
    samples=exp(samples)
  }

  m=mean(samples)
  s=sd(samples)
  md=median(samples)
  q1=as.numeric(quantile(samples,0.25))
  q2=as.numeric(quantile(samples,0.75))
  lower=as.numeric(quantile(samples,0.025))
  upper=as.numeric(quantile(samples,0.975))
  prob.pos=mean(samples>0)

  if(i==1)
    print.srcref(sprintf("%20s %7s %7s %7s %7s %7s %7s %7s %7s",
      "Name","mean","sd","q1","median","q3","lower","upper", "prob.pos"))
  
  print.srcref(sprintf("%20s %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f",
    nn,m,s,q1,md,q2,lower,upper,prob.pos))
    
}

#                Name    mean      sd      q1  median      q3   lower   upper  prob.pos
# ...
#         use.o18.bin -0.1032  0.9931 -0.7673 -0.1449  0.5614 -2.0049  1.8959  0.4351
#        use.o18.occu -0.3970  0.9170 -0.9659 -0.4293  0.0019 -2.1276  1.6497  0.2504
#      global.o18.bin -0.1579  0.9864 -0.8111 -0.2072  0.4901 -2.0337  1.8599  0.4017
#     global.o18.occu  0.0253  0.9996 -0.6498  0.0381  0.6992 -1.9425  1.9715  0.5157
#          use.ca.bin -0.7092  0.7059 -1.1160 -0.6285 -0.2611 -2.2210  0.7604  0.0560
#         use.ca.occu -0.6630  0.7485 -1.0960 -0.6040 -0.2289 -2.2055  1.0491  0.0843
#       global.ca.bin  0.0200  1.0005 -0.6537  0.0306  0.6935 -1.9484  1.9740  0.5124
#      global.ca.occu  0.0316  0.9985 -0.6412  0.0484  0.7045 -1.9399  1.9711  0.5201

nn=names(par0)
u1=which(nn=="use.o18.bin")
u2=which(nn=="use.o18.occu")
u3=which(nn=="use.ca.bin")
u4=which(nn=="use.ca.occu")

s0=mean(res2[,u1]<=0 & res2[,u2]<=0 & res2[,u3]<=0 & res2[,u4]<=0)
s1=mean(res2[,u1]<=0 & res2[,u2]<=0 & res2[,u3]<=0 & res2[,u4]>0)
s2=mean(res2[,u1]<=0 & res2[,u2]<=0 & res2[,u3]>0 & res2[,u4]<=0)
s3=mean(res2[,u1]<=0 & res2[,u2]<=0 & res2[,u3]>0 & res2[,u4]>0)
s4=mean(res2[,u1]<=0 & res2[,u2]>0 & res2[,u3]<=0 & res2[,u4]<=0)
s5=mean(res2[,u1]<=0 & res2[,u2]>0 & res2[,u3]<=0 & res2[,u4]>0)
s6=mean(res2[,u1]<=0 & res2[,u2]>0 & res2[,u3]>0 & res2[,u4]<=0)
s7=mean(res2[,u1]<=0 & res2[,u2]>0 & res2[,u3]>0 & res2[,u4]>0)
s8=mean(res2[,u1]>0 & res2[,u2]<=0 & res2[,u3]<=0 & res2[,u4]<=0)
s9=mean(res2[,u1]>0 & res2[,u2]<=0 & res2[,u3]<=0 & res2[,u4]>0)
s10=mean(res2[,u1]>0 & res2[,u2]<=0 & res2[,u3]>0 & res2[,u4]<=0)
s11=mean(res2[,u1]>0 & res2[,u2]<=0 & res2[,u3]>0 & res2[,u4]>0)
s12=mean(res2[,u1]>0 & res2[,u2]>0 & res2[,u3]<=0 & res2[,u4]<=0)
s13=mean(res2[,u1]>0 & res2[,u2]>0 & res2[,u3]<=0 & res2[,u4]>0)
s14=mean(res2[,u1]>0 & res2[,u2]>0 & res2[,u3]>0 & res2[,u4]<=0)
s15=mean(res2[,u1]>0 & res2[,u2]>0 & res2[,u3]>0 & res2[,u4]>0)

s=c(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15)
s
# [1] 0.34567000 0.03615750 0.02455750 0.00203375 0.13487250 0.01274625
# [7] 0.00819750 0.00063750 0.30054625 0.02465375 0.01489500 0.00113250
#[13] 0.08268750 0.00663000 0.00423875 0.00034375

# Most likely, except null model, s8, o18.bin. 

p.null=s[1]
p.else=mean(s[2:16])

# If o18, what's the status for global vs species-wise?
i8.2=which(nn=="global.o18.bin")
mean(res2[,u1]>0 & res2[,u2]<=0 & res2[,u3]<=0 & res2[,u4]<=0 & res2[,i8.2]<=0)
#0.2271738
# Species-wise o18-dependency preferred.

