
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
N.runs=146

# all
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

  if(i==1)
    print.srcref(sprintf("%20s %7s %7s %7s %7s %7s %7s %7s",
      "Name","mean","sd","q1","median","q3","lower","upper"))
  
  print.srcref(sprintf("%20s %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f",
    nn,m,s,q1,md,q2,lower,upper))
    
}
#                Name    mean      sd      q1  median      q3   lower   upper
#          betabin.s1  0.0931  0.0266  0.0742  0.0901  0.1086  0.0502  0.1540
#          betabin.s2  0.0516  0.0229  0.0353  0.0476  0.0635  0.0184  0.1071
#          betabin.s3  0.0401  0.0185  0.0269  0.0374  0.0504  0.0123  0.0839
#          betabin.s4  0.0727  0.0138  0.0630  0.0716  0.0812  0.0490  0.1031
#          const.bin1  0.0543  0.0256  0.0373  0.0496  0.0655  0.0200  0.1167
#          const.bin2  0.0285  0.0373  0.0065  0.0149  0.0377  0.0010  0.1169
#          const.bin3  0.0207  0.0144  0.0115  0.0174  0.0261  0.0044  0.0560
#          const.bin4  0.4204  0.0796  0.3702  0.4146  0.4637  0.2766  0.6009
#         const.occu1  0.9722  0.0443  0.9623  0.9871  0.9974  0.8728  1.0000
#         const.occu2  0.7300  0.3004  0.6087  0.8724  0.9373  0.0214  0.9969
#         const.occu3  0.8764  0.1735  0.8571  0.9350  0.9801  0.3018  0.9997
#           sd.f.occu  0.6361  1.6066  0.0190  0.1332  0.6462  0.0005  3.8903
#            sd.f.bin  0.8545  0.3421  0.6212  0.7782  1.0059  0.4229  1.7185
#         sd.sf.occu1  0.6730  1.6720  0.0162  0.1222  0.6242  0.0002  4.5130
#         sd.sf.occu2  5.6211 34.6212  0.0436  0.4732  3.6859  0.0005 33.3192
#         sd.sf.occu3  4.2798 25.3213  0.0300  0.2691  1.7213  0.0004 28.1846
#          sd.sf.bin1  0.8070  0.4797  0.5281  0.7990  1.0769  0.0022  1.8384
#          sd.sf.bin2  2.0488  1.4841  1.1384  1.8605  2.6929  0.0060  5.5672
#          sd.sf.bin3  1.0520  0.6661  0.6559  0.9872  1.3771  0.0056  2.5926
#          sd.sf.bin4  0.2305  0.2971  0.0130  0.0881  0.3758  0.0002  0.9938
#             f.occu1  0.1482  1.3664 -0.0468  0.0005  0.0976 -1.4683  2.6411
#             f.occu2  0.1220  1.5619 -0.0621  0.0001  0.0761 -1.7974  2.8345
#             f.occu3  0.2087  1.6324 -0.0451  0.0005  0.1017 -1.2806  3.0105
#             f.occu4 -0.0658  1.4497 -0.0971 -0.0004  0.0476 -2.7487  1.7768
#             f.occu5  0.0755  1.5298 -0.0708 -0.0000  0.0663 -1.8699  2.3689
#             f.occu6  0.0562  1.5069 -0.0684 -0.0000  0.0683 -1.8984  2.1433
#             f.occu7 -0.1487  1.3760 -0.1129 -0.0007  0.0411 -3.2320  1.4962
#             f.occu8  0.3939  1.8933 -0.0251  0.0022  0.1676 -0.8199  4.3373
#             f.occu9  0.2997  1.3343 -0.0231  0.0026  0.1735 -0.8797  3.4974
#              f.bin1  0.6437  0.4840  0.3222  0.6178  0.9375 -0.2427  1.6820
#              f.bin2 -1.0262  0.4355 -1.2574 -0.9881 -0.7491 -2.0326 -0.2586
#              f.bin3 -0.0894  0.4098 -0.3293 -0.0735  0.1686 -0.9648  0.6889
#              f.bin4 -1.2278  0.4691 -1.4729 -1.1834 -0.9296 -2.3212 -0.4186
#              f.bin5 -0.3383  0.4793 -0.5625 -0.2613 -0.0248 -1.5334  0.4090
#	      f.bin6  0.2413  0.4314  0.0133  0.2880  0.5228 -0.7707  0.9734
#              f.bin7  0.1742  0.4713 -0.1164  0.1440  0.4272 -0.6866  1.2524
#              f.bin8  0.0696  0.3856 -0.1399  0.0956  0.3114 -0.8005  0.7700
#              f.bin9  0.9544  0.4885  0.6362  0.8681  1.1898  0.1735  2.1357
#            sf.occu1  0.1654  1.7529 -0.0489  0.0001  0.0803 -1.7197  3.0306
#            sf.occu2  0.0981  1.6972 -0.0594  0.0000  0.0666 -2.1181  2.7299
#            sf.occu3  0.1789  1.6959 -0.0469  0.0002  0.0833 -1.6258  3.1199
#            sf.occu4  0.0252  1.5816 -0.0674 -0.0000  0.0588 -2.6218  2.5170
#            sf.occu5 -0.0651  1.5349 -0.0805 -0.0001  0.0497 -3.3117  2.1710
#            sf.occu6  0.1643  1.5565 -0.0482  0.0001  0.0815 -1.6736  3.0582
#            sf.occu7  0.2169  1.5013 -0.0391  0.0004  0.0979 -1.4275  3.4230
#            sf.occu8  0.1339  1.5860 -0.0517  0.0001  0.0757 -1.8096  2.9111
#            sf.occu9  0.1071  1.4659 -0.0568  0.0000  0.0685 -1.8408  2.6978
#           sf.occu10 -3.2510 31.0404 -1.3527 -0.0125  0.0342 -23.6799  4.3870
#           sf.occu11  4.0479 31.6824 -0.0209  0.0211  1.5061 -1.7376 26.7345
#           sf.occu12  3.7159 32.0621 -0.0307  0.0138  1.1731 -2.0638 27.2732
#           sf.occu13 -2.7490 21.9808 -1.3363 -0.0115  0.0368 -24.0779  5.0646
#           sf.occu14  3.2528 43.2227 -0.2039  0.0001  0.2690 -4.2212 21.1035
#           sf.occu15 -4.1602 30.1106 -2.3642 -0.0283  0.0157 -27.7435  2.4459
#           sf.occu16 -3.6098 30.9266 -1.7579 -0.0190  0.0236 -25.7161  3.3135
#           sf.occu17  4.4526 32.3316 -0.0078  0.0487  2.6128 -1.0964 28.1602
#           sf.occu18  3.1234 24.9420 -0.0182  0.0244  1.5805 -1.9925 17.9477
#           sf.occu19  3.1295 23.1667 -0.0288  0.0060  0.5219 -1.6609 21.7665
#           sf.occu20 -2.8938 25.3471 -0.5324 -0.0047  0.0347 -20.7099  2.6648
#           sf.occu21  2.4946 27.0281 -0.1355  0.0000  0.1433 -3.9864 15.5693
#           sf.occu22  2.5383 22.2236 -0.0815  0.0007  0.2308 -3.0496 19.1899
#           sf.occu23  2.7013 22.8653 -0.0556  0.0019  0.3192 -2.2102 19.5197
#           sf.occu24  2.7433 22.3703 -0.0465  0.0028  0.3631 -2.0485 20.2074
#           sf.occu25 -3.0475 21.3298 -0.6902 -0.0073  0.0254 -22.2616  2.1512
#           sf.occu26  2.6095 19.4250 -0.0522  0.0023  0.3354 -2.1530 19.3505
#           sf.occu27  2.5473 22.9607 -0.0366  0.0042  0.4222 -2.1074 17.3346
#             sf.bin1  0.5649  0.6626  0.0220  0.4871  0.9957 -0.5348  2.0034
#             sf.bin2 -0.3235  0.5670 -0.6577 -0.2323  0.0121 -1.5994  0.6748
#             sf.bin3  0.2917  0.5206 -0.0104  0.2248  0.6144 -0.6808  1.4113
#             sf.bin4 -0.5483  0.7067 -0.9548 -0.4265 -0.0073 -2.1964  0.5462
#             sf.bin5 -0.9133  0.7237 -1.3816 -0.8759 -0.3331 -2.4706  0.1189
#             sf.bin6 -0.3424  0.5023 -0.6579 -0.2840  0.0000 -1.4264  0.5570
#             sf.bin7  0.8191  0.6334  0.3011  0.8090  1.2529 -0.1261  2.1138
#             sf.bin8 -0.3547  0.4680 -0.6497 -0.3070 -0.0045 -1.3613  0.4751
#             sf.bin9  0.6438  0.5398  0.2035  0.6689  1.0129 -0.2881  1.6936
#	    sf.bin10 -1.2599  2.2559 -2.2072 -0.8120  0.0063 -6.6478  2.0546
#            sf.bin11  0.6455  1.1949 -0.0744  0.4574  1.3421 -1.4077  3.3293
#            sf.bin12  0.2056  1.1822 -0.4705  0.0343  0.8458 -2.0011  2.8234
#            sf.bin13 -1.2975  2.2098 -2.1963 -0.8240  0.0012 -6.6477  1.7578
#            sf.bin14 -0.7240  1.2685 -1.4539 -0.5848  0.0077 -3.4911  1.6473
#            sf.bin15 -1.7678  2.3078 -2.8233 -1.4121 -0.1184 -7.1996  1.4842
#            sf.bin16 -1.4695  2.2729 -2.4524 -1.0438 -0.0058 -6.8698  1.7491
#            sf.bin17  1.8342  1.3294  0.8148  1.8269  2.6922 -0.1417  4.6006
#            sf.bin18  1.9355  1.3940  0.8531  1.9527  2.8569 -0.1842  4.7890
#            sf.bin19  0.6645  0.8061  0.0280  0.5673  1.1716 -0.6689  2.4375
#            sf.bin20 -0.9157  1.1523 -1.4797 -0.7107 -0.0688 -3.6682  0.7108
#            sf.bin21 -0.4973  0.8043 -0.9466 -0.3771  0.0030 -2.3348  0.8890
#            sf.bin22 -0.2067  0.8023 -0.6343 -0.0971  0.2266 -1.9893  1.3177
#            sf.bin23  0.1108  0.6491 -0.2523  0.0393  0.4669 -1.1488  1.5211
#            sf.bin24 -0.1053  0.6346 -0.4740 -0.0599  0.2281 -1.3931  1.2170
#            sf.bin25 -0.9238  1.1691 -1.5106 -0.7322 -0.0636 -3.6831  0.7864
#            sf.bin26  0.0738  0.5963 -0.2697  0.0108  0.3952 -1.0567  1.3821
#            sf.bin27  1.4041  0.8151  0.8908  1.4696  1.9416 -0.0145  2.9527
#            sf.bin28 -0.0650  0.3046 -0.0804 -0.0008  0.0216 -0.9206  0.4548
#            sf.bin29  0.0282  0.2748 -0.0338  0.0001  0.0537 -0.5284  0.7453
#            sf.bin30  0.0051  0.2530 -0.0385  0.0000  0.0453 -0.5780  0.6048
#            sf.bin31  0.0191  0.2932 -0.0420 -0.0000  0.0443 -0.5824  0.7787
#            sf.bin32  0.1518  0.3447 -0.0051  0.0069  0.1935 -0.2396  1.1452
#            sf.bin33  0.1313  0.2957 -0.0044  0.0079  0.1863 -0.2282  0.9665
#            sf.bin34 -0.0908  0.3378 -0.0984 -0.0013  0.0173 -1.0729  0.4235
#            sf.bin35  0.0425  0.2393 -0.0210  0.0008  0.0743 -0.4198  0.6714
#            sf.bin36 -0.2302  0.4262 -0.3601 -0.0133  0.0022 -1.3499  0.1910


i7=which(substr(names(par0),1,12)=="logit.const.")

png("scatter.png",height=1200,width=1600)
par(mfrow=c(2,2),cex=2)
plot(res2[,i7[1]],res2[,i7[5]],pch=".",
 xlab="logit(const.detection,sp1)",ylab="logit(const.occupancy.sp1")
plot(res2[,i7[2]],res2[,i7[6]],pch=".",
 xlab="logit(const.detection,sp2)",ylab="logit(const.occupancy.sp2")
plot(res2[,i7[3]],res2[,i7[7]],pch=".",
 xlab="logit(const.detection,sp3)",ylab="logit(const.occupancy.sp3")
dev.off()

# Look at species2
x=seq(min(res2[,i7[2]]),max(res2[,i7[2]]),0.5)
nx=length(x)
y=seq(min(res2[,i7[6]]),max(res2[,i7[6]]),0.5)
ny=length(y)
A=array(0,c(nx,ny))
for(i in 1:nx)
  for(j in 1:ny)
    A[i,j]=mean(res2[,i7[2]]>=x[i] & res2[,i7[2]]<x[i]+0.1 &
      res2[,i7[6]]>=y[j] & res2[,i7[6]]<y[j]+0.1)
logA=log(A)
logA[is.infinite(logA)]=NA

png("heatmap2.png",height=1200,width=1600)
par(cex=2.4)
image(x,y,logA,
 xlab="logit(const.detection,sp2)",ylab="logit(const.occupancy.sp2")
dev.off()


c(cor(res2[,i7[1]],res2[,i7[5]]), cor(res2[,i7[2]],res2[,i7[6]]),
 cor(res2[,i7[3]],res2[,i7[7]]))
# -0.02042208 -0.44218690 -0.17245282

library(bayesplot)
mcmc_hex(res2, pars=c("logit.const.bin2","logit.const.occu2"))

# Theoretical entropy: 0.5*log(2*pi*exp(1)*sigma^2)
ent.norm=function(sigma) 0.5*log2(2*pi)+0.5+log(sigma)

for(i in 1:length(par0))
{
  nn=names(par0)[i]
  samples=res2[,i]
  m=0; s=0
  if(substr(nn,1,13)=="log.betabin.s")
    { m=hyper$log.betabin.s.mu; s= hyper$log.betabin.s.sd }
  if(substr(nn,1,11)=="logit.const")
    { m=hyper$logit.mu; s= hyper$logit.sd }
  if(substr(nn,1,6)=="log.sd")
    { m=hyper$logsd.mu; s= hyper$logsd.sd }
  if(substr(nn,1,6)=="log.dt")
    { m=hyper$log.dt.mu; s= hyper$log.dt.sd }
  if(substr(nn,1,10)=="o18.lwidth")
    { m=hyper$o18.lwidth.mu; s= hyper$o18.lwidth.sd }
  if(substr(nn,1,9)=="o18.optim")
    { m=hyper$o18.optim.mu; s= hyper$o18.optim.sd }
  
  if(s>0)
  {
    par(ask=T)
    #readline(prompt="Press <Enter> to continue")
    
    xx=seq(m-5*s,m+5*s,s/20)
    d=density(samples, from=min(xx),to=max(xx))
    plot(xx,dnorm(xx,m,s),type="l",ylim=c(0,max(max(dnorm(xx,m,s)),max(d$y))),
         xlab=nn,ylab="prob.density")
    lines(d,col="red")

    fname=sprintf("%s_%s.png",model,nn)
    png(sprintf("hists/%s",fname),height=800,width=1000)
    par(cex=1.4)
    plot(xx,dnorm(xx,m,s),type="l",ylim=c(0,max(max(dnorm(xx,m,s)),max(d$y))),
         xlab=nn,ylab="prob.density")
    lines(d,col="red")
    dev.off()

    dx=d$x[10]-d$x[9]
    ent.new=ent.norm(sd(samples))   # sum(-log2(d$y[d$y>0])*d$y[d$y>0]*dx)
    ent.old=ent.norm(s)
    ent.diff=ent.old-ent.new

    print.srcref(sprintf("%20s: diff.ent=%7.3f", nn, ent.diff))
  }
}


# Specific plots for showing over-dispersion:

#                Name    mean      sd      q1  median      q3   lower   upper
#          betabin.s1  0.0931  0.0266  0.0742  0.0901  0.1086  0.0502  0.1540
#          betabin.s2  0.0516  0.0229  0.0353  0.0476  0.0635  0.0184  0.1071
#          betabin.s3  0.0401  0.0185  0.0269  0.0374  0.0504  0.0123  0.0839
#          betabin.s4  0.0727  0.0138  0.0630  0.0716  0.0812  0.0490  0.1031

kk=0:60

png("overdisp_superspecies.png",height=1600,width=2000)
plot(kk,dbinom(kk,60,.4),col="red",type="l",xlab="counts",ylab="Prob.")
lines(kk,dbetabin(kk,60,0.4,0.07))
dev.off()

png("overdisp1.png",height=1600,width=2000)
plot(kk,dbinom(kk,60,.05),col="red",type="l",xlab="counts",ylab="Prob.")
lines(kk,dbetabin(kk,60,0.05,0.09))
dev.off()

png("overdisp2.png",height=1600,width=2000)
plot(kk,dbinom(kk,60,.03),col="red",type="l",xlab="counts",ylab="Prob.")
lines(kk,dbetabin(kk,60,0.03,0.05))
dev.off()

png("overdisp3.png",height=1600,width=2000)
plot(kk,dbinom(kk,60,.02),col="red",type="l",xlab="counts",ylab="Prob.")
lines(kk,dbetabin(kk,60,0.02,0.04))
dev.off()


png("overdisp_all.png",height=1600,width=2000)
par(mfrow=c(2,2),cex=1.5)
plot(kk,dbinom(kk,60,.05),col="red",type="l",xlab="counts",ylab="Prob.",
  main="Antarctothoa_tongima",lwd=3)
lines(kk,dbetabin(kk,60,0.05,0.09),lwd=3)
plot(kk,dbinom(kk,60,.03),col="red",type="l",xlab="counts",ylab="Prob.",
 main="Escharoides_excavata",lwd=3)
lines(kk,dbetabin(kk,60,0.03,0.05),lwd=3)
plot(kk,dbinom(kk,60,.02),col="red",type="l",xlab="counts",ylab="Prob.",
  main="Arachnopusia_unicornis",lwd=3)
lines(kk,dbetabin(kk,60,0.02,0.04),lwd=3)
plot(kk,dbinom(kk,60,.4),col="red",type="l",xlab="counts",ylab="Prob.",
  main="Superspecies",lwd=3)
lines(kk,dbetabin(kk,60,0.4,0.07),lwd=3)
dev.off()

