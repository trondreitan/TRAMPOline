source("read_data.R")

models=c("f_sf","lambda_f_sf","f_sf_ou","f_sf_o18bin","f_sf_cabin")
n.mod=length(models)

ra.mean=array(NA,c(n.mod,Nf,N.sp))
for(i in 1:n.mod)
{
  load(sprintf("mean_ra_%s.Rdata",models[i]))
  ra.mean[i,,]=mean.rel.abundance	
}

ra.lower=array(NA,c(n.mod,Nf,N.sp))
for(i in 1:n.mod)
{
  load(sprintf("lower_ra_%s.Rdata",models[i]))
  ra.lower[i,,]=lower.rel.abundance	
}

ra.upper=array(NA,c(n.mod,Nf,N.sp))
for(i in 1:n.mod)
{
  load(sprintf("upper_ra_%s.Rdata",models[i]))
  ra.upper[i,,]=upper.rel.abundance	
}

sp=names(y)[22:25]


cols=c("black","purple","red","green","blue")

png("rel_abundance_compare2.png",height=1200,width=1600)
par(cex=3,cex.lab=1.9,cex.sub=1.9,cex.axis=1.9,mar=c(10,5,1,1))
par(mfrow=c(2,2))
plot(K, ra.mean[1,,1], type="b", ylim=c(min(ra.mean[,,1]),max(ra.mean[,,1])),
  xlab="time", ylab=sprintf("Relative abundance: %s",sp[1]),lwd=4,log="y")
lines(K, ra.mean[2,,1],lty=2,col=cols[2],lwd=3)
lines(K, ra.mean[3,,1],lty=3,col=cols[3],lwd=3)
lines(K, ra.mean[4,,1],lty=4,col=cols[4],lwd=3)
lines(K, ra.mean[5,,1],lty=5,col=cols[5],lwd=3)

plot(K, ra.mean[1,,2], type="b", ylim=c(min(ra.mean[,,2]),max(ra.mean[,,2])),
  xlab="time", ylab=sprintf("Relative abundance: %s",sp[2]),lwd=4,log="y")
lines(K, ra.mean[2,,2],lty=2,col=cols[2],lwd=3)
lines(K, ra.mean[3,,2],lty=3,col=cols[3],lwd=3)
lines(K, ra.mean[4,,2],lty=4,col=cols[4],lwd=3)
lines(K, ra.mean[5,,2],lty=5,col=cols[5],lwd=3)

plot(K, ra.mean[1,,3], type="b", ylim=c(min(ra.mean[,,3]),max(ra.mean[,,3])),
  xlab="time", ylab=sprintf("Relative abundance: %s",sp[3]),lwd=4,log="y")
lines(K, ra.mean[2,,3],lty=2,col=cols[2],lwd=3)
lines(K, ra.mean[3,,3],lty=3,col=cols[3],lwd=3)
lines(K, ra.mean[4,,3],lty=4,col=cols[4],lwd=3)
lines(K, ra.mean[5,,3],lty=5,col=cols[5],lwd=3)

plot(K, ra.mean[1,,4], type="b", ylim=c(min(ra.mean[,,4]),max(ra.mean[,,4])),
  xlab="time", ylab=sprintf("Relative abundance: %s","Superspecies"),lwd=4,log="y")
lines(K, ra.mean[2,,4],lty=2,col=cols[2],lwd=3)
lines(K, ra.mean[3,,4],lty=3,col=cols[3],lwd=3)
lines(K, ra.mean[4,,4],lty=4,col=cols[4],lwd=3)
lines(K, ra.mean[5,,4],lty=5,col=cols[5],lwd=3)
dev.off()




png("rel_abundance_compare_uncert.png",height=1200,width=1600)
par(cex=3,cex.lab=1.9,cex.sub=1.9,cex.axis=1.9,mar=c(10,5,1,1))
par(mfrow=c(2,2))
plot(K, ra.mean[1,,1], type="b", ylim=c(min(ra.lower[,,1]),max(ra.upper[,,1])),
  xlab="time", ylab=sprintf("Relative abundance: %s",sp[1]),lwd=4,log="y")
lines(K, ra.mean[2,,1],lty=1,col=cols[2],lwd=3)
#lines(K, ra.mean[3,,1],lty=1,col=cols[3],lwd=3)
#lines(K, ra.mean[4,,1],lty=1,col=cols[4],lwd=3)
#lines(K, ra.mean[5,,1],lty=1,col=cols[5],lwd=3)
lines(K, ra.lower[1,,1],lty=2,col=cols[1])
lines(K, ra.lower[2,,1],lty=2,col=cols[2])
#lines(K, ra.lower[3,,1],lty=2,col=cols[3])
#lines(K, ra.lower[4,,1],lty=2,col=cols[4])
#lines(K, ra.lower[5,,1],lty=2,col=cols[5])
lines(K, ra.upper[1,,1],lty=2,col=cols[1])
lines(K, ra.upper[2,,1],lty=2,col=cols[2])
#lines(K, ra.upper[3,,1],lty=2,col=cols[3])
#lines(K, ra.upper[4,,1],lty=2,col=cols[4])
#lines(K, ra.upper[5,,1],lty=2,col=cols[5])

plot(K, ra.mean[1,,2], type="b", ylim=c(min(ra.lower[,,2]),max(ra.upper[,,2])),
  xlab="time", ylab=sprintf("Relative abundance: %s",sp[2]),lwd=4,log="y")
lines(K, ra.mean[2,,2],lty=1,col=cols[2],lwd=3)
#lines(K, ra.mean[3,,2],lty=1,col=cols[3],lwd=3)
#lines(K, ra.mean[4,,2],lty=1,col=cols[4],lwd=3)
#lines(K, ra.mean[5,,2],lty=1,col=cols[5],lwd=3)
lines(K, ra.lower[1,,2],lty=2,col=cols[1])
lines(K, ra.lower[2,,2],lty=2,col=cols[2])
#lines(K, ra.lower[3,,2],lty=2,col=cols[3])
#lines(K, ra.lower[4,,2],lty=2,col=cols[4])
#lines(K, ra.lower[5,,2],lty=2,col=cols[5])
lines(K, ra.upper[1,,2],lty=2,col=cols[1])
lines(K, ra.upper[2,,2],lty=2,col=cols[2])
#lines(K, ra.upper[3,,2],lty=2,col=cols[3])
#lines(K, ra.upper[4,,2],lty=2,col=cols[4])
#lines(K, ra.upper[5,,2],lty=2,col=cols[5])

plot(K, ra.mean[1,,3], type="b", ylim=c(min(ra.lower[,,3]),max(ra.upper[,,3])),
  xlab="time", ylab=sprintf("Relative abundance: %s",sp[3]),lwd=4,log="y")
lines(K, ra.mean[2,,3],lty=1,col=cols[2],lwd=3)
#lines(K, ra.mean[3,,3],lty=1,col=cols[3],lwd=3)
#lines(K, ra.mean[4,,3],lty=1,col=cols[4],lwd=3)
#lines(K, ra.mean[5,,3],lty=1,col=cols[5],lwd=3)
lines(K, ra.lower[1,,3],lty=2,col=cols[1])
lines(K, ra.lower[2,,3],lty=2,col=cols[2])
#lines(K, ra.lower[3,,3],lty=2,col=cols[3])
#lines(K, ra.lower[4,,3],lty=2,col=cols[4])
#lines(K, ra.lower[5,,3],lty=2,col=cols[5])
lines(K, ra.upper[1,,3],lty=2,col=cols[1])
lines(K, ra.upper[2,,3],lty=2,col=cols[2])
#lines(K, ra.upper[3,,3],lty=2,col=cols[3])
#lines(K, ra.upper[4,,3],lty=2,col=cols[4])
#lines(K, ra.upper[5,,3],lty=2,col=cols[5])

plot(K, ra.mean[1,,4], type="b", ylim=c(min(ra.lower[,,4]),max(ra.upper[,,4])),
  xlab="time", ylab=sprintf("Relative abundance: %s","Superspecies"),lwd=4,log="y")
lines(K, ra.mean[2,,4],lty=1,col=cols[2],lwd=3)
#lines(K, ra.mean[3,,4],lty=1,col=cols[3],lwd=3)
#lines(K, ra.mean[4,,4],lty=1,col=cols[4],lwd=3)
#lines(K, ra.mean[5,,4],lty=1,col=cols[5],lwd=3)
lines(K, ra.lower[1,,4],lty=2,col=cols[1])
lines(K, ra.lower[2,,4],lty=2,col=cols[2])
#lines(K, ra.lower[3,,4],lty=2,col=cols[3])
#lines(K, ra.lower[4,,4],lty=2,col=cols[4])
#lines(K, ra.lower[5,,4],lty=2,col=cols[5])
lines(K, ra.upper[1,,4],lty=2,col=cols[1])
lines(K, ra.upper[2,,4],lty=2,col=cols[2])
#lines(K, ra.upper[3,,4],lty=2,col=cols[3])
#lines(K, ra.upper[4,,4],lty=2,col=cols[4])
#lines(K, ra.upper[5,,4],lty=2,col=cols[5])
dev.off()


