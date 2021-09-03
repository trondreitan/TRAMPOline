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

sp[1]="Antarctothoa tongima"
sp[2]="Escharoides excavata"
sp[3]="Arachnopusia unicornis"

cols=c("black","purple","red","green","blue")

png("R_compare2.png",height=1200,width=1600)
par(cex=3,cex.main=2.4,cex.lab=2.4,cex.sub=1.9,cex.axis=1.9,mar=c(8,5,7,7))
par(mfrow=c(2,2))
plot(K, ra.mean[1,,1], type="b", ylim=c(min(ra.mean[,,1]),max(ra.mean[,,1])),ylab="",
  axes=F,xlab="", main=sprintf("%s\n",sp[1]),lwd=4,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=seq(2, 0, -0.5))
axis(2, at=c(0.001, 0.01, 0.1, 1), labels=c(0.001,0.01, 0.1, 1))
box()
#lines(K, ra.mean[2,,1],lty=2,col=cols[2],lwd=3)
lines(K, ra.mean[3,,1],lty=3,col=cols[3],lwd=3)
lines(K, ra.mean[4,,1],lty=4,col=cols[4],lwd=3)
lines(K, ra.mean[5,,1],lty=5,col=cols[5],lwd=3)

plot(K, ra.mean[1,,2], type="b", ylim=c(min(ra.mean[,,2]),max(ra.mean[,,2])),ylab="",
  axes=F,xlab="", main=sprintf("%s\n",sp[2]),lwd=4,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=seq(2, 0, -0.5))
axis(2, at=c(0.001, 0.01, 0.1, 1), labels=c(0.001,0.01, 0.1, 1))
box()
#lines(K, ra.mean[2,,2],lty=2,col=cols[2],lwd=3)
lines(K, ra.mean[3,,2],lty=3,col=cols[3],lwd=3)
lines(K, ra.mean[4,,2],lty=4,col=cols[4],lwd=3)
lines(K, ra.mean[5,,2],lty=5,col=cols[5],lwd=3)

plot(K, ra.mean[1,,3], type="b", ylim=c(min(ra.mean[,,3]),max(ra.mean[,,3])),ylab="",
  axes=F,xlab="", main=sprintf("%s\n",sp[3]),lwd=4,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=seq(2, 0, -0.5))
axis(2, at=c(0.001, 0.01, 0.1, 1), labels=c(0.001,0.01, 0.1, 1))
box()
#lines(K, ra.mean[2,,3],lty=2,col=cols[2],lwd=3)
lines(K, ra.mean[3,,3],lty=3,col=cols[3],lwd=3)
lines(K, ra.mean[4,,3],lty=4,col=cols[4],lwd=3)
lines(K, ra.mean[5,,3],lty=5,col=cols[5],lwd=3)

plot(K, ra.mean[1,,4], type="b", ylim=c(min(ra.mean[,,4]),max(ra.mean[,,4])),ylab="",
  axes=F,xlab="", main=sprintf("%s","Superspecies\n"),lwd=4,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=seq(2, 0, -0.5))
axis(2, at=c(0.001, 0.01, 0.1, 1), labels=c(0.001,0.01, 0.1, 1))
box()
#lines(K, ra.mean[2,,4],lty=2,col=cols[2],lwd=3)
lines(K, ra.mean[3,,4],lty=3,col=cols[3],lwd=3)
lines(K, ra.mean[4,,4],lty=4,col=cols[4],lwd=3)
lines(K, ra.mean[5,,4],lty=5,col=cols[5],lwd=3)
mtext(side=1, "millions of years ago", outer=T, cex=2.8, padj=1)
mtext(side=2, "R\n", outer=T, cex=3.8, padj=1)
dev.off()





png("rel_abundance_compare_uncert.png",height=1200,width=1600)
par(cex=3,cex.lab=2.8,cex.sub=2.8,cex.axis=2.8,cex.main=2.8, mar=c(10,8,1,1))
par(mfrow=c(2,2))
plot(K, ra.mean[1,,1], type="b", ylim=c(min(ra.lower[,,1]),max(ra.upper[,,1])),
   xlab="time", ylab=sprintf("Relative abundance\n%s",sp[1]),lwd=3,log="y")
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
  xlab="time", ylab=sprintf("Relative abundance\n%s",sp[2]),lwd=4,log="y")
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
  xlab="time", ylab=sprintf("Relative abundance\n%s",sp[3]),lwd=4,log="y")
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
  xlab="time", ylab=sprintf("Relative abundance\n%s","Superspecies"),lwd=4,log="y")
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




# Re-done

cols=c("black","darkgrey")

png("R_compare_uncert.png",height=1200,width=1600)
par(mar=c(5,10,5,5), oma=c(7,4,0,0))
par(cex=3,cex.lab=2.8,cex.sub=2.8,cex.axis=2.8,cex.main=2.8)
par(mfrow=c(2,2))
plot(K, ra.mean[1,,1], type="b",    
   ylim=c(0.001,1), xlab="", ylab="",
   main=substitute(paste(italic('A. tongima'))),   
   axes=F, pch=20, cex=2.5, col=cols[1],
   lwd=3,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=seq(2, 0, -0.5))
axis(2, at=c(0.001, 0.01, 0.1, 1), labels=c(0.001,0.01, 0.1, 1))
box()
lines(K, ra.mean[2,,1],lty=1,col=cols[2],lwd=3)
lines(K, ra.lower[1,,1],lty=2,col=cols[1])
lines(K, ra.lower[2,,1],lty=2,col=cols[2])
lines(K, ra.upper[1,,1],lty=2,col=cols[1])
lines(K, ra.upper[2,,1],lty=2,col=cols[2])

plot(K, ra.mean[1,,2], type="b",   
   ylim=c(1e-14,1),xlab="", ylab="",
   main=substitute(paste(italic('E. excavata'))),
   axes=F, pch=20, cex=2.5, col=cols[1],
   lwd=3,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=seq(2, 0, -0.5))
axis(2, at=c(1e-14, 1e-9, 0.0001, 1), labels=c(1e-14, 1e-9, 0.0001, 1))
box()
lines(K, ra.mean[2,,2],lty=1,col=cols[2],lwd=3)
lines(K, ra.lower[1,,2],lty=2,col=cols[1])
lines(K, ra.lower[2,,2],lty=2,col=cols[2])
lines(K, ra.upper[1,,2],lty=2,col=cols[1])
lines(K, ra.upper[2,,2],lty=2,col=cols[2])

plot(K, ra.mean[1,,3], type="b",    
   ylim=c(1e-10,1),xlab="", ylab="",
   main=substitute(paste(italic('A. unicornis'))), 
   axes=F, pch=20, cex=2.5, col=cols[1],
   lwd=3,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=seq(2, 0, -0.5))
axis(2, at=c(1e-10, 1e-6, 0.01, 1), labels=c(1e-10, 1e-6,  0.01, 1))
box()
lines(K, ra.mean[2,,3],lty=1,col=cols[2],lwd=3)
lines(K, ra.lower[1,,3],lty=2,col=cols[1])
lines(K, ra.lower[2,,3],lty=2,col=cols[2])
lines(K, ra.upper[1,,3],lty=2,col=cols[1])
lines(K, ra.upper[2,,3],lty=2,col=cols[2])

plot(K, ra.mean[1,,4], type="b",    
   ylim=c(0.5,1),xlab="", ylab="",
   main="Super Species",
   axes=F, pch=20, cex=2.5, col=cols[1],
   lwd=3,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=seq(2, 0, -0.5))
axis(2, at=seq(0.5,1.0,0.1), labels=seq(0.5,1.0, 0.1))
box()
lines(K, ra.mean[2,,4],lty=1,col=cols[2],lwd=3)
lines(K, ra.lower[1,,4],lty=2,col=cols[1])
lines(K, ra.lower[2,,4],lty=2,col=cols[2])
lines(K, ra.upper[1,,4],lty=2,col=cols[1])
lines(K, ra.upper[2,,4],lty=2,col=cols[2])
mtext(side=1, "millions of years ago", outer=T, cex=3.8, padj=1)
mtext(side=2, "R\n", outer=T, cex=3.8, padj=1)
dev.off()




# 4 models


cols=c("black","purple","red","green","blue")

png("R_compare.png",height=1200,width=1600)
par(mar=c(5,10,5,5), oma=c(7,4,0,0))
par(cex=3,cex.lab=2.8,cex.sub=2.8,cex.axis=2.8,cex.main=2.8)
par(mfrow=c(2,2))
plot(K, ra.mean[1,,1], type="b",    
   ylim=c(0.01,1), xlab="", ylab="",
   main=substitute(paste(italic('A. tongima'))),   
   axes=F, pch=20, cex=2.5, col=cols[1],
   lwd=3,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=seq(2, 0, -0.5),padj=1)
axis(2, at=c(0.001, 0.01, 0.1, 1), labels=c(0.001,0.01, 0.1, 1))
box()
#lines(K, ra.mean[2,,1],lty=1,col=cols[2],lwd=3)
lines(K, ra.mean[3,,1],lty=3,col=cols[3],lwd=3)
lines(K, ra.mean[4,,1],lty=4,col=cols[4],lwd=3)
lines(K, ra.mean[5,,1],lty=5,col=cols[5],lwd=3)


plot(K, ra.mean[1,,2], type="b",   
   ylim=c(0.0001,1),xlab="", ylab="",
   main=substitute(paste(italic('E. excavata'))),
   axes=F, pch=20, cex=2.5, col=cols[1],
   lwd=3,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=seq(2, 0, -0.5),padj=1)
axis(2, at=c(1e-4, 1e-2,  1), labels=c("0.0001", "0.01", "1"))
box()
#lines(K, ra.mean[2,,2],lty=2,col=cols[2],lwd=3)
lines(K, ra.mean[3,,2],lty=3,col=cols[3],lwd=3)
lines(K, ra.mean[4,,2],lty=4,col=cols[4],lwd=3)
lines(K, ra.mean[5,,2],lty=5,col=cols[5],lwd=3)


plot(K, ra.mean[1,,3], type="b",    
   ylim=c(0.001,1),xlab="", ylab="",
   main=substitute(paste(italic('A. unicornis'))), 
   axes=F, pch=20, cex=2.5, col=cols[1],
   lwd=3,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=seq(2, 0, -0.5),padj=1)
axis(2, at=c(1e-3, 1e-2, 0.1, 1), labels=c(1e-3, 1e-2,  0.1, 1))
box()
#lines(K, ra.mean[2,,3],lty=2,col=cols[2],lwd=3)
lines(K, ra.mean[3,,3],lty=3,col=cols[3],lwd=3)
lines(K, ra.mean[4,,3],lty=4,col=cols[4],lwd=3)
lines(K, ra.mean[5,,3],lty=5,col=cols[5],lwd=3)

plot(K, ra.mean[1,,4], type="b",    
   ylim=c(0.5,1),xlab="", ylab="",
   main="Super Species",
   axes=F, pch=20, cex=2.5, col=cols[1],
   lwd=3,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=seq(2, 0, -0.5),padj=1)
axis(2, at=seq(0.5,1.0,0.1), labels=seq(0.5,1.0, 0.1))
box()
#lines(K, ra.mean[2,,4],lty=2,col=cols[2],lwd=3)
lines(K, ra.mean[3,,4],lty=3,col=cols[3],lwd=3)
lines(K, ra.mean[4,,4],lty=4,col=cols[4],lwd=3)
lines(K, ra.mean[5,,4],lty=5,col=cols[5],lwd=3)


mtext(side=1, "millions of years ago", outer=T, cex=3.8, padj=1)
mtext(side=2, "Relative species abundance, R\n", outer=T, cex=3.8, padj=1)
dev.off()






