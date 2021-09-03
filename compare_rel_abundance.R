source("read_data.R")

models=c("f_sf","lambda_f_sf")
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

# Re-done

F=FALSE
cols=c("black","darkgrey")

png("rel_abundance_compare_uncert.png",height=1200,width=1600)
par(mar=c(5,10,5,5), oma=c(7,4,0,0))
par(cex=3,cex.lab=3,cex.sub=3,cex.axis=3,cex.main=3)
par(mfrow=c(2,2))
plot(K, ra.mean[1,,1], type="b",    
   ylim=c(0.001,1), xlab="", ylab="",
   main=substitute(paste(italic('A. tongima'))),   
   axes=F, pch=20, cex=2.5, col=cols[1],
   lwd=6,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=NA)
axis(2, at=c(0.001, 0.01, 0.1, 1), labels=c(0.001,0.01, 0.1, 1))
box()
lines(K, ra.mean[2,,1],lty=1,col=cols[2],lwd=6)
lines(K, ra.lower[1,,1],lty=2,col=cols[1])
lines(K, ra.lower[2,,1],lty=2,col=cols[2])
lines(K, ra.upper[1,,1],lty=2,col=cols[1])
lines(K, ra.upper[2,,1],lty=2,col=cols[2])

plot(K, ra.mean[1,,2], type="b",   
   ylim=c(1e-14,1),xlab="", ylab="",
   main=substitute(paste(italic('E. excavata'))),
   axes=F, pch=20, cex=2.5, col=cols[1],
   lwd=6,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=NA)
axis(2, at=c(1e-14, 1e-9, 0.0001, 1), labels=c(1e-14, 1e-9, 0.0001, 1))
box()
lines(K, ra.mean[2,,2],lty=1,col=cols[2],lwd=6)
lines(K, ra.lower[1,,2],lty=2,col=cols[1])
lines(K, ra.lower[2,,2],lty=2,col=cols[2])
lines(K, ra.upper[1,,2],lty=2,col=cols[1])
lines(K, ra.upper[2,,2],lty=2,col=cols[2])

plot(K, ra.mean[1,,3], type="b",    
   ylim=c(1e-10,1),xlab="", ylab="",
   main=substitute(paste(italic('A. unicornis'))), 
   axes=F, pch=20, cex=2.5, col=cols[1],
   lwd=6,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=seq(2, 0, -0.5),padj=1)
axis(2, at=c(1e-10, 1e-6, 0.01, 1), labels=c(1e-10, 1e-6,  0.01, 1))
box()
lines(K, ra.mean[2,,3],lty=1,col=cols[2],lwd=6)
lines(K, ra.lower[1,,3],lty=2,col=cols[1])
lines(K, ra.lower[2,,3],lty=2,col=cols[2])
lines(K, ra.upper[1,,3],lty=2,col=cols[1])
lines(K, ra.upper[2,,3],lty=2,col=cols[2])

plot(K, ra.mean[1,,4], type="b",    
   ylim=c(0.6,1),xlab="", ylab="",
   main="Superspecies",
   axes=F, pch=20, cex=2.5, col=cols[1],
   lwd=6,log="y")
axis(1, at=seq(-2, 0, 0.5), labels=seq(2, 0, -0.5),padj=1)
axis(2, at=seq(0.6,1.0,0.1), labels=seq(0.6,1.0, 0.1))
box()
lines(K, ra.mean[2,,4],lty=1,col=cols[2],lwd=6)
lines(K, ra.lower[1,,4],lty=2,col=cols[1])
lines(K, ra.lower[2,,4],lty=2,col=cols[2])
lines(K, ra.upper[1,,4],lty=2,col=cols[1])
lines(K, ra.upper[2,,4],lty=2,col=cols[2])
mtext(side=1, "millions of years ago", outer=T, cex=3.8, padj=1)
mtext(side=2, "Relative species abundance, R\n", outer=T, cex=3.8, padj=1)
dev.off()



