
#load("/Users/leehl/Desktop/Occu 20.8.2021/plotting_variables.Rdata")
load("plotting_variables.Rdata")

sp=c( "Antarctothoa tongima", "Escharoides excavata",
   "Arachnopusia unicornis" ,"Superspecies")

sp2=c( "A. tongima", "E. excavata",
   "A. unicornis" ,"Superspecies")

species=c(substitute(paste(italic('A. tongima'))),
          substitute(paste(italic('E. excavata'))),
          substitute(paste(italic('A. unicornis'))),
          substitute(paste('Superspecies')))


#####################
# Summary plots:
#####################


# Plot raw detection rate data:

png("raw_data.png",height=1200,width=1600)
par(cex=3) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
par(mfrow=c(2,2)) 
plot(K[y$form.nr], y[,names(y)==sp[1]]/y$Total,xlab="millions of years ago",ylab=sprintf("Raw ratio: %s",sp[1]),ylim=c(0,1))
lines(K, mean.raw[,1])
plot(K[y$form.nr], y[,names(y)==sp[2]]/y$Total,xlab="millions of years ago",ylab=sprintf("Raw ratio: %s",sp[2]),ylim=c(0,1))
lines(K, mean.raw[,2])
plot(K[y$form.nr], y[,names(y)==sp[3]]/y$Total,xlab="millions of years ago",ylab=sprintf("Raw ratio: %s",sp[3]),ylim=c(0,1))
lines(K, mean.raw[,3])
plot(K[y$form.nr], y[,names(y)==sp[4]]/y$Total,xlab="millions of years ago",ylab=sprintf("Raw ratio: %s",sp[4]),ylim=c(0,1))
lines(K, mean.raw[,4])
dev.off()




# Formation-wise adjustments to overall occupancy:
png("log_odds_occupancy_formation.png",height=1200,width=1600)
par(cex=3) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
plot(K, f.occu.mean, type="b", ylim=c(-5,7),
  xlab="millions of years ago", ylab="Log odds ratio, occupancy",lwd=5)
lines(K, f.occu.median,type="b",lty=2,lwd=2)
lines(K, f.occu.lower,type="b",lty=2,lwd=2)
lines(K, f.occu.upper,type="b",lty=2,lwd=2)
dev.off()

png("odds_occupancy_formation.png",height=1200,width=1600)
par(cex=3) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
plot(K, exp(f.occu.mean), type="b", ylim=c(0.01,500),
  xlab="millions of years ago", ylab="Odds ratio, occupancy",log="y",lwd=5)
lines(K, exp(f.occu.median),type="b",lty=2,lwd=2)
lines(K, exp(f.occu.lower),type="b",lty=2,lwd=2)
lines(K, exp(f.occu.upper),type="b",lty=2,lwd=2)
dev.off()

# Formation-wise adjustments to detection:
png("log_odds_detection_formation.png",height=1200,width=1600)
par(cex=3) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
plot(K, f.bin.mean, type="b", ylim=c(-2.4,2),
  xlab="millions of years ago", ylab="Log-ddds ratio, detection",lwd=5)
lines(K, f.bin.median,type="b",lty=2,lwd=2)
lines(K, f.bin.lower,type="b",lty=2,lwd=2)
lines(K, f.bin.upper,type="b",lty=2,lwd=2)
dev.off()

png("odds_detection_formation.png",height=1200,width=1600)
par(cex=3) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
plot(K, exp(f.bin.mean), type="b", ylim=c(0.1,10),
  xlab="millions of years ago", ylab="Odds ratio, detection",log="y",lwd=5)
lines(K, exp(f.bin.median),type="b",lty=2,lwd=2)
lines(K, exp(f.bin.lower),type="b",lty=2,lwd=2)
lines(K, exp(f.bin.upper),type="b",lty=2,lwd=2)
dev.off()



# Formation-species-wise adjustments to overall occupancy:
png("log_odds_occupancy_species_formation.png",height=1200,width=1600)
par(cex=3, mfrow=c(2,2)) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
for(s in 1:(N.sp-1))
{
 plot(K, sf.occu.mean[s,], type="b", ylim=c(min(sf.occu.lower[s,]),
                                            max(sf.occu.upper[s,])),
  xlab="millions of years ago", ylab="Log odds ratio, occupancy",lwd=5,main=sp[s])
 lines(K, sf.occu.median[s,],type="b",lty=2,lwd=2)
 lines(K, sf.occu.lower[s,],type="b",lty=2,lwd=2)
 lines(K, sf.occu.upper[s,],type="b",lty=2,lwd=2)
}
dev.off()

png("odds_occupancy_species_formation.png",height=1200,width=1600)
par(cex=3, mfrow=c(2,2)) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
for(s in 1:(N.sp-1))
{
 plot(K, exp(sf.occu.mean[s,]), type="b", ylim=c(min(exp(sf.occu.lower[s,])),
                                            max(exp(sf.occu.upper[s,]))),
  xlab="millions of years ago", ylab="Log odds ratio, occupancy",lwd=5,main=sp[s])
 lines(K, exp(sf.occu.median[s,]),type="b",lty=2,lwd=2)
 lines(K, exp(sf.occu.lower[s,]),type="b",lty=2,lwd=2)
 lines(K, exp(sf.occu.upper[s,]),type="b",lty=2,lwd=2)
}
dev.off()

# Formation-species-wise adjustments to overall occupancy:
png("log_odds_detection_species_formation.png",height=1200,width=1600)
par(cex=3, mfrow=c(2,2)) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
for(s in 1:N.sp)
{
 plot(K, sf.bin.mean[s,], type="b", ylim=c(min(sf.bin.lower[s,]),
                                            max(sf.bin.upper[s,])),
  xlab="millions of years ago", ylab="Log odds ratio, detection",lwd=5,main=sp[s])
 lines(K, sf.bin.median[s,],type="b",lty=2,lwd=2)
 lines(K, sf.bin.lower[s,],type="b",lty=2,lwd=2)
 lines(K, sf.bin.upper[s,],type="b",lty=2,lwd=2)
}
dev.off()

png("odds_detection_species_formation.png",height=1200,width=1600)
par(cex=3, mfrow=c(2,2)) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
for(s in 1:N.sp)
{
 plot(K, exp(sf.bin.mean[s,]), type="b", ylim=c(min(exp(sf.bin.lower[s,])),
                                            max(exp(sf.bin.upper[s,]))),
  xlab="millions of years ago", ylab="Log odds ratio, detection",lwd=5,main=sp[s])
 lines(K, exp(sf.bin.median[s,]),type="b",lty=2,lwd=2)
 lines(K, exp(sf.bin.lower[s,]),type="b",lty=2,lwd=2)
 lines(K, exp(sf.bin.upper[s,]),type="b",lty=2,lwd=2)
}
dev.off()






# Per formation per species occupancy

# All species:
png("occupancy_species_formation.png",height=1200,width=1600)
par(cex=3,cex.lab=2.9,cex.sub=1.9,cex.axis=1.9,cex.main=1.9,mar=c(10,5,6,1))
cols=c("red","green","blue", "black")
tstr=sprintf("Red=%s,Green=%s,Blue=%s",sp[1],sp[2],sp[3])
par(mfrow=c(1,1))
plot(K, occu.sf.mean[1,], type="b", ylim=c(0.4,1),
  xlab="millions of years ago", ylab="Occupancy",main=tstr,col=cols[1],lwd=4,log="y")
lines(K, occu.sf.mean[2,],col=cols[2],lwd=4)
lines(K, occu.sf.mean[3,], type="b",col=cols[3],lwd=4)
dev.off()

png("occupancy_species_formation_uncert.png",height=1200,width=1600)
par(cex=2,cex.lab=2.9,cex.sub=1.9,cex.axis=1.9,mar=c(10,5,1,1))
par(mfrow=c(2,2))
for(s in 1:(N.sp-1))
{
 plot(K, occu.sf.mean[s,], type="b", ylim=c(0,1),
   xlab="millions of years ago", ylab=sprintf("Occupancy: %s",species[s]),lwd=4)
 lines(K, occu.raw[,s],lty=2,lwd=2,type="b",pch=2)
 lines(K, occu.sf.lower[s,],lty=3,lwd=3,type="b",pch=4)
 lines(K, occu.sf.upper[s,],lty=3,lwd=3,type="b",pch=4)
}
dev.off()




# Per formation per species detection
png("detection_species_formation.png",height=1200,width=1600)
par(cex=3,cex.lab=2.9,cex.sub=1.9,cex.axis=1.9,cex.main=1.9,mar=c(10,5,6,1))
cols=c("red","green","blue", "black")
tstr=sprintf("red=%s,green=%s\nblue=%s,black=superspecies",sp[1],sp[2],sp[3])
par(mfrow=c(1,1))
plot(K, bin.sf.mean[1,], type="b", ylim=c(0.005,1),
  xlab="millions of years ago", ylab="Detection",main=tstr,col=cols[1],lwd=4,log="y")
lines(K, bin.sf.mean[2,],col=cols[2],lwd=4)
lines(K, bin.sf.mean[3,], type="b",col=cols[3],lwd=4)
lines(K, bin.sf.mean[4,], type="b",col=cols[4],lwd=4)
dev.off()

png("detection_species_formation_uncert.png",height=1200,width=1600)
par(cex=3,cex.lab=2.9,cex.sub=1.9,cex.axis=1.9,mar=c(10,5,1,1))
par(mfrow=c(2,2))
for(s in 1:N.sp)
{
 plot(K, bin.sf.mean[s,], type="b", ylim=c(0.0001,1),
   xlab="millions of years ago", ylab=sprintf("Detection: %s",species[s]),lwd=4,log="y")
 lines(K, bin.sf.median[s,],lty=2)
 lines(K, bin.sf.lower[s,],lty=2)
 lines(K, bin.sf.upper[s,],lty=2)
}
dev.off()



png("occupancy_detection_uncert.png",height=1800,width=3200)
par(cex=3,cex.lab=9.4,cex.sub=12.9,cex.main=7.4,cex.axis=5.4,mar=c(17,22,8,1),tcl=-2,font.main=3)
par(mfrow=c(2,4))
for(s in 1:(N.sp-1))
{
 ytext=""
 if(s==1)
   ytext="Occupancy prob.\n"
 plot(K, occu.sf.mean[s,], type="b", ylim=c(0,1),
   xlab="", ylab=ytext,axes=FALSE)
 title(species[s],line=5)
 axis(1, at=c(-2, -1.5, -1,-0.5), labels=NA,lwd=9)
 if(s==1)
   axis(2, at=seq(0,1,0.2), 
     labels=seq(0,1,0.2),lwd=9,padj=-1.5)
 if(s>1)
   axis(2, at=seq(0,1,0.2), 
     labels=NA,lwd=9,padj=-1.5)
 box(lwd=9)
 #lines(K, occu.raw[,s],lty=2,lwd=2,type="b",pch=2)
 lines(K, occu.sf.lower[s,],lty=4,lwd=3,type="b",pch=4)
 lines(K, occu.sf.upper[s,],lty=4,lwd=3,type="b",pch=4)
}
plot(0,0,axes=FALSE,xlab="",ylab="",col="white",pch=NA)
title(species[4],line=5)
for(s in 1:N.sp)
{
 ytext=""
 if(s==1)
   ytext="Detection prob.\n"
 plot(K, bin.sf.mean[s,], type="b", ylim=c(0.0001,1),
   xlab="", ylab=ytext,lwd=6,log="y",axes=FALSE)
 axis(1, at=c(-2, -1.5, -1,-0.5),
   labels=c(-2,-1.5,-1,-0.5),lwd=9,padj=1.5)
 if(s==1)
   axis(2, at=c(0.001, 0.01, 0.1,1), 
     labels=c("0.001",  "0.01", "0.1","1"),lwd=9,padj=-1.5)
 if(s>1)
   axis(2, at=c(0.001, 0.01, 0.1,1), 
     labels=NA,lwd=9,padj=-1.5)
 box(lwd=9)
 #title(species[s],line=5)
 #lines(K, bin.sf.median[s,],lty=3)
 lines(K, bin.sf.lower[s,],lty=4)
 lines(K, bin.sf.upper[s,],lty=4)
}

mtext(side=1, "millions of years ago", outer=T, cex=5.5, padj=-0.5)

dev.off()






R_plotdata = data.frame(x=rep(K,each=4), y=as.vector(t(mean.rel.abundance)), Species= rep(sp2,times=9))
plot1 = ggplot(R_plotdata, aes(x=x, y=y, fill=Species)) +
  geom_area() +
  xlab("Millions of years ago") +
  ylab("Relative species abundance, R")
 
plot2 = plot1 + coord_cartesian(ylim = c(0.6, 1))
plot3 = plot2 + scale_fill_manual(values = alpha(c("grey50","grey75","black","grey90"),0.2))
plot4 = plot3 + scale_x_discrete(labels=c("-2.0" = "2", "-1.5" = "1.5",   "-1" = "1", "-0.5" = "0.5"))
png("R_stackced_area.png", res=300, height=2200, width=1200)
plot4 + theme(legend.position="bottom", legend.direction = "vertical")
dev.off()



png("R.png",height=1200,width=1600)
par(cex=3,cex.lab=2.9,cex.sub=1.9,cex.axis=2.2,cex.main=1.9,mar=c(10,5,6,1))
cols=c("red","green","blue", "black")
tstr=sprintf("red=%s,green=%s\nblue=%s,black=superspecies",sp[1],sp[2],sp[3])
par(mfrow=c(1,1))
plot(K, mean.rel.abundance[,1], type="b", ylim=c(0.001,1),
  xlab="millions of years ago", ylab="mean R",main=tstr,col=cols[1],lwd=4,log="y")
lines(K, mean.rel.abundance[,2], type="b",col=cols[2],lwd=4)
lines(K, mean.rel.abundance[,3], type="b",col=cols[3],lwd=4)
lines(K, mean.rel.abundance[,4], type="b",col=cols[4],lwd=4)
dev.off()





png("R_uncert.png",height=2200,width=3200)
par(cex=3,cex.lab=6.4,cex.sub=12.9,cex.main=7.4,cex.axis=5.4,mar=c(18,32,10,1),tcl=-2,font.main=3)
par(mfrow=c(2,2))
plot(K, mean.rel.abundance[,1], type="l", ylim=c(0.0001,1),
  main=species[1],lwd=9,log="y",axes=F,xlab="",ylab="",cex=6)
points(K, mean.rel.abundance[,1], cex=6, col="white",lwd=26)
points(K, mean.rel.abundance[,1], cex=3, col="black",lwd=15)
axis(1, at=c(-2, -1.5, -1,-0.5), labels=NA,lwd=9)
axis(2, at=c(0.001, 0.01, 0.1, 1), labels=c("0.001", "0.01",  "0.1", "1"),lwd=9,padj=-1.5)
box(lwd=9)
#lines(K, median.rel.abundance[,1],lty=2)
lines(K, lower.rel.abundance[,1],lwd=9,col="grey")
lines(K, upper.rel.abundance[,1],lwd=9,col="grey")

plot(K, mean.rel.abundance[,2], type="b", ylim=c(0.0001,1),
  main=species[2],lwd=9,log="y",axes=F,xlab="",ylab="",cex=6)
points(K, mean.rel.abundance[,2], cex=6, col="white",lwd=26)
points(K, mean.rel.abundance[,2], cex=3, col="black",lwd=15)
axis(1, at=c(-2, -1.5, -1,-0.5), labels=NA,lwd=9)
axis(2, at=c(0.001, 0.01, 0.1, 1), labels=c("0.001", "0.01",  "0.1", "1"),lwd=9,padj=-1.5)
box(lwd=9)
#lines(K, median.rel.abundance[,2],lty=2)
lines(K, lower.rel.abundance[,2],lwd=9,col="grey")
lines(K, upper.rel.abundance[,2],lwd=9,col="grey")

plot(K, mean.rel.abundance[,3], type="b", ylim=c(0.0001,1),
  main=species[3],lwd=9,log="y",axes=F,xlab="",ylab="",cex=6)
points(K, mean.rel.abundance[,3], cex=6, col="white",lwd=26)
points(K, mean.rel.abundance[,3], cex=3, col="black",lwd=15)
axis(1, at=c(-2, -1.5, -1,-0.5), labels=c("2", "1.5", "1","0.5"),lwd=9,padj=1.5)
axis(2, at=c(0.001, 0.01, 0.1, 1), labels=c("0.001", "0.01",  "0.1", "1"),lwd=9,padj=-1.5)
box(lwd=9)
#lines(K, median.rel.abundance[,3],lty=2)
lines(K, lower.rel.abundance[,3],lwd=9,col="grey")
lines(K, upper.rel.abundance[,3],lwd=9,col="grey")

plot(K, mean.rel.abundance[,4], type="b", ylim=c(0.5,1),
  main=species[4],lwd=9,log="y",axes=F,xlab="",ylab="",cex=6)
points(K, mean.rel.abundance[,4], cex=6, col="white",lwd=26)
points(K, mean.rel.abundance[,4], cex=3, col="black",lwd=15)
axis(1, at=c(-2, -1.5, -1,-0.5), labels=c("2", "1.5", "1","0.5"),lwd=9,padj=1.5)
axis(2, at=c(0.5, 0.75, 1), labels=c("0.5", "0.75", "1"),lwd=9,padj=-1.5)
box(lwd=9)
#lines(K, median.rel.abundance[,4],lty=2)
lines(K, lower.rel.abundance[,4],lwd=9,col="grey")
lines(K, upper.rel.abundance[,4],lwd=9,col="grey")

mtext(side=1, "millions of years ago", outer=T, cex=7.5, padj=-0.5)
mtext(side=2, "Relative species abundance, R", outer=T, cex=7.5, padj=2)

dev.off()





png("rel_abundance_formation.png",height=1200,width=1600)
par(cex=3,cex.lab=2, cex.sub=2,cex.main=2,cex.axis=2)
par(mfrow=c(3,3))
for(i in 1:Nf)
   barplot(A[1:(N.sp-1),i],names.arg=species[1:(N.sp-1)],
     main=sprintf("Formation age %2.2f-%2.2f",-Ke[i], -Ks[i]))
dev.off()




png("Q.png",height=1200,width=1600)
par(cex=3,cex.lab=2.9,cex.sub=1.9,cex.axis=2.2,cex.main=1.9,mar=c(10,5,6,1))
cols=c("red","green","blue", "black")
tstr=sprintf("red=%s,green=%s\nblue=%s,black=superspecies",sp[1],sp[2],sp[3])
par(mfrow=c(1,1))
plot(K, mean.Q[,1], type="b", ylim=c(0,4),
  xlab="millions of years ago",   
  ylab="mean Q",main=tstr,col=cols[1],lwd=4)
lines(K, mean.Q[,2], type="b",col=cols[2],lwd=4)
lines(K, mean.Q[,3], type="b",col=cols[3],lwd=4)
lines(K, mean.Q[,4], type="b",col=cols[4],lwd=4)
abline(c(1,0),lty=2)
dev.off()





png("Q_uncert.png",height=2200,width=3200)
par(cex=3,cex.lab=6.4,cex.sub=12.9,cex.main=7.4,cex.axis=5.4,mar=c(18,32,10,1),tcl=-2,font.main=3)
par(mfrow=c(2,2))
plot(K, mean.Q[,1], type="l", ylim=c(0,6),
  main=species[1],lwd=9,axes=F,xlab="",ylab="",cex=6)
points(K, mean.Q[,1], cex=6, col="white",lwd=26)
points(K, mean.Q[,1], cex=3, col="black",lwd=15)
axis(1, at=c(-2, -1.5, -1,-0.5), labels=NA,lwd=9)
axis(2, at=c(0, 2, 4, 6), labels=c("0", "2",  "4", "6"),lwd=9,padj=-1.5)
box(lwd=9)
#lines(K, median.Q[,1],lty=2)
lines(K, lower.Q[,1],lwd=9,col="grey")
lines(K, upper.Q[,1],lwd=9,col="grey")
abline(c(1,0),lty=2,lwd=4)

plot(K, mean.Q[,2], type="b", ylim=c(0,6),
  main=species[2],lwd=9,axes=F,xlab="",ylab="",cex=6)
points(K, mean.Q[,2], cex=6, col="white",lwd=26)
points(K, mean.Q[,2], cex=3, col="black",lwd=15)
axis(1, at=c(-2, -1.5, -1,-0.5), labels=NA,lwd=9)
axis(2, at=c(0, 2, 4, 6), labels=c("0", "2",  "4", "6"),lwd=9,padj=-1.5)
box(lwd=9)
#lines(K, median.Q[,2],lty=2)
lines(K, lower.Q[,2],lwd=9,col="grey")
lines(K, upper.Q[,2],lwd=9,col="grey")
abline(c(1,0),lty=2,lwd=4)


plot(K, mean.Q[,3], type="b", ylim=c(0,6),
  main=species[3],lwd=9,axes=F,xlab="",ylab="",cex=6)
points(K, mean.Q[,3], cex=6, col="white",lwd=26)
points(K, mean.Q[,3], cex=3, col="black",lwd=15)
axis(1, at=c(-2, -1.5, -1,-0.5), labels=c("2", "1.5", "1","0.5"),lwd=9,padj=1.5)
axis(2, at=c(0, 2, 4, 6), labels=c("0", "2",  "4", "6"),lwd=9,padj=-1.5)
box(lwd=9)
#lines(K, median.Q[,3],lty=2)
lines(K, lower.Q[,3],lwd=9,col="grey")
lines(K, upper.Q[,3],lwd=9,col="grey")
abline(c(1,0),lty=2,lwd=4)


plot(K, mean.Q[,4], type="b", ylim=c(0,6),
  main=species[4],lwd=9,axes=F,xlab="",ylab="",cex=6)
points(K, mean.Q[,4], cex=6, col="white",lwd=26)
points(K, mean.Q[,4], cex=3, col="black",lwd=15)
axis(1, at=c(-2, -1.5, -1,-0.5), labels=c("2", "1.5", "1","0.5"),lwd=9,padj=1.5)
axis(2, at=c(0, 2, 4, 6), labels=c("0", "2",  "4", "6"),lwd=9,padj=-1.5)
box(lwd=9)
#lines(K, median.Q[,4],lty=2)
lines(K, lower.Q[,4],lwd=9,col="grey")
lines(K, upper.Q[,4],lwd=9,col="grey")
abline(c(1,0),lty=2,lwd=4)


mtext(side=1, "millions of years ago", outer=T, cex=7.5, padj=-0.4)
mtext(side=2, "Relative population density, Q", outer=T, cex=7.5, padj=2)

dev.off()





# Plot raw detection rate data with ocpcuancy*detection

png("raw_data_with_occu_times_det.png",height=1200,width=1600)
par(cex=4,cex.lab=1.5,cex.sub=1.5,cex.axis=1.5,
  mar=c(5,8,4,2))
par(mfrow=c(2,2))
for(s in 1:N.sp)
{
 plot(K[y$form.nr], y[,names(y)==sp[s]]/y$Total,xlab="millions of years ago",ylab=sprintf("Ratio: %s",sp[s]),ylim=c(0,max(y[,names(y)==sp[s]]/y$Total)),lwd=2)
 lines(K, mean.raw[,s],lwd=2)
 lines(K, detection.sf.mean[s,], lty=3, lwd=3)
 lines(K, detection.sf.lower[s,], lty=2, lwd=1)
 lines(K, detection.sf.upper[s,], lty=2, lwd=1)
}
dev.off()


####LH's try 
postscript("try.eps",width=3.56,height=4)
y_values=t(mean.rel.abundance)
x_values=K

par(mfrow=c(2,1), mar=c(1,2,1,1),cex.main=0.9,cex.axis=0.7)
plot(x_values, y_values[4,],  xlim=c(-2.3, 0), ylim=c(0.6,1), type="n", xaxs="i", yaxs="i", axes=F, xlab="Millions of years ago", ylab="R", main=substitute(paste("Combined")))
polygon(x=c(rev(x_values), x_values), y=c(rep(0,9), y_values[4,]), col="black", border=NA)
polygon(x=c(rev(x_values), x_values), y=c(rev(y_values[4,]+y_values[3,]), (y_values[4,])), col="grey30", border=NA)
polygon(x=c(rev(x_values), x_values), y=c(rev(y_values[4,]+y_values[3,]+y_values[2,]), (y_values[4,]+y_values[3,])), col="grey50", border=NA)
polygon(x=c(rev(x_values), x_values), y=c(rep(1,9), (y_values[4,]+y_values[3,]+y_values[2,])), col="gray70", border=NA)
axis(1, label=NA, at= c(-2,-1.5, -1,-0.5), tck=-0.035)
axis(2, label=c(0.6,0.8,1), at=c(0.6,0.8,1),line=-0.73)


plot(x_values, y_values[4,],  xlim=c(-2.3, 0), ylim=c(0.6,1),  axes=F, 
     xlab="", ylab="",type="n", xaxs="i", yaxs="i")
legend(-2,0.9,legend=c(expression(italic("A. tonigma")), expression(italic("A. unicornis")), expression(italic("E. excavata")), "Superspecies"), 
       col=c("gray80", "gray50", "gray70", "black"), bty="n", pch=c(15,15,15,15), pt.cex=2)

dev.off()
