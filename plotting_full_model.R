
load("plotting_variables.Rdata")


#####################
# Summary plots:
#####################


# Plot raw detection rate data:

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
plot(K, f.occu.mean, type="b", ylim=c(-5,7),
  xlab="time", ylab="Log odds ratio, occupancy",lwd=5)
lines(K, f.occu.median,type="b",lty=2,lwd=2)
lines(K, f.occu.lower,type="b",lty=2,lwd=2)
lines(K, f.occu.upper,type="b",lty=2,lwd=2)
dev.off()

png("odds_occupancy_formation.png",height=1200,width=1600)
par(cex=3) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
plot(K, exp(f.occu.mean), type="b", ylim=c(0.01,500),
  xlab="time", ylab="Odds ratio, occupancy",log="y",lwd=5)
lines(K, exp(f.occu.median),type="b",lty=2,lwd=2)
lines(K, exp(f.occu.lower),type="b",lty=2,lwd=2)
lines(K, exp(f.occu.upper),type="b",lty=2,lwd=2)
dev.off()

# Formation-wise adjustments to detection:
png("log_odds_detection_formation.png",height=1200,width=1600)
par(cex=3) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
plot(K, f.bin.mean, type="b", ylim=c(-2.4,2),
  xlab="time", ylab="Log-ddds ratio, detection",lwd=5)
lines(K, f.bin.median,type="b",lty=2,lwd=2)
lines(K, f.bin.lower,type="b",lty=2,lwd=2)
lines(K, f.bin.upper,type="b",lty=2,lwd=2)
dev.off()

png("odds_detection_formation.png",height=1200,width=1600)
par(cex=3) #,cex.lab=1.3,cex.sub=1.3,cex.axis=1.3)
plot(K, exp(f.bin.mean), type="b", ylim=c(0.1,10),
  xlab="time", ylab="Odds ratio, detection",log="y",lwd=5)
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
  xlab="time", ylab="Log odds ratio, occupancy",lwd=5,main=sp[s])
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
  xlab="time", ylab="Log odds ratio, occupancy",lwd=5,main=sp[s])
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
  xlab="time", ylab="Log odds ratio, detection",lwd=5,main=sp[s])
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
  xlab="time", ylab="Log odds ratio, detection",lwd=5,main=sp[s])
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
  xlab="time", ylab="Occupancy",main=tstr,col=cols[1],lwd=4,log="y")
lines(K, occu.sf.mean[2,],col=cols[2],lwd=4)
lines(K, occu.sf.mean[3,], type="b",col=cols[3],lwd=4)
dev.off()

png("occupancy_species_formation_uncert.png",height=1200,width=1600)
par(cex=2,cex.lab=2.9,cex.sub=1.9,cex.axis=1.9,mar=c(10,5,1,1))
par(mfrow=c(2,2))
for(s in 1:(N.sp-1))
{
 plot(K, occu.sf.mean[s,], type="b", ylim=c(0,1),
   xlab="time", ylab=sprintf("Occupancy: %s",species[s]),lwd=4)
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
  xlab="time", ylab="Detection",main=tstr,col=cols[1],lwd=4,log="y")
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
   xlab="time", ylab=sprintf("Detection: %s",species[s]),lwd=4,log="y")
 lines(K, bin.sf.median[s,],lty=2)
 lines(K, bin.sf.lower[s,],lty=2)
 lines(K, bin.sf.upper[s,],lty=2)
}
dev.off()





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



png("rel_abundance_formation.png",height=1200,width=1600)
par(cex=3,cex.lab=2, cex.sub=2,cex.main=2,cex.axis=2)
par(mfrow=c(3,3))
for(i in 1:Nf)
   barplot(A[1:(N.sp-1),i],names.arg=species[1:(N.sp-1)],
     main=sprintf("Formation age %2.2f-%2.2f",-Ke[i], -Ks[i]))
dev.off()






# Plot raw detection rate data with ocpcuancy*detection

png("raw_data_with_occu_times_det.png",height=1200,width=1600)
par(cex=4,cex.lab=1.5,cex.sub=1.5,cex.axis=1.5,
  mar=c(5,8,4,2))
par(mfrow=c(2,2))
for(s in 1:N.sp)
{
 plot(K[y$form.nr], y[,names(y)==sp[s]]/y$Total,xlab="Time",ylab=sprintf("Ratio: %s",sp[s]),ylim=c(0,max(y[,names(y)==sp[s]]/y$Total)),lwd=2)
 lines(K, mean.raw[,s],lwd=2)
 lines(K, detection.sf.mean[s,], lty=3, lwd=3)
 lines(K, detection.sf.lower[s,], lty=2, lwd=1)
 lines(K, detection.sf.upper[s,], lty=2, lwd=1)
}
dev.off()

