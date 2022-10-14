
# Install the package. Check if this goes well before proceeding!
# Might require installation of Rcpp first.
#install.packages("https://github.com/trondreitan/layeranalyzer/raw/master/layeranalyzer_0.1.0.tar.gz",type="source",verbose=T)

# If this went well, load the package:
library(layeranalyzer)


# Number of results. Should ideally not exceed the number 
# of cores on the machine, or the analysis will go slower
N=4


# Fetch N results in parallel:
res=list() # Define the results list
i=1:N
# Define a function that performs the analysis 
# and returns the result. We perform a causal connection
# analysis on bivalve extinctionr ate to brachiopod 
# extinction rate, using a 2-layered process for each 
# time series. 
doanalysis=function(j) 
{
  library(layeranalyzer)  
  
  # 2-layered search
  bi.lext.struct2=layer.series.structure(bi.lext, 
    numlayers=2,prior=lrate.pr)
  br.lext.struct2=layer.series.structure(br.lext, 
    numlayers=2,prior=lrate.pr)

  layer.analyzer(bi.lext.struct2,br.lext.struct2,
	causal=cbind(c(1,1,2,1)), num.MCMC=400, spacing=10,
     use.stationary.stdev=TRUE, silent.mode=FALSE, num.temp=1,
     mcmc=TRUE)
}

# Perform the N parallel analyses:
library(parallel)
cl <- makeCluster(getOption("cl.cores", N))
res=clusterApply(cl, 1:N, doanalysis)


# Put the MCMC resulst into a list:
library(coda)
mcmcres=mcmc.list(res[[1]]$mcmc)
for(i in 2:N)
  mcmcres[[i]]=res[[i]]$mcmc

# Plot all the parameters (PS: may simply scroll 
# through most of them fast):
num.par=dim(mcmcres[[1]])[2]
plot(mcmcres)

# Perform the Gelman statistics. The closer it is to 1.0, 
# the better:
gelman.diag(mcmcres)
#                                 Point est. Upper C.I.
#mu_bivalve_lext                        1.00       1.00
#dt_bivalve_lext_1                      1.15       1.34
#sd_bivalve_lext_1                      1.12       1.31
#dt_bivalve_lext_2                      1.13       1.31
#sd_bivalve_lext_2                      1.08       1.23
#mu_brach_lext                          1.00       1.00
#dt_brach_lext_1                        1.02       1.07
#sd_brach_lext_1                        1.01       1.03
#dt_brach_lext_2                        1.12       1.16
#sd_brach_lext_2                        1.00       1.00
#beta_bivalve_lext,1_brach_lext,1       1.01       1.03
#Multivariate psrf
#1.09

# Not too shabby. I would call this convergence.
# Thus no parallel tempering needed in this case.
# Let's still take a look at the "worst" parameter,
# parameter number 2, "dt_bivalve_lext_1":

j=2
pname=res[[1]]$parameter.names[j]

# Plot the MCMC chains for this parameter:
mcmcres.worst=array(0,c(dim(mcmcres[[1]])[1],N))
for(i in 1:N)
  mcmcres.worst[,i]=res[[i]]$mcmc[,j]
matplot(mcmcres.worst, type="l", ylab=pname)
# Also looks like convergence. However, there's 
# quite a bit of auto-correlation, apparently

# collate into one array, and fetch the auto-correlation:
a=c(mcmcres.worst)
rho=acf(a)$acf[2]
rho

# Quite a lot of auto-correlation

# Number of independent samples:
length(a)*(1-rho)/(1+rho)
# 84

# Number of samples for each independent sample:
(rho+1)/(1-rho)
# 19

# That's at least some independent sample. 
# Still, It would probably be better 
# to use a much larger spacing than the default.
# spacing=200 should work well, but the analysis will take 
# 20 times longer than for spacing=10 (the default).







