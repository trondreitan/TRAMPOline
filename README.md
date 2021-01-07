This is a set of scripts for analyzing temporal multispecies occupancy models targeted towards estimating relative abundance estimates. Relative abundance is the abundance of a species compared to the other species in question for different time intervals (formations). For our case, we singled out 3 focus bryozoan species, and let the colonies of other species fall into a "superspecies" category. We made the scripts so that it would be easy to switch between models, using the variable called “model” to do this. These scripts fall into several categories.
1.	Model specification scripts. These have file names of the type "model_(modelname).R", where (modelname) is the name of the model. Each model file contains two functions for initializing parameters (either with or without parameter estimates from a simpler model), a hyper-parameter specification, a log-prior specification, and a log-likelihood-specification. Here is a list of out models in ascending model complexity order up to the full model (the more complex models are not nested inside each other):

a) global – Treats all data points as coming from the same formation, i.e. assumes a global occupancy and detection probability for each species. Note that the occupancy probability for the superspecies is set to 1). Also has an over-dispersion parameter for each species. Thus, the parameter set is 3xS-1, where S stands for the number of focus species plus the superspecies (in our case, S=4). 

b) foccu – In addition to the global parameters, it also has formation-dependent random effects for the occupancy probabilities. Thus, for each formation to it allows for having some species-independent deviances in the occupancy probability. It also has a parameter for the standard deviation of these random effects. Thus, it has 3S+F parameters (top parameters and random effects), where F is the number of formations.  

c) f – Has formation-dependent (but species-independent) random effects for both detection and occupancy probabilities, as well as standard deviations for both these random effects, so 3S+1+2F parameters. The formation-dependent random effects are there to carry away pure formation-dependent deviations from the relative abundance estimation.

d) f_sf – Has formation-dependent (but species-independent) plus formation- and species-dependent random effects. This is also called the full model, as it has all the necessary components for estimating relative abundance. 

e) f_sf_o18bin/f_sf_o18occu/f_sf_o18sep/f_sf_o18both – – In addition to the elements of the full  model, this model has ∂18O regression terms, either for detection probabilities (bin), occupancy probabilities (occu), both detection and occupancy probabilities separately (sep), or common random factors for detection and occupancy probabilities (both). 

f) f_sf_cabin/f_sf_caoccu/f_sf_casep/f_sf_caboth – In addition to the elements of the full  model, this model has Mg/Ca regression terms, either for detection probabilities (bin), occupancy probabilities (occu), both detection and occupancy probabilities separately (sep), or common random factors for detection and occupancy probabilities (both). 

g) f_sf_modelselect – Performs model selection analysis using MCMC sampling of indicator variables for the various regression ∂18O and Mg/Ca terms (see the two previous model sets).

h) f_sf_ou – An extended version of the full model where the random effects are not independent, but rather auto-correlated. We impose an Ornstein-Uhlenbeck (OU) process for describing how the random effects vary through time. The extra parameter introduced is characteristic time (the time for the correlation to drop to 1/e). Both the pure formation-dependent and the formation- and species-dependent random effects are given this treatment, for both occupancy and detection probability variables, with one characteristic time for each type of random effects (thus 4 in total). 

i) lambda_f_sf – This is not an extension of the full model, but rather an alternative to it which puts increased emphasis on relative abundance. The abundance estimates rather than the detection probabilities are used as a basis for an additive model with species-dependent global parameters, formation-dependent random effects and formations- and species-dependent random effects. The detection probabilities are then derived from these abundance estimates, for use in the likelihood. Since this model has a slightly different structure, it has its own relative abundance script (real_abundance_lambda.R), see below.

2.	Run scripts. These have file names of the type "run_(modelname).R" (see list of models above). These are simple files, specifying the model to be examined, the run number and the number of MCMC-samples, before calling the MCMC routine and then storing the result. Since the script only stores one result file, while MCMC analysis scripts rely on multiple result files, you should copy the run file multiple time and change the run number in each copy. 

3.	Support scripts. These are scripts providing useful functions and code for the run and model scripts and are generally made so that the model used can be changed to a different model easily. These are:

a) 3stage.R/3stage_cpu.R – These two scripts run the MCMC sampling using the LaplacesDemon package. This is done in three stages. The first runs "UESS" method, and the parameter set sample with the best likelihood is then chosen as the start for the next MCMC stage, which uses the "AMWG" method. The best parameter set according to the likelihood from that MCMC sampling is then used as the start for the third stage, using the "CHARM" method. it is this latter MCMC sampling which is then saved in a file named (modelname).Rdata". The script calls the other support scripts listed below. The "3stage_cpu.R" script allows for setting the number of CPUs (or cores) used by LaplacesDemon and is called by the run scripts of the more complex models.

b) read_data.R - Reads the data from the data file, in our case named "allsamples_with_counts_and_metainfo.csv", and creates the data variables the other scripts use.

c) init_rerun.R - The init functions should be defined in each "model_(modelname).R" script, but the functions in the init_reruns.R script call these multiple times and chooses the one with the best combined prior and likelihood, in order to start the first MCMC run from a good place.

d) make_laplace_wrapper.R - Creates the structures that the LaplacesDemon package expects, from the data and prior/likelihood-specifications.

e) find_best_par.R - Runs through the samples returned by LaplacesDemon and finds the parameter set with the highest likelihood.

4.	Analysis scripts. These are also made so that one model can be replaced by another easily, using the “model” variable.

a) lookat.R - Reads multiple output files for a given model (specified at the start of the script), lists the model log-likelihood estimates, makes an overall model likelihood from all runs, creates parameter mean, median and standard deviation files and gives parameter estimate summaries. The parameter mean, median and standard deviation files are important, as more complex models rely on such files from simpler models. (The "foccu" model relies on the "global" model results, the "f" model relies on "foccu" results and more complex models rely on the "f" results). These files have been provided for our dataset, so you don't have to rerun the simplest model over again.

b) lookat_parameters.R - Same as lookat.R but with more code focusing on parameter estimates.

c) lookat_modelselect.R - Special "lookat" file for the model selection case.

d) rel_abundance2.R - Calculates relative abundance right from the MCMC output files.

e) rel_abundance_start.R - Also calculates relative abundance, but stores the estimates in an output file called "plotting_variables" rather than plot directly.

f) rel_abundance_lambda.R - Calculates, saves and plots relative abundance estimates for the "lambda_f_sf" model.

g) compare_rel_abundance.R - Plots relative abundance estimates for multiple models. (One needs to use the rel_abundance_start.R and rel_abundance_lambda.R scripts 
for making the relative abundance output files used here).

h) plotting_full_model.R - Reads the plotting_variables.RData file from "rel_abundance_start.R" and plots from that.

5.	Data files: allsamples_with_counts_and_metainfo.csv. Bryozoan site-wise data for 3 focus species and 1 superspeces from the Wanganui basin. Contains ∂18O- or Mg/Ca data in addition, for regression models.

6.	Parameter estimate files: par_mean_global.RData, par_mean_foccu.RData, par_mean_f.RData, par_median_global.RData, par_median_foccu.RData, par_median_f.RData, par_sd_global.RData, par_sd_foccu.RData, par_sd_f.RData. Contains mean, median, and standard deviation from the MCMC samples of the simplest models ("global", "foccu" and "f").

