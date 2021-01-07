This is a set of scripts for analyzing temporal multispecies occupancy models targeted towards estimating relative abundance estimates. Relative abundance is the abundance of a species compared to the other species in question for different time intervals (formations). For our case, we singled out 3 focus bryozoan species, and let the colonies of other species fall into a "superspecies" category. We made the scripts so that it would be easy to switch between models, using the variable called “model” to do this. These scripts fall into several categories.

1.	Model specification scripts. These have file names of the type "model_(modelname).R", where (modelname) is the name of the model. Our full model is called "f_sf" ("f" for formation-dependent random effects and "sf" species- and formation-dependent random effects). There are also simpler models ("global", "foccu" (formation-dependent random effects for occupancy probabilities) and "f" (formation-dependent random effects for occupancy and detection probabilities)). There are also more advanced models: ∂18O- or Mg/Ca-regression ("f_sf_o18" and "f_sf_ca" variants respectively), either only on occupancy or detection probability, both separately or with the same effect for both, and an auto-regressive model named "f_sf_ou". There is also a larger model where model selection between the various regression elements is allowed, called "f_sf_modelselect". Lastly there is a model with the same structure as the full model, but with focus on local abundance rather than detection probability, called "lambda_f_sf". Each model file contains two functions for initializing parameters (either with or without parameter estimates from a simpler model), a hyper-parameter specification, a log-prior specification, and a log-likelihood-specification.

2.	Run scripts. These have file names of the type "run_(modelname).R". These are simple files, specifying the model to be examined, the run number and the number of MCMC-samples, before calling the MCMC routine and then storing the result. Since the script only stores one result file, while MCMC analysis scripts rely on multiple result files, you should copy the run file multiple time and change the run number in each copy.

3.	Support scripts. These are scripts providing useful functions and code for the run and model scripts and are generally made so that the model used can be changed to a different model easily. These are: 

      a)	3stage.R/3stage_cpu.R – These two scripts run the MCMC sampling using the LaplacesDemon package. This is done in three stages. The first runs "UESS" method, and the parameter set sample with the best likelihood is then chosen as the start for the next MCMC stage, which uses the "AMWG" method. The best parameter set according to the likelihood from that MCMC sampling is then used as the start for the third stage, using the "CHARM" method. it is this latter MCMC sampling which is then saved in a file named (modelname).Rdata". The script calls the other support scripts listed below. The "3stage_cpu.R" script allows for setting the number of CPUs (or cores) used by LaplacesDemon and is called by the run scripts of the more complex models. 

      b)	read_data.R - Reads the data from the data file, in our case named "allsamples_with_counts_and_metainfo.csv", and creates the data variables the other scripts use. 

      c)	init_rerun.R - The init functions should be defined in each "model_(modelname).R" script, but the functions in the init_reruns.R script call these multiple times and chooses the one with the best combined prior and likelihood, in order to start the first MCMC run from a good place. 

      d)	make_laplace_wrapper.R - Creates the structures that the LaplacesDemon package expects, from the data and prior/likelihood-specifications. 

      e)	find_best_par.R - Runs through the samples returned by LaplacesDemon and finds the parameter set with the highest likelihood.

4.	Analysis files. These are also made so that one model can be replaced by another easily, using the “model” variable.

      a)	lookat.R - Reads multiple output files for a given model (specified at the start of the script), lists the model log-likelihood estimates, makes an overall model likelihood from all runs, creates parameter mean, median and standard deviation files and gives parameter estimate summaries. The parameter mean, median and standard deviation files are important, as more complex models rely on such files from simpler models. (The "foccu" model relies on the "global" model results, the "f" model relies on "foccu" results and more complex models rely on the "f" results). These files have been provided for our dataset, so you don't have to rerun the simplest model over again. 
      b)	lookat_parameters.R - Same as lookat.R but with more code focusing on parameter estimates. 
      c)	lookat_modelselect.R - Special "lookat" file for the model selection case.  
      d)	rel_abundance2.R - Calculates relative abundance right from the MCMC output files. 
      e)	rel_abundance_start.R - Also calculates relative abundance, but stores the estimates in an output file called "plotting_variables" rather than plot directly. 
      f)	rel_abundance_lambda.R - Calculates, saves and plots relative abundance estimates for the "lambda_f_sf" model. 
      g)	compare_rel_abundance.R - Plots relative abundance estimates for multiple models. (One needs to use the rel_abundance_start.R and rel_abundance_lambda.R scripts for making the relative abundance output files used here). 
      h)	plotting_full_model.R - Reads the plotting_variables.RData file from "rel_abundance_start.R" and plots from that.

5.	Data file: allsamples_with_counts_and_metainfo.csv. Bryozoan site-wise data for 3 focus species and 1 superspeces from the Wanganui basin. Contains ∂18O- or Mg/Ca data in addition, for regression models.

6.	Parameter estimate files: par_mean_global.RData, par_mean_foccu.RData, par_mean_f.RData, par_median_global.RData, par_median_foccu.RData, par_median_f.RData, par_sd_global.RData, par_sd_foccu.RData, par_sd_f.RData. Contains mean, median, and standard deviation from the MCMC samples of the simplest models ("global", "foccu" and "f").

