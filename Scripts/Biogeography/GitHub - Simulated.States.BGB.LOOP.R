# Start by running through the DEC/DECj analyses for the empirical dataset
# we need to get the 'res' object to estimate our null model q matrix from

setwd("/YOUR_DIRECTORY/BioGeoBEARS")
library(optimx)   
library(FD)       
library(parallel)
library(BioGeoBEARS)
library(snow)
library(ggplot2)
library(ape)
library(rexpokit)
library(cladoRcpp)
library(progress)
source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R") # (needed now that traits model added; source FIRST!)
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R") # added traits model
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
# slight speedup hopefully

#extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
wd <- getwd()

#############################################################################################
#############################################################################################
## Step 1: we are running the empirical data through BGB to get the empirical estimates
## of cladogenetic dispersal types and frequencies. We're going to use DEC and DEC+j models
## and then based on our model (DEC+j) and tree (user input), we'll create a stochastic set
## of possible alternatives (Step 2), to incorporate some method of uncertainty in the ASR
#############################################################################################
#############################################################################################

# Data Inputs
#######################################################
## drop in your tree
trfn = np(paste(addslash(wd), "Data.TREES/BGB.Meliphagides.tre", sep="")) #import the tree file
moref(trfn) # Look at the raw Newick file
tr = read.tree(trfn)
#plot(tr)

## if using empirical data
geogfn = np(paste(addslash(wd), "Data.GEOG.files/BGB.Meliphagides.geog.txt", sep="")) #empirical data
moref(geogfn) # Look at the raw geography text file:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn) # Look at your geographic range data:
max_range_size = 5 # Set the maximum number of areas any species may occupy

# Run DEC
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE    # get ancestral states from optim run

BioGeoBEARS_run_object$max_range_size = max_range_size # Input the maximum range size
BioGeoBEARS_run_object$num_cores_to_use=8 # Multicore processing if desired
BioGeoBEARS_run_object$force_sparse=FALSE
BioGeoBEARS_run_object$geogfn = geogfn # Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$trfn = trfn # Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object) # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object # Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object$BioGeoBEARS_model_object # This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table # This table contains the parameters of the model 
check_BioGeoBEARS_run(BioGeoBEARS_run_object) # Run this to check inputs. Read the error messages if you get them!

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "Sim_DEC_M0_unconstrained_v1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}
# Run DEC+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$speedup=TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use=1
BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$calc_ancprobs=TRUE    # get ancestral states from optim run

BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001
# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Sim_DEC+J_M0_unconstrained_v1.Rdata"
runslow = TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}

# Set up empty tables to hold the statistical results
restable = NULL
teststable = NULL

# Statistics -- DEC vs. DEC+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)
numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats
# DEC, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
# confer the same likelihood on the data. See: Brian O'Meara's webpage:
# http://www.brianomeara.info/tutorials/aic
# ...for an intro to LRT, AIC, and AICc
rbind(res2, res1)
tmp_tests = conditional_format_table(stats)
restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)


#############################################################################################
#############################################################################################
## Step 2: Biogeographic stochastic maps are going to include some (here 50 bsm) estimate of
## uncertainty from the ancestral state reconstructions of the biogeographic history
## we are creating at each node. We're going to pull out the range of possible scenarios in
## Step 3, so we can compare our empirical and simulated data
#############################################################################################
#############################################################################################



########################################
## Run BSMs on Empirical data
########################################
clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0

BSM_inputs_fn = "BSM_inputs_file.Rdata"
runInputsSlow = TRUE
if (runInputsSlow)
{
  stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
  save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
} else {
  # Loads to "stochastic_mapping_inputs_list"
  load(BSM_inputs_fn)
} # END if (runInputsSlow)

# Check inputs (doesn't work the same on unconstr)
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))

runBSMslow = TRUE
if (runBSMslow == TRUE)
{
  # Saves to: RES_clado_events_tables.Rdata
  # Saves to: RES_ana_events_tables.Rdata
  BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=100, nummaps_goal=50, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
  
  RES_clado_events_tables = BSM_output$RES_clado_events_tables
  RES_ana_events_tables = BSM_output$RES_ana_events_tables
} else {
  # Load previously saved...
  
  # Loads to: RES_clado_events_tables
  load(file="RES_clado_events_tables.Rdata")
  # Loads to: RES_ana_events_tables
  load(file="RES_ana_events_tables.Rdata")
  BSM_output = NULL
  BSM_output$RES_clado_events_tables = RES_clado_events_tables
  BSM_output$RES_ana_events_tables = RES_ana_events_tables
} # END if (runBSMslow == TRUE)

# Extract/Summarize BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables
areanames = names(tipranges@df)
actual_names = areanames

# Get the dmat and times (if any)
dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
dmat_times

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables

# Simulate the source areas
BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
ana_events_tables = BSMs_w_sourceAreas$ana_events_tables

# Count all anagenetic and cladogenetic events
counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)

summary_counts_BSMs = counts_list$summary_counts_BSMs
print(conditional_format_table(summary_counts_BSMs))

# Histogram of event counts
hist_event_counts(counts_list, pdffn=paste0("DEC+J", "_histograms_of_event_counts.pdf"))


#############################################################################################
#############################################################################################
## Step 3: From BSM results object, we're going to focus on the clado_events_tables, and
## summarizing the frequencies of different dispersal types. The four types of cladogenetic
## dispersal types (founder, vicariance, sympatry, subset symp) are pulled out
## alongside the total # of events, according to a time frame and time interval of interest
### and stored in the 'empirical' object.
#############################################################################################
#############################################################################################

########################################################################################
## Let's compare dispersal frequencies of different types within the same time period ##
########################################################################################

# pull out all the Vicariance (v) dispersal events from a single BSM
#vicariance.events<-subset(clado_events_tables[[1]], clado_events_tables[[1]]$clado_event_type == "vicariance (v)")
#...or just loop over all the BSMs and pull out all the vicariance events (object 'vicariance')
vicariance <- NULL
for (i in 1:length(clado_events_tables)) {
  vic<-subset(clado_events_tables[[i]], clado_events_tables[[i]]$clado_event_type == "vicariance (v)")
  vicariance<-rbind(vicariance, vic)
}

# pull out all the Founder (j) dispersal events from a single BSM
#founder<-subset(clado_events_tables[[1]], clado_events_tables[[1]]$clado_event_type == "founder (f)")
#...or just loop over all the BSMs and pull out all the founder events (object 'founder')
founder <- NULL
for (i in 1:length(clado_events_tables)) {
  found<-subset(clado_events_tables[[i]], clado_events_tables[[i]]$clado_event_type == "founder (j)")
  founder<-rbind(founder, found)
}

# bind together the 'vicariance' and 'founder' tables
vj<-rbind(vicariance, founder)

# pull out the Sympatric (y) dispersal events
sym.y<-NULL
for (i in 1:length(clado_events_tables)) {
  symy<-subset(clado_events_tables[[i]], clado_events_tables[[1]]$clado_event_type == "sympatric (y)")
  sym.y<-rbind(sym.y, symy)
}
# pull out the Subset (s) dispersal events
sub.s<-NULL
for (i in 1:length(clado_events_tables)) {
  sub<-subset(clado_events_tables[[i]], clado_events_tables[[1]]$clado_event_type == "subset (s)")
  sub.s<-rbind(sub.s, sub)
}

# bind together both sympatric events (Sympatric 'y' and Subset 's')
ys<-rbind(sym.y, sub.s)

# bind together the vicariance (v) and sympatric events (Sympatric 'y' and Subset 's')
vys<-rbind(ys, vicariance)

# do the same steps for ALL events
# start by pulling ALL of the nonterminal dispersal events
total <- NULL
for (i in 1:length(clado_events_tables)) {
  non.terminal<-subset(clado_events_tables[[i]], !clado_events_tables[[i]]$clado_event_type == "")
  total<-rbind(total, non.terminal)
}

# now let's loop through both sets and establish a ratio of dispersal frequencies within a time period
# make sure to define your "window size" with the 'time.max' object
empirical <- NULL
nums <- seq(0, (max(node.height(tr))), .1) #set your sliding window by (start, end, distance of window movement)
for (i in nums) { #use this if you want the total age
  #for (i in 0:14) { #use this if you want just a subset of the age
  time.min <- i-0.5 #sets the window start date (here, 500,000 years)
  time.max <- i+0.5 #and sets the window close date (in total, the window width is 1 million years)
  
  vic.found <- subset(vj, time_bp>time.min & time_bp<time.max)
  num.vj <- nrow(vic.found)
  
  sym.sub <- subset(ys, time_bp>time.min & time_bp<time.max)
  num.ys <- nrow(sym.sub)
  
  vic.sym <- subset(vys, time_bp>time.min & time_bp<time.max)
  num.vys <- nrow(vic.sym)
  
  all <- subset(total, time_bp>time.min & time_bp<time.max)
  num.all <- nrow(all)
  
  ratio.vj <- num.vj/num.all
  ratio.ys <- num.ys/num.all
  ratio.vys <- num.vys/num.all
  
  ratio.vjys <- (num.vj+num.ys)/num.all
  
  x <- matrix(nrow=1, ncol=8)
  x[,1] <- num.all
  x[,2] <- num.vj
  x[,3] <- ratio.vj
  x[,4] <- num.ys
  x[,5] <- ratio.ys
  x[,6] <- ratio.vjys
  x[,7] <- num.vys
  x[,8] <- ratio.vys
  
  
  row.names(x)<-paste(time.min,'to', time.max)
  empirical <- rbind(empirical, x)
}
colnames(empirical)<-c('num.all', 'num.vj', 'ratio.vj', 'num.ys', 
                       'ratio.ys', 'ratio.vjys', 'num.vys', 'ratio.vys')
empirical<-as.data.frame(empirical)
empirical[is.na(empirical)] <- 0
write.table(empirical, file="BGB.Skinks.MCC.REFINED.Dispersal.Frequencies.txt", sep="\t", row.names=T, quote=F)

### the "empirical" table above is important, because we'll call from it later!



#############################################################################################
#############################################################################################
## Step 4: Now we need to (repeatedly) simulate data based on our tree and model. First is to
## build our Q-matrix based off the 'res' object from the empirical data. 
## Then, a giant loop which runs the simulated data under the DEC & DEC+j models (pseudo Step 1), 
## follows this with a 50 BSMs (pseudo Step 2), and then summarizes the results of the BSMs
## (pseudo Step 3). Remember, this is all on the simulated data! The loop will iterate through 
## 100 times, after which we'll have to extract our 95% confidence intervals (Step 5). This
## loop can take a WHILE, so check progress by opening the "Progress.txt" document.
#############################################################################################
#############################################################################################

# Now we can move on to the simulations
########################################
# *below each task are the original functions saved for posterity

# first step: simulating the data

#######################################################
# Get the Q matrix and cladogenesis mode from the
# BioGeoBEARS_results_object "res"
#######################################################
#get_Qmat_COOmat_from_BioGeoBEARS_run_object <- function(BioGeoBEARS_run_object, numstates=NULL, include_null_range=TRUE)
q.mat <- get_Qmat_COOmat_from_res(res, numstates=5, include_null_range=F)

#check to make sure you're using the correct inputs for your qmat
res$inputs$geogfn
res$inputs$trfn

#let's make a loop for all the simulated things
frame <- NULL
full.table.vj <- NULL; full.table.ys <- NULL; full.table.vys <- NULL
bsm.vicariance <- NULL
bsm.founder <- NULL
bsm.sympatry <- NULL
bsm.subset <- NULL
all.events <- NULL
progress <- NULL
#pb <- progress_bar$new(total = 100)

for (t in 1:100) {
  #######################################################
  # simulate_biogeog_history
  #######################################################
  #simulate_biogeog_history <- function(phy, Qmat, COO_probs_columnar, index_Qmat_0based_of_starting_state)
  sim.history <- simulate_biogeog_history(tr, q.mat$Qmat, q.mat$COO_weights_columnar, 1)
  
  ##############################################################
  # Convert simulated Qmat 0-based indexes to a tipranges object
  ##############################################################
  #simulated_indexes_to_tipranges_object <- function(simulated_states_by_node, areas_list, states_list, trfn)
  areanames = c("T", "U", "S", "F", "A")
  list_0based_states = rcpp_areas_list_to_states_list(areas=areanames, maxareas=max_range_size, include_null_range=F)
  sim.tips.object <- simulated_indexes_to_tipranges_object(sim.history, areanames, list_0based_states, trfn) 
  
  ##############################################################
  # Convert simulated Qmat 0-based indexes to a tipranges file
  ##############################################################
  #simulated_indexes_to_tipranges_file <- function(simulated_states_by_node, areas_list, states_list, trfn, out_geogfn="lagrange_area_data_file.data")
  sim.tips.file <- simulated_indexes_to_tipranges_file(sim.history, areanames, list_0based_states, trfn, out_geogfn="SIMULATED_area_data_file.data")
  
  # Make a counter for the loops so we know how far we are
  #pb$tick()
  #Sys.sleep(1 / 100)
  progress <- print(t)
  write(progress, file="Sim.Progress.txt", append=T) #remember to trash the old progress file before starting!
  
  
  # next step: analyzing our simulated data
  
  # Data Inputs
  #######################################################
  
  ## if using simulated data
  tipranges = sim.tips.object #simulated data
  geogfn = geogfn = np(paste(addslash(wd), "SIMULATED_area_data_file.data", sep="")) #empirical data
  moref(geogfn)
  max_range_size = 5 # Set the maximum number of areas any species may occupy
  
  
  # Run DEC
  #######################################################
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
  BioGeoBEARS_run_object$speedup=TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE
  BioGeoBEARS_run_object$calc_ancprobs=TRUE    # get ancestral states from optim run
  
  BioGeoBEARS_run_object$max_range_size = max_range_size # Input the maximum range size
  BioGeoBEARS_run_object$num_cores_to_use=8 # Multicore processing if desired
  BioGeoBEARS_run_object$force_sparse=FALSE
  BioGeoBEARS_run_object$geogfn = geogfn # Give BioGeoBEARS the location of the geography text file
  BioGeoBEARS_run_object$trfn = trfn # Give BioGeoBEARS the location of the phylogeny Newick file
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object) # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE
  BioGeoBEARS_run_object # Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
  BioGeoBEARS_run_object$BioGeoBEARS_model_object # This contains the model object
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table # This table contains the parameters of the model 
  check_BioGeoBEARS_run(BioGeoBEARS_run_object) # Run this to check inputs. Read the error messages if you get them!
  
  # For a slow analysis, run once, then set runslow=FALSE to just 
  # load the saved result.
  runslow = TRUE
  resfn = "SIMULATED_DEC_M0_unconstrained_v1.Rdata"
  if (runslow)
  {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    resDEC = res
  } else {
    # Loads to "res"
    load(resfn)
    resDEC = res
  }
  # Run DEC+J
  #######################################################
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = trfn
  BioGeoBEARS_run_object$geogfn = geogfn
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$speedup=TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE
  BioGeoBEARS_run_object$num_cores_to_use=1
  BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
  BioGeoBEARS_run_object$calc_ancprobs=TRUE    # get ancestral states from optim run
  
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE
  # Set up DEC+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDEC$outputs@params_table["d","est"]
  estart = resDEC$outputs@params_table["e","est"]
  jstart = 0.0001
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  # Add j as a free parameter
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  resfn = "SIMULATED_DEC+J_M0_unconstrained_v1.Rdata"
  runslow = TRUE
  if (runslow)
  {
    #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
    
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    
    resDECj = res
  } else {
    # Loads to "res"
    load(resfn)
    resDECj = res
  }
  
  # Set up empty tables to hold the statistical results
  restable = NULL
  teststable = NULL
  
  # Statistics -- DEC vs. DEC+J
  #######################################################
  # We have to extract the log-likelihood differently, depending on the 
  # version of optim/optimx
  LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
  LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)
  numparams1 = 3
  numparams2 = 2
  stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
  stats
  # DEC, null model for Likelihood Ratio Test (LRT)
  res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  # DEC+J, alternative model for Likelihood Ratio Test (LRT)
  res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  # The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
  # confer the same likelihood on the data. See: Brian O'Meara's webpage:
  # http://www.brianomeara.info/tutorials/aic
  # ...for an intro to LRT, AIC, and AICc
  rbind(res2, res1)
  tmp_tests = conditional_format_table(stats)
  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)
  
  # Stochastic mapping 
  #######################################################
  clado_events_tables = NULL
  ana_events_tables = NULL
  lnum = 0
  
  BSM_inputs_fn = "BSM_inputs_file.Rdata"
  runInputsSlow = TRUE
  if (runInputsSlow)
  {
    stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
    save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
  } else {
    # Loads to "stochastic_mapping_inputs_list"
    load(BSM_inputs_fn)
  } # END if (runInputsSlow)
  
  # Check inputs (doesn't work the same on unconstr)
  names(stochastic_mapping_inputs_list)
  stochastic_mapping_inputs_list$phy2
  stochastic_mapping_inputs_list$COO_weights_columnar
  stochastic_mapping_inputs_list$unconstr
  set.seed(seed=as.numeric(Sys.time()))
  
  runBSMslow = TRUE
  if (runBSMslow == TRUE)
  {
    # Saves to: RES_clado_events_tables.Rdata
    # Saves to: RES_ana_events_tables.Rdata
    BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=100, nummaps_goal=50, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
    
    RES_clado_events_tables = BSM_output$RES_clado_events_tables
    RES_ana_events_tables = BSM_output$RES_ana_events_tables
  } else {
    # Load previously saved...
    
    # Loads to: RES_clado_events_tables
    load(file="RES_clado_events_tables.Rdata")
    # Loads to: RES_ana_events_tables
    load(file="RES_ana_events_tables.Rdata")
    BSM_output = NULL
    BSM_output$RES_clado_events_tables = RES_clado_events_tables
    BSM_output$RES_ana_events_tables = RES_ana_events_tables
  } # END if (runBSMslow == TRUE)
  
  # Extract BSM output
  clado_events_tables = BSM_output$RES_clado_events_tables
  ana_events_tables = BSM_output$RES_ana_events_tables
  
  #######################################################
  # Summarize stochastic map tables
  #######################################################
  areanames = names(tipranges@df)
  actual_names = areanames
  actual_names
  
  # Get the dmat and times (if any)
  dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)

  # Extract BSM output
  clado_events_tables = BSM_output$RES_clado_events_tables
  ana_events_tables = BSM_output$RES_ana_events_tables
  
  # Simulate the source areas
  BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
  clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
  ana_events_tables = BSMs_w_sourceAreas$ana_events_tables
  
  # Count all anagenetic and cladogenetic events
  counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)
  
  summary_counts_BSMs = counts_list$summary_counts_BSMs
  print(conditional_format_table(summary_counts_BSMs))
  
  
  ########################################################################################
  ## Let's compare dispersal frequencies of different types within the same time period ##
  ########################################################################################
  
  # pull out all the Vicariance (v) dispersal events from a single BSM
  #vicariance.events<-subset(clado_events_tables[[1]], clado_events_tables[[1]]$clado_event_type == "vicariance (v)")
  #...or just loop over all the BSMs and pull out all the vicariance events (object 'vicariance')
  vicariance <- NULL
  for (i in 1:length(clado_events_tables)) {
    vic<-subset(clado_events_tables[[i]], clado_events_tables[[i]]$clado_event_type == "vicariance (v)")
    vicariance<-rbind(vicariance, vic)
  }

  # pull out all the Founder (j) dispersal events from a single BSM
  #founder<-subset(clado_events_tables[[1]], clado_events_tables[[1]]$clado_event_type == "founder (f)")
  #...or just loop over all the BSMs and pull out all the founder events (object 'founder')
  founder <- NULL
  for (i in 1:length(clado_events_tables)) {
    found<-subset(clado_events_tables[[i]], clado_events_tables[[i]]$clado_event_type == "founder (j)")
    founder<-rbind(founder, found)
  }

  # bind together the 'vicariance' and 'founder' tables
  vj<-rbind(vicariance, founder)
  
  # pull out the Sympatric (y) dispersal events
  sym.y<-NULL
  for (i in 1:length(clado_events_tables)) {
    symy<-subset(clado_events_tables[[i]], clado_events_tables[[i]]$clado_event_type == "sympatry (y)")
    sym.y<-rbind(sym.y, symy)
  }

  # pull out the Subset (s) dispersal events
  sub.s<-NULL
  for (i in 1:length(clado_events_tables)) {
    sub<-subset(clado_events_tables[[i]], clado_events_tables[[i]]$clado_event_type == "subset (s)")
    sub.s<-rbind(sub.s, sub)
  }

  # bind together both sympatric events (Sympatric 'y' and Subset 's')
  ys<-rbind(sym.y, sub.s)
  
  # bind together the sympatric events (y & s) with the vicariance (v), ignoring the founder events
  vys<-rbind(ys, vicariance)
  
  # do the same steps for ALL events
  # start by pulling ALL of the nonterminal dispersal events
  total <- NULL
  for (i in 1:length(clado_events_tables)) {
    non.terminal<-subset(clado_events_tables[[i]], !clado_events_tables[[i]]$clado_event_type == "")
    total<-rbind(total, non.terminal)
  }
  
  # now let's loop through both sets and establish a ratio of dispersal frequencies within a time period
  # make sure to define your "window size" with the 'time.max' object
  frame <- NULL
  extra.butt <- NULL
  sim.nums <- seq(0, (max(node.height(tr))), .1) #set your sliding window by (start, end, distance of window movement)
  for (t in sim.nums) { #use this if you want the total age
  #for (i in 0:23) { #use this if you want just a subset of the age
    time.min <- t
    time.max <- t+1
    
    vic.counts <- subset(vicariance, time_bp>time.min & time_bp<time.max)
    num.vic <- nrow(vic.counts)
    
    founder.counts <- subset(founder, time_bp>time.min & time_bp<time.max)
    num.founder <- nrow(founder.counts)
    
    sympatry.counts <- subset(sym.y, time_bp>time.min & time_bp<time.max)
    num.sympatry <- nrow(sympatry.counts)
    
    subset.counts <- subset(sub.s, time_bp>time.min & time_bp<time.max)
    num.subset <- nrow(subset.counts)
    
    vic.found <- subset(vj, time_bp>time.min & time_bp<time.max)
    num.vj <- nrow(vic.found)
    
    sym.sub <- subset(ys, time_bp>time.min & time_bp<time.max)
    num.ys <- nrow(sym.sub)
    
    vic.sym <- subset(vys, time_bp>time.min & time_bp<time.max)
    num.vys <- nrow(vic.sym)
    
    all <- subset(total, time_bp>time.min & time_bp<time.max)
    num.all <- nrow(all)
    
    ratio.vj <- num.vj/num.all
    ratio.ys <- num.ys/num.all
    ratio.vys <- num.vys/num.all
    
    x <- matrix(nrow=1, ncol=3)
    #x[,1] <- num.all
    #x[,2] <- num.vj
    x[,1] <- ratio.vj
    x[,2] <- ratio.ys
    x[,3] <- ratio.vys
    #x[,4] <- num.ys
    #x[,5] <- ratio.ys
    #row.names(x)<-paste(time.min,'to', time.max)
    frame <- rbind(frame, x)
    
    c <- matrix(nrow=1, ncol=5)
    c[,1] <- num.all
    c[,2] <- num.vic
    c[,3] <- num.founder
    c[,4] <- num.sympatry
    c[,5] <- num.subset
    #row.names(x)<-paste(time.min,'to', time.max)
    extra.butt <- rbind(extra.butt, c)
  }
  #colnames(frame)<-c('num.all', 'num.vj', 'ratio.vj', 'num.ys', 'ratio.ys')
  frame <- as.data.frame(frame)
  full.table.vj  <- cbind(full.table.vj,  frame[,1]) # save the vj ratios
  full.table.ys  <- cbind(full.table.ys,  frame[,2]) # save the ys ratios
  full.table.vys <- cbind(full.table.vys, frame[,3]) # save the vys ratios
  full.table.vj  <- as.data.frame(full.table.vj)
  full.table.ys  <- as.data.frame(full.table.ys)
  full.table.vys <- as.data.frame(full.table.vys)
  
  save(full.table.vj, file="SIMULATED.full.table.vj.RData") # save the vj progress externally
  save(full.table.ys, file="SIMULATED.full.table.ys.RData") # save the ys progress externally
  save(full.table.vys,file="SIMULATED.full.table.vys.RData") # save the vys progress externally
  
  extra.butt <- as.data.frame(extra.butt)
  all.events <- cbind(all.events, extra.butt[,1])
  bsm.vicariance <- cbind(bsm.vicariance, extra.butt[,2])
  bsm.founder <- cbind(bsm.founder, extra.butt[,3])
  bsm.sympatry <- cbind(bsm.sympatry, extra.butt[,4])
  bsm.subset <- cbind(bsm.subset, extra.butt[,5])
  
  save(all.events,     file="SIMULATED.all.events.RData") # save the all events progress externally
  save(bsm.vicariance, file="SIMULATED.bsm.vicariance.RData") # save the bsm V progress externally
  save(bsm.founder,    file="SIMULATED.bsm.founder.RData") # save the bsm J progress externally
  save(bsm.sympatry,   file="SIMULATED.bsm.sympatry.RData") # save the bsm S progress externally
  save(bsm.subset,     file="SIMULATED.bsm.subset.RData") # save the bsm Y progress externally
}

colnames(full.table.vj)  <- c(1:(length(full.table.vj[1,])))
colnames(full.table.ys)  <- c(1:(length(full.table.ys[1,])))
colnames(full.table.vys) <- c(1:(length(full.table.vys[1,])))
colnames(all.events) <-c(1:(length(all.events[1,])))
colnames(bsm.vicariance) <- c(1:(length(bsm.vicariance[1,])))
colnames(bsm.founder) <- c(1:(length(bsm.founder[1,])))
colnames(bsm.sympatry) <- c(1:(length(bsm.sympatry[1,])))
colnames(bsm.subset) <- c(1:(length(bsm.subset[1,])))

rnames <- NULL
for (i in sim.nums) {
  tmin <- i-0.5
  tmax <- i+0.5
  clock <- paste(tmin, 'to', tmax)
  rnames <- append(rnames, clock)
}
row.names(full.table.vj)<-rnames
row.names(full.table.ys)<-rnames
row.names(full.table.vys)<-rnames
row.names(all.events)<-rnames
row.names(bsm.vicariance)<-rnames
row.names(bsm.founder)<-rnames
row.names(bsm.sympatry)<-rnames
row.names(bsm.subset)<-rnames

full.table.vj[is.na(full.table.vj)] <- 0
full.table.ys[is.na(full.table.ys)] <- 0
full.table.vys[is.na(full.table.vys)] <- 0
all.events[is.na(all.events)] <- 0
bsm.vicariance[is.na(bsm.vicariance)] <- 0
bsm.founder[is.na(bsm.founder)] <- 0
bsm.sympatry[is.na(bsm.sympatry)] <- 0
bsm.subset[is.na(bsm.subset)] <- 0

write.table(full.table.vj,      file="BGB.Meliphagides.UPDATED.SIMULATIONS.VJ.Dispersal.Freq.txt",  sep="\t", row.names=T, quote=F)
write.table(full.table.ys,      file="BGB.Meliphagides.UPDATED.SIMULATIONS.YS.Dispersal.Freq.txt",  sep="\t", row.names=T, quote=F)
write.table(full.table.vys,     file="BGB.Meliphagides.UPDATED.SIMULATIONS.VYS.Dispersal.Freq.txt",  sep="\t", row.names=T, quote=F)
write.table(all.events,         file="BGB.Meliphagides.UPDATED.SIMULATIONS.all.events.txt", sep="\t", row.names=T, quote=F)
write.table(bsm.vicariance,     file="BGB.Meliphagides.UPDATED.SIMULATIONS.vicariance.txt", sep="\t", row.names=T, quote=F)
write.table(bsm.founder,        file="BGB.Meliphagides.UPDATED.SIMULATIONS.founder.txt", sep="\t", row.names=T, quote=F)
write.table(bsm.sympatry,       file="BGB.Meliphagides.UPDATED.SIMULATIONS.sympatry.txt", sep="\t", row.names=T, quote=F)
write.table(bsm.subset,         file="BGB.Meliphagides.UPDATED.SIMULATIONS.subset.txt", sep="\t", row.names=T, quote=F)


#############################################################################################
#############################################################################################
## Step 5: This is going to pull out our 95% confidence intervals, and store them into a DF.
## Then we can add the empirical data, and plot them both together, to get an idea of how 
## our empirical data compares to the simulated data through time. Our focus is on the Miocene.
#############################################################################################
#############################################################################################

########################################################################################
## Let's get our 95% confidence intervals to compare against our empirical data ##
########################################################################################

# Start with Vicariance (v) and Founder (j)
pass <- NULL
butt <- NULL
for (i in 1:(length(full.table.vj[,1]))) { #use this if you want the total age
  pass <- t.test(full.table.vj[i,])
  lower <- pass$conf.int[1]
  upper <- pass$conf.int[2]
  ci.mean <- pass$estimate
  
  x <- matrix(nrow=1, ncol=3)
  x[,1] <- lower
  x[,2] <- upper
  x[,3] <- ci.mean
  
  butt <- rbind(butt, x)
  butt[is.na(butt)] <- 0
}
colnames(butt) <- c("lowerCI", "upperCI", "meanCI")
complete.vj <- cbind(full.table.vj, butt)
#complete.vj <- cbind(butt2, empirical$ratio.vj) #call in the empirical data to the frame
write.table(complete.vj, file="BGB.Meliphagides.UPDATED.SIMULATIONS.VJ.with.CIs.txt", sep="\t", row.names=T, quote=F)

# Next we'll do Sympatry with (y) and (s)
pass <- NULL
butt <- NULL
for (i in 1:(length(full.table.ys[,1]))) { #use this if you want the total age
  pass <- t.test(full.table.ys[i,])
  lower <- pass$conf.int[1]
  upper <- pass$conf.int[2]
  ci.mean <- pass$estimate
  
  x <- matrix(nrow=1, ncol=3)
  x[,1] <- lower
  x[,2] <- upper
  x[,3] <- ci.mean
  
  butt <- rbind(butt, x)
  butt[is.na(butt)] <- 0
}
colnames(butt) <- c("lowerCI", "upperCI", "meanCI")
complete.ys <- cbind(full.table.ys, butt)
#complete.vj <- cbind(butt2, empirical$ratio.vj) #call in the empirical data to the frame
write.table(complete.ys, file="BGB.Meliphagides.UPDATED.SIMULATIONS.YS.with.CIs.txt", sep="\t", row.names=T, quote=F)

# And finally the Vicariance (v) and Sympatry (y & s)
pass <- NULL
butt <- NULL
for (i in 1:(length(full.table.vys[,1]))) { #use this if you want the total age
  pass <- t.test(full.table.vys[i,])
  lower <- pass$conf.int[1]
  upper <- pass$conf.int[2]
  ci.mean <- pass$estimate
  
  x <- matrix(nrow=1, ncol=3)
  x[,1] <- lower
  x[,2] <- upper
  x[,3] <- ci.mean
  
  butt <- rbind(butt, x)
  butt[is.na(butt)] <- 0
}
colnames(butt) <- c("lowerCI", "upperCI", "meanCI")
complete.vys <- cbind(full.table.vys, butt)
#complete.vj <- cbind(butt2, empirical$ratio.vj) #call in the empirical data to the frame
write.table(complete.vys, file="BGB.Meliphagides.UPDATED.SIMULATIONS.VYS.with.CIs.txt", sep="\t", row.names=T, quote=F)


# Set the ages for the windows you ran
simulation.dates <- seq(0, 20, 0.1) #here I've limited it to 20 million years, because that's our focal depth

# plot the ratio of V+J events to all events, through time
#frame$ratio.vj
#setlower <- complete.vj[1:20,"lowerCI"]
#setupper <- complete.vj[1:20, "upperCI"]
#setmean <- complete.vj[1:20, "meanCI"]
#set.all <- complete.vj[0:23, c("lowerCI", "upperCI", "meanCI")]
all.vj <- complete.vj[0:length(simulation.dates), c("lowerCI", "upperCI", "meanCI")]
all.ys <- complete.ys[0:length(simulation.dates), c("lowerCI", "upperCI", "meanCI")]
all.vys <- complete.vys[0:length(simulation.dates), c("lowerCI", "upperCI", "meanCI")]

#attempt<-data.frame(x=(1:((max(node.height(tr)))+1)), y=sety) # adjust the max x value to fit the age of the tree!
#lower.ci <- data.frame(x=(1:20), y=setlower) # adjust the max x value to fit the age of the tree!
#upper.ci <- data.frame(x=(1:20), y=setupper)
#mean.ci <- data.frame(x=(1:20), y=setmean)
#all.ci <- data.frame(x=(0:22), y=set.all)
sim.vj.ci <- data.frame(x=(simulation.dates), y=all.vj)
sim.ys.ci <- data.frame(x=(simulation.dates), y=all.ys)
sim.vys.ci <- data.frame(x=(simulation.dates), y=all.vys)


## Shade between the LOESS smoothed lines for UPPER and LOWER CIs
#################################################################
(ggplot(data=sim.vj.ci, mapping=aes(x=x, y=sim.vj.ci["y.lowerCI"])) 
  + geom_point()
  + geom_smooth()
  + geom_line()
  + scale_x_reverse())
par()

(ggplot(data=all.emp.ci)
  + geom_ribbon(aes(x=x, ymin=all.emp.ci$y.lowerCI, ymax=all.emp.ci$y.upperCI))
  #+ geom_smooth()
  + scale_x_reverse()
  + geom_smooth(method="auto", aes(x=x, y=all.emp.ci$y.meanCI), se=T))
####################################################################

## Shade between the LOESS smoothed lines for UPPER and LOWER CIs
#################################################################
g1 <- (ggplot(sim.vj.ci)
  + geom_smooth(aes(x=x, y=y.lowerCI), se=F)
  + geom_smooth(aes(x=x, y=y.upperCI), se=F))
gg1 <- ggplot_build(g1)
df2 <- data.frame(x=gg1$data[[1]]$x,
                  ymin=gg1$data[[1]]$y,
                  ymax=gg1$data[[2]]$y)
(g1 +  geom_ribbon(data=df2, aes(x=x, ymin=ymin, ymax=ymax),fill="grey")
  #+ geom_smooth(aes(x=x, y=y.empirical.ratio.vj), se=F)
  + scale_x_reverse()) # this second line is the empirical data
