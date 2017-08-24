# to run the BSMs you'll have to go all the way through the BGB analysis. So speed it up, you
# can just do the favored model (DEC or whatever) and get it into a 'res' element for the BSM
# then the BSMs start down towards line 760

library(optimx)   
library(FD)       
library(parallel)
library(BioGeoBEARS)
library(snow)
library(ggplot2)
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
setwd("/YOUR_DIRECTORY")
wd <- getwd()

# Just a little TREE preparation:
#################################
#BGB doesn't like using multiphylo objects, so we'll make each tree from the 
#posterior set an individual tree, and then run the analysis
#################################
trees = read.nexus("Data.TREES/PB.Meliphagides.100.trees") #pick out our set of posterior trees
#trees <- lapply(trees, drop.tip, "A_Physignathus_cocincinus") #if you have to drop tips from the set of trees
#class(trees)<-"multiPhylo"

# Run a loop that writes each one to an individual file in a named folder
dir.create(np(paste(addslash(wd), "Meliphagides.PB.Trees", sep=""))) #create folder for the trees
for (i in 1:100){
  name <- paste("/YOUR_DIRECTORY/PB.Meliphagides.",i,".tre", sep="")
  write.tree(trees[[i]], file=name)
} #you should now have a folder with 100 tree separate tree files



#############################################################################################
#############################################################################################
## Step 1: Start with a giant loop which runs the simulated data under the DEC & DEC+j models,
## follows this with a 50 BSMs, and then summarizes the results of the BSMs.
## The loop will iterate through
## 100 times, after which we'll have to extract our 95% confidence intervals. This
## loop can take a WHILE, so check emp.progress by opening the "emp.progress.txt" document.
#############################################################################################
#############################################################################################

#let's make a loop for all the simulated things
emp.frame <- NULL
emp.full.table.vj <- NULL; emp.full.table.ys <- NULL; emp.full.table.vys <- NULL
emp.bsm.vicariance <- NULL
emp.bsm.founder <- NULL
emp.bsm.sympatry <- NULL
emp.bsm.subset <- NULL
emp.all.events <- NULL
emp.progress <- NULL
#pb <- emp.progress_bar$new(total = 100)

# Input Data
################################
# before we start the loop, we can designate the tip states (they won't change)
geogfn = np(paste(addslash(wd), "Data.GEOG.files/BGB.Meliphagides.geog.txt", sep="")) #empirical data
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn) # Look at your geographic range data:
max_range_size = 5 # Set the maximum number of areas any species may occupy

# we have to wait to drop in the trees because of BGB's tree set-up

for (t in 1:100) {
  # Make a counter for the loops so we know how far we are
  ########################################################
  #pb$tick()
  #Sys.sleep(1 / 100)
  emp.progress <- print(t)
  write(emp.progress, file="emp.progress.txt", append=T) #remember to trash the old emp.progress file before starting!
  
  # and we can drop in the trees
  next.tree <- paste("Meliphagides.PB.Trees/PB.Meliphagides.",t,".tre", sep="") #direct to the folder and tree files
  #next.tree <- paste("BGB.Liolaemidae.MCC.tre")
  trfn = np(paste(addslash(wd), next.tree, sep="")) #import the tree file
  tr = read.tree(trfn) #can't use read.nexus, so make sure to change nexus to newick trees

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
  resfn = "Emp_DEC_M0_unconstrained_v1.Rdata"
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
  
  resfn = "sim_DEC+J_M0_unconstrained_v1.Rdata"
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
  sym.y <- NULL
  for (i in 1:length(clado_events_tables)) {
    symy<-subset(clado_events_tables[[i]], clado_events_tables[[i]]$clado_event_type == "sympatry (y)")
    sym.y<-rbind(sym.y, symy)
  }
  
  # pull out the Subset (s) dispersal events
  sub.s <- NULL
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
  emp.frame <- NULL
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
    emp.frame <- rbind(emp.frame, x)
    
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
  emp.frame <- as.data.frame(emp.frame)
  emp.full.table.vj <- cbind(emp.full.table.vj, emp.frame[,1]) # save the vj ratios
  emp.full.table.ys <- cbind(emp.full.table.ys, emp.frame[,2]) # save the ys ratios
  emp.full.table.vys <- cbind(emp.full.table.vys, emp.frame[,3]) # save the vys ratios
  emp.full.table.vj <- as.data.frame(emp.full.table.vj)
  emp.full.table.ys <- as.data.frame(emp.full.table.ys)
  emp.full.table.vys <- as.data.frame(emp.full.table.vys)
  
  save(emp.full.table.vj, file="EMPIRICAL.emp.full.table.vj.RData") # save the vj progress externally
  save(emp.full.table.ys, file="EMPIRICAL.emp.full.table.ys.RData") # save the ys progress externally
  save(emp.full.table.vys, file="EMPIRICAL.emp.full.table.vys.RData") # save the vys progress externally
  
  extra.butt <- as.data.frame(extra.butt)
  emp.all.events <- cbind(emp.all.events, extra.butt[,1])
  emp.bsm.vicariance <- cbind(emp.bsm.vicariance, extra.butt[,2])
  emp.bsm.founder <- cbind(emp.bsm.founder, extra.butt[,3])
  emp.bsm.sympatry <- cbind(emp.bsm.sympatry, extra.butt[,4])
  emp.bsm.subset <- cbind(emp.bsm.subset, extra.butt[,5])
  
  save(emp.all.events, file="EMPIRICAL.emp.all.events.RData") # save the all events progress externally
  save(emp.bsm.vicariance, file="EMPIRICAL.emp.bsm.vicariance.RData") # save the bsm V progress externally
  save(emp.bsm.founder, file="EMPIRICAL.emp.bsm.founder.RData") # save the bsm J progress externally
  save(emp.bsm.sympatry, file="EMPIRICAL.emp.bsm.sympatry.RData") # save the bsm S progress externally
  save(emp.bsm.subset, file="EMPIRICAL.emp.bsm.subset.RData") # save the bsm Y progress externally
}

colnames(emp.full.table.vj) <- c(1:(length(emp.full.table.vj[1,])))
colnames(emp.full.table.ys) <- c(1:(length(emp.full.table.ys[1,])))
colnames(emp.full.table.vys) <- c(1:(length(emp.full.table.vys[1,])))
colnames(emp.all.events) <-c(1:(length(emp.all.events[1,])))
colnames(emp.bsm.vicariance) <- c(1:(length(emp.bsm.vicariance[1,])))
colnames(emp.bsm.founder) <- c(1:(length(emp.bsm.founder[1,])))
colnames(emp.bsm.sympatry) <- c(1:(length(emp.bsm.sympatry[1,])))
colnames(emp.bsm.subset) <- c(1:(length(emp.bsm.subset[1,])))

rnames <- NULL
for (i in sim.nums) {
  tmin <- i-0.5
  tmax <- i+0.5
  clock <- paste(tmin, 'to', tmax)
  rnames <- append(rnames, clock)
}
row.names(emp.full.table.vj)<-rnames
row.names(emp.full.table.ys)<-rnames
row.names(emp.full.table.vys)<-rnames
row.names(emp.all.events)<-rnames
row.names(emp.bsm.vicariance)<-rnames
row.names(emp.bsm.founder)<-rnames
row.names(emp.bsm.sympatry)<-rnames
row.names(emp.bsm.subset)<-rnames

emp.full.table.vj[is.na(emp.full.table.vj)] <- 0
emp.full.table.ys[is.na(emp.full.table.ys)] <- 0
emp.full.table.vys[is.na(emp.full.table.vys)] <- 0
emp.all.events[is.na(emp.all.events)] <- 0
emp.bsm.vicariance[is.na(emp.bsm.vicariance)] <- 0
emp.bsm.founder[is.na(emp.bsm.founder)] <- 0
emp.bsm.sympatry[is.na(emp.bsm.sympatry)] <- 0
emp.bsm.subset[is.na(emp.bsm.subset)] <- 0

write.table(emp.full.table.vj,  file="BGB.Meliphagides.UPDATED.EMPIRICAL.VJ.Dispersal.Freq.txt",  sep="\t", row.names=T, quote=F)
write.table(emp.full.table.ys,  file="BGB.Meliphagides.UPDATED.EMPIRICAL.YS.Dispersal.Freq.txt",  sep="\t", row.names=T, quote=F)
write.table(emp.full.table.vys, file="BGB.Meliphagides.UPDATED.EMPIRICAL.VYS.Dispersal.Freq.txt",  sep="\t", row.names=T, quote=F)
write.table(emp.all.events,     file="BGB.Meliphagides.UPDATED.EMPIRICAL.emp.all.events.txt", sep="\t", row.names=T, quote=F)
write.table(emp.bsm.vicariance, file="BGB.Meliphagides.UPDATED.EMPIRICAL.vicariance.txt", sep="\t", row.names=T, quote=F)
write.table(emp.bsm.founder,    file="BGB.Meliphagides.UPDATED.EMPIRICAL.founder.txt", sep="\t", row.names=T, quote=F)
write.table(emp.bsm.sympatry,   file="BGB.Meliphagides.UPDATED.EMPIRICAL.sympatry.txt", sep="\t", row.names=T, quote=F)
write.table(emp.bsm.subset,     file="BGB.Meliphagides.UPDATED.EMPIRICAL.subset.txt", sep="\t", row.names=T, quote=F)


#############################################################################################
#############################################################################################
## Step 5: This is going to pull out our 95% confidence intervals, and store them into a DF.
## Then we can add the empirical data, and plot them both together, to get an idea of how 
## our empirical data compares to the simulated data through time. Our focus is on the Miocene.
#############################################################################################
#############################################################################################

#############################################################################################################
## Let's get our 95% confidence intervals for VJ, YS, and VYS ratios to compare against our empirical data ##
#############################################################################################################

# Start with Vicariance (v) and Founder (j)
###########################################
# in case you need to reload a file:
#emp.full.table.vj <- read.csv("BGB.____________.UPDATED.EMPIRICAL.VJ.Dispersal.Freq.txt", header=T, sep="\t", row.names=1) # row.names may or may not be necessary
emp.full.table.vj <- read.csv("/YOUR_DIRECTORY/BGB.All.Radiations.UPDATED.SIMULATED.vj.ratio.txt", header=T, sep="\t", row.names=1) # row.names may or may not b
#emp.full.table.vj[is.na(emp.full.table.vj)] <- 0 #use this if too many NAs are introduced
###########################################
pass <- NULL
butt <- NULL
for (i in 1:(length(emp.full.table.vj[,1]))) { #use this if you want the total age
  pass <- t.test(emp.full.table.vj[i,])
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
complete.vj <- cbind(emp.full.table.vj, butt)
#complete.vj <- cbind(butt2, empirical$ratio.vj) #call in the empirical MCC data to the frame
write.table(complete.vj, file="BGB.Meliphagides.UPDATED.EMPIRICAL.VJ.with.CIs.txt", sep="\t", row.names=T, quote=F)

# Next we'll do Sympatry with (y) and (s)
###########################################
# in case you need to reload the file:
#emp.full.table.vj <- read.csv("BGB.__________.UPDATED.EMPIRICAL.YS.Dispersal.Freq.txt", header=T, sep="\t")
###########################################
pass <- NULL
butt <- NULL
for (i in 1:(length(emp.full.table.ys[,1]))) { #use this if you want the total age
  pass <- t.test(emp.full.table.ys[i,])
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
complete.ys <- cbind(emp.full.table.ys, butt)
#complete.vj <- cbind(butt2, empirical$ratio.vj) #call in the empirical MCC data to the frame
write.table(complete.ys, file="BGB.Meliphagides.UPDATED.EMPIRICAL.YS.with.CIs.txt", sep="\t", row.names=T, quote=F)

# And finally the Vicariance (v) and Sympatry (y & s)
# Next we'll do Sympatry with (y) and (s)
###########################################
# in case you need to reload the file:
#emp.full.table.vj <- read.csv("BGB.__________.UPDATED.EMPIRICAL.VYS.Dispersal.Freq.txt", header=T, sep="\t")
###########################################
pass <- NULL
butt <- NULL
for (i in 1:(length(emp.full.table.vys[,1]))) { #use this if you want the total age
  pass <- t.test(emp.full.table.vys[i,])
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
complete.vys <- cbind(emp.full.table.vys, butt)
#complete.vj <- cbind(butt2, empirical$ratio.vj) #call in the empirical MCC data to the frame
write.table(complete.vys, file="BGB.Meliphagides.UPDATED.EMPIRICAL.VYS.with.CIs.txt", sep="\t", row.names=T, quote=F)

# Set the ages for the windows you ran
simulation.dates <- seq(0, 20, 0.1) #here I've limited it to 20 million years, because that's our focal depth

# plot the ratio of V+J events to all events, through time
#frame$ratio.vj
#setlower <- complete.vj[1:20,"lowerCI"]
#setupper <- complete.vj[1:20, "upperCI"]
#setmean <- complete.vj[1:20, "meanCI"]
#set.all <- complete.vj[0:23, c("lowerCI", "upperCI", "meanCI")]
all.emp.vj <- complete.vj[0:length(simulation.dates), c("lowerCI", "upperCI", "meanCI")]
all.emp.ys <- complete.ys[0:length(simulation.dates), c("lowerCI", "upperCI", "meanCI")]
all.emp.vys <- complete.vys[0:length(simulation.dates), c("lowerCI", "upperCI", "meanCI")]

#attempt<-data.frame(x=(1:((max(node.height(tr)))+1)), y=sety) # adjust the max x value to fit the age of the tree!
#lower.ci <- data.frame(x=(1:20), y=setlower) # adjust the max x value to fit the age of the tree!
#upper.ci <- data.frame(x=(1:20), y=setupper)
#mean.ci <- data.frame(x=(1:20), y=setmean)
#all.ci <- data.frame(x=(0:22), y=set.all)
all.emp.ci.vj <- data.frame(x=(simulation.dates), y=all.emp.vj)
all.emp.ci.ys <- data.frame(x=(simulation.dates), y=all.emp.ys)
all.emp.ci.vys <- data.frame(x=(simulation.dates), y=all.emp.vys)

## Shade between the LOESS smoothed lines for UPPER and LOWER CIs
#################################################################
(ggplot(data=all.emp.ci.vj, mapping=aes(x=x, y=all.emp.ci.vj["y.lowerCI"])) 
 + geom_point()
 + geom_smooth()
 + geom_line()
 + scale_x_reverse())
par()

(ggplot(data=all.emp.ci.vj)
  + geom_ribbon(aes(x=x, ymin=all.emp.ci.vj$y.lowerCI, ymax=all.emp.ci.vj$y.upperCI))
  #+ geom_smooth()
  + scale_x_reverse()
  + geom_smooth(method="auto", aes(x=x, y=all.emp.ci.vj$y.meanCI), se=T))

## Shade between the LOESS smoothed lines for UPPER and LOWER CIs
#################################################################
g1 <- (ggplot(all.emp.ci.ys)
       + geom_smooth(aes(x=x, y=y.lowerCI), se=F)
       + geom_smooth(aes(x=x, y=y.upperCI), se=F))
gg1 <- ggplot_build(g1)
df2 <- data.frame(x=gg1$data[[1]]$x,
                  ymin=gg1$data[[1]]$y,
                  ymax=gg1$data[[2]]$y)
(g1 +  geom_ribbon(data=df2, aes(x=x, ymin=ymin, ymax=ymax),fill="grey")
  #+ geom_smooth(aes(x=x, y=y.empirical.ratio.vj), se=F) # this is 2nd line is from the MCC tree
  + scale_x_reverse()) # uncomment the above if you want to plot it


#########################################################################
# Now let's plot the Empirical CI against the Simulated CI
#########################################################################
simulation.vj <- read.csv("BGB.Agamidae.UPDATED.SIMULATIONS.vj.with.CIs.xls", header=T, sep="\t")
emp.v.sim <- cbind(complete.vj, simulation.vj$lowerCI) #call in the empirical data to the frame
emp.v.sim <- cbind(emp.v.sim, simulation.vj$upperCI)

# Set the ages for the windows you ran
simulation.dates <- seq(0, 20, 0.1) #here I've limited it to 20 million years, because that's our focal depth

# combine the SIMULATED and EMPIRICAL data
all.emp.sim <- emp.v.sim[0:length(simulation.dates), c("lowerCI", "upperCI", "meanCI", 
                                                       "empirical$ratio.vj", "simulation.vj$lowerCI",
                                                       "simulation.vj$upperCI")]

#attempt<-data.frame(x=(1:((max(node.height(tr)))+1)), y=sety) # adjust the max x value to fit the age of the tree!
all.emp.sim.ci <- data.frame(x=(simulation.dates), y=all.emp.sim)

## Shade between the LOESS smoothed lines for UPPER and LOWER CIs
#################################################################
g1 <- (ggplot(all.emp.sim.ci) #this is for the empirical data
       + geom_smooth(aes(x=x, y=y.lowerCI), se=F)
       + geom_smooth(aes(x=x, y=y.upperCI), se=F))
gg1 <- ggplot_build(g1)

g2 <- (ggplot(all.emp.sim.ci) #this is for the simulated data
       + geom_smooth(aes(x=x, y=y.simulation.vj.lowerCI), se=F)
       + geom_smooth(aes(x=x, y=y.simulation.vj.upperCI), se=F))
gg2 <- ggplot_build(g2)

df2 <- data.frame(x.emp=gg1$data[[1]]$x,
                  ymin.emp=gg1$data[[1]]$y,
                  ymax.emp=gg1$data[[2]]$y,
                  x.sim=gg2$data[[1]]$x,
                  ymin.sim=gg2$data[[1]]$y,
                  ymax.sim=gg2$data[[2]]$y)

(g2
  + geom_ribbon(data=df2, aes(x=x.sim, ymin=ymin.sim, ymax=ymax.sim),fill="navy", alpha=0.75) #simulated first
  + geom_ribbon(data=df2, aes(x=x.emp, ymin=ymin.emp, ymax=ymax.emp),fill="pink", alpha=0.75) #empirical second
  + geom_smooth(aes(x=x, y=y.empirical.ratio.vj), se=F, color="salmon")
  + scale_x_reverse()) # this second line is the empirical data
