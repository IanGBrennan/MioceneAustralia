## Scripts for Analysis of Phenotypic Evolution
   + _New Models Adapted from Slater 2013_:  
   Slater 2013 code which includes the 'fitContinuous_paleo' function, with new models (SRC, TRC) from this study added.   
   
   + _ModelFitting_PRSB_:  
   R code for loop of model fitting and comparison across a set of empirical trees drawn from the posterior of a dating analysis. You can then summarize and plot the results, to determine the best model-fit, and comparative AICcWt. Finally, you can pull out parameter estimates from the best fitting model, and plot these (including: rate/sigma.sq, theta/optimum trait value, alpha/constraint, and shift time).  
   
   + _Test Alpha and Timing_:  
   R code for simulating trait data to test the reliability of model estimates of the shift time and the alpha parameter. Given input tree(s), you can simulate trait data, then run loops to estimate the parameter of interest, and plot the correspence between simulated and estimated parameter values. This includes a function to plot the linear model outputs (slope, r.squared, p-value) on the figure. 
   
   + _Sim.Fossil.Trees.SOURCE.R_ and _Simulating.Extinct.Trees.Data_PRSB.R_:  
   R code ('Sim.Fossil.Trees.SOURCE.R' is source code for the 'Simulating.Extinct...R' script) for simulating extinction regimes onto empirical trees, then simulating data onto the trees including extinct taxa, and finally dropping all extinct taxa from the trees and data. If interested, we left the original code for simulating data under the BMOU and BMOUi models in the file _Simulating.Extinct.Trees.Data.LOOP.R_
