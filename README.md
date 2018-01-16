# Miocene Australia
Associated files and scripts for manuscript on Miocene biome turnover, and phenotypic evolution of Australian vertebrates. This project focuses on shifts in biome dispersal history and body size evolution across five widely distributed Australian vertebrate radiations: marsupial mammals, meliphagoid birds, pygopodoid geckos, sphenomorphine skinks, and agamid lizards. If you're not sure what you need, check the "GitHub Directory..." excel spreadsheet for descriptions of the included files. And if you are looking for something that I've forgotten to put up, please email me at ian.brennan.anu.edu.au

Otherwise, this repository includes the following (organized under folders of the same names):  
## Data:
  ### Trees:  
   + Maximum clade credibility trees. (5x)  
   + 100 trees from dating analysis posteriors. (5x)     
  ### Phenotypic and Distributional Traits
   + Body size data (snout-vent length, body length, mass). (5x)
   + Distributional data, in binary format for BioGeoBEARS analysis. (5x)
   + Spreadsheet of raw body size and distributional data.
## Scripts:
  ### Biogeographic Analyses:
   + R code for BioGeoBEARS and Biogeographic Stochastic Mapping that iterates across a set of empirical trees, and summarizes and plots temporal trends in dispersal histories. 
   + R code for BioGeoBEARS and Biogeographic Stochastic Mapping that simulates data onto a specified tree under a specified historical biogeographic model, and summarizes and plots temporal trends in dispersal histories. 
  ### Modelling Phenotypic Evolution:
   + R code which includes the two novel mode-variable evolutionary models (Single-Rate Constraint, Two-Rate Constraint), built into the 'fitContinuous.paleo' function of Slater (2013). 
   + R code which iterates across set of empirical trees, comparatively fitting trait evolution models, saving information on fit and parameter estimates, and summarizes and plots results.
   + R code for simulating (1) trait data to test and plot the reliability of model estimates of shift time and alpha, and (2) the recoverability of novel mode-variable (SRC, TRC) models.  
## Figures:
   + Figure 1. Estimated shifts and biogeographic dispersal patterns through time, highlighting changes in morphological and biogeographic evolution in the Late Miocene.
   + Figure 2. Stacked barplots of comparative model fitting for each radiation, with focus on preference of mode-variable models.
   + Supplemental Figures 1-10.
## Outputs:
   + Comparative model fitting output for 100 posterior trees. This includes ALL compared models. (5x)
   + Comparative model fitting output for 100 posterior trees. This is a trimmed set of models (those in Fig.2), with recalculated AICcWt. (5x)  
   + 'Best' models for each tree in each radiation, as a result of comparative model fitting. (5x)
   + ALL empirical and simulated results of biogeographic analyses. (many)


