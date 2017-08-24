library(ggplot2)

#########################################################################
# Now let's plot the Empirical CI against the Simulated CI [VICARIANCE/JUMP]
#########################################################################

## read in the simulated data first
#sim.vj <- read.csv("Data.SIMULATED/BGB.Meliphagides.UPDATED.SIMULATIONS.VJ.with.CIs.txt", header=T, sep="\t")
sim.vj <- read.csv("YOUR_DIRECTORY/BGB.All.Radiations.UPDATED.SIMULATED.vj.ratio.CIs.txt", header=T, sep="\t")

## read in the empirical data next
#emp.vj <- as.data.frame(read.csv("Data.EMPIRICAL/BGB.Meliphagides.UPDATED.EMPIRICAL.VJ.with.CIs.txt", header=T, sep="\t"))
emp.vj <- as.data.frame(read.csv("YOUR_DIRECTORY/BGB.All.Radiations.UPDATED.EMPIRICAL.vj.ratio.CIs.txt", header=T, sep="\t"))

# combine the SIMULATED and EMPIRICAL data
all.emp.sim <- NULL
all.emp.sim <- cbind.data.frame(sim.vj$lowerCI, sim.vj$upperCI, sim.vj$meanCI, 
                                emp.vj$lowerCI, emp.vj$upperCI, emp.vj$meanCI)

# Set the ages for the windows you ran, and trim the data set
simulation.dates <- seq(0, 20, 0.1) #here I've limited it to 20 million years, because that's our focal depth
all.emp.sim <- all.emp.sim[0:length(simulation.dates), c("sim.vj$lowerCI", 
                                                         "sim.vj$upperCI",
                                                         "sim.vj$meanCI",
                                                         "emp.vj$lowerCI", 
                                                         "emp.vj$upperCI",
                                                         "emp.vj$meanCI")]
                                                         #"emp.vj$empirical.ratio.vj")]

# adjust the max x value to fit the age of the tree!
all.emp.sim.ci <- data.frame(x=(simulation.dates), y=all.emp.sim)

## Shade between the LOESS smoothed lines for UPPER and LOWER CIs
#################################################################
g1 <- (ggplot(all.emp.sim.ci) #this is for the empirical data
       + geom_smooth(aes(x=x, y=y.emp.vj.lowerCI), se=F)
       + geom_smooth(aes(x=x, y=y.emp.vj.upperCI), se=F))
gg1 <- ggplot_build(g1)

g2 <- (ggplot(all.emp.sim.ci) #this is for the simulated data
       + geom_smooth(aes(x=x, y=y.sim.vj.lowerCI), se=F, col="navy", alpha=0.6)
       + geom_smooth(aes(x=x, y=y.sim.vj.upperCI), se=F, col="navy", alpha=0.6))
gg2 <- ggplot_build(g2)

df2 <- data.frame(x.emp=gg1$data[[1]]$x,
                  ymin.emp=gg1$data[[1]]$y,
                  ymax.emp=gg1$data[[2]]$y,
                  x.sim=gg2$data[[1]]$x,
                  ymin.sim=gg2$data[[1]]$y,
                  ymax.sim=gg2$data[[2]]$y)
(g2
  + geom_ribbon(data=df2, aes(x=x.sim, ymin=ymin.sim, ymax=ymax.sim),fill="navy", alpha=0.6) #simulated first
  + geom_ribbon(data=df2, aes(x=x.emp, ymin=ymin.emp, ymax=ymax.emp),fill="pink", alpha=0.75) #empirical second
  #+ geom_smooth(aes(x=x, y=y.emp.vj.empirical.ratio.vj), se=F, color="salmon")
  + scale_x_reverse()) # this second line is the empirical data




(ggplot(data=all.emp.sim.ci)
  + geom_ribbon(aes(x=x, ymin=all.emp.sim.ci$y.sim.vj.lowerCI, ymax=all.emp.sim.ci$y.sim.vj.upperCI))
  #+ geom_smooth()
  + scale_x_reverse())
  #+ geom_smooth(method="auto", aes(x=x, y=all.emp.ci.vj$y.meanCI), se=T))

(ggplot(data=all.emp.sim.ci)
  + geom_ribbon(aes(x=x, ymin=all.emp.sim.ci$y.emp.vj.lowerCI, ymax=all.emp.sim.ci$y.emp.vj.upperCI))
  + geom_ribbon(aes(x=x, ymin=all.emp.sim.ci$y.sim.vj.lowerCI, ymax=all.emp.sim.ci$y.sim.vj.upperCI))
  #+ geom_smooth()
  + scale_x_reverse()
  + geom_smooth(method="auto", aes(x=x, y=all.emp.sim.ci$y.emp.vj.meanCI), se=T)
  + geom_smooth(method="auto", aes(x=x, y=all.emp.sim.ci$y.sim.vj.meanCI), se=T)
  + geom_line(aes(x=x, y=all.emp.sim.ci$y.emp.vj.meanCI)))







#########################################################################
# Now let's plot the Empirical CI against the Simulated CI [SYMPATRY]
#########################################################################

## read in the simulated data first
sim.ys <- read.csv("Data.SIMULATED/BGB.Meliphagides.UPDATED.SIMULATIONS.YS.with.CIs.txt", header=T, sep="\t")
#sim.vj <- read.csv("Data.SIMULATED/BGB.All.Radiations.UPDATED.SIMULATIONS.VJ.with.CIs.txt", header=T, sep="\t")

## read in the empirical data next
emp.ys <- as.data.frame(read.csv("Data.EMPIRICAL/BGB.Meliphagides.UPDATED.EMPIRICAL.YS.with.CIs.txt", header=T, sep="\t"))
#emp.vj <- as.data.frame(read.csv("Data.EMPIRICAL/BGB.All.Radiations.ASGMB.UPDATED.EMPIRICAL.VJ.with.CIs.txt", header=T, sep="\t"))

# combine the SIMULATED and EMPIRICAL data
all.emp.sim <- NULL
all.emp.sim <- cbind.data.frame(sim.ys$lowerCI, sim.ys$upperCI, emp.ys$lowerCI, emp.ys$upperCI)

# Set the ages for the windows you ran, and trim the data set
simulation.dates <- seq(0, 20, 0.1) #here I've limited it to 20 million years, because that's our focal depth
all.emp.sim <- all.emp.sim[0:length(simulation.dates), c("sim.ys$lowerCI", 
                                                         "sim.ys$upperCI", 
                                                         "emp.ys$lowerCI", 
                                                         "emp.ys$upperCI")]
#"emp.vj$empirical.ratio.vj")]

# adjust the max x value to fit the age of the tree!
all.emp.sim.ci <- data.frame(x=(simulation.dates), y=all.emp.sim)

## Shade between the LOESS smoothed lines for UPPER and LOWER CIs
#################################################################
g1 <- (ggplot(all.emp.sim.ci) #this is for the empirical data
       + geom_smooth(aes(x=x, y=y.emp.ys.lowerCI), se=F)
       + geom_smooth(aes(x=x, y=y.emp.ys.upperCI), se=F))
gg1 <- ggplot_build(g1)

g2 <- (ggplot(all.emp.sim.ci) #this is for the simulated data
       + geom_smooth(aes(x=x, y=y.sim.ys.lowerCI), se=F, col="navy", alpha=0.6)
       + geom_smooth(aes(x=x, y=y.sim.ys.upperCI), se=F, col="navy", alpha=0.6))
gg2 <- ggplot_build(g2)

df2 <- data.frame(x.emp=gg1$data[[1]]$x,
                  ymin.emp=gg1$data[[1]]$y,
                  ymax.emp=gg1$data[[2]]$y,
                  x.sim=gg2$data[[1]]$x,
                  ymin.sim=gg2$data[[1]]$y,
                  ymax.sim=gg2$data[[2]]$y)
(g2
  + geom_ribbon(data=df2, aes(x=x.sim, ymin=ymin.sim, ymax=ymax.sim),fill="navy", alpha=0.6) #simulated first
  + geom_ribbon(data=df2, aes(x=x.emp, ymin=ymin.emp, ymax=ymax.emp),fill="pink", alpha=0.75) #empirical second
  #+ geom_smooth(aes(x=x, y=y.emp.vj.empirical.ratio.vj), se=F, color="salmon")
  + scale_x_reverse()) # this second line is the empirical data
