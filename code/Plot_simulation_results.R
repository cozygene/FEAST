library("reshape2")
library("ggplot2")
library("Rcpp")
library("vegan")
library("dplyr")
library("doParallel")
library("foreach")
library("mgcv")
library("cowplot")

setwd('/data/')
load(file = "plot_table_simulation.RData")
comp_plot<-ggplot(plot_table_simulation, aes(js_values_new)) +
  geom_errorbar(aes(ymin=st_R2_averages_new-st_sd, ymax=st_R2_averages_new+st_sd),
                width=.01) +
  geom_point(aes(y=st_R2_averages_new, colour="SourceTracker"), size = 1, shape = 16) +
  geom_line(aes(y=st_R2_averages_new, colour="SourceTracker")) +
  geom_point(aes(y=emnoise_R2_averages_new, colour="FEAST"), size = 1, shape = 16) +
  geom_errorbar(aes(ymin=emnoise_R2_averages_new-emnoise_sd, 
                    ymax=emnoise_R2_averages_new+emnoise_sd), 
                width=.01) +
  geom_line(aes(y=emnoise_R2_averages_new, colour="FEAST")) +  theme_bw() +ylim (0,1) + xlim(0,1)
comp_plot<- comp_plot + labs(colour="Method", y="R2 of proportion estimates", 
                                     x="Jensen shannon divergence")


ggsave(filename="/results/Simulation_plot.png", plot = comp_plot , dpi = 600, width = 8.75, height = 6.1, units = "in")
