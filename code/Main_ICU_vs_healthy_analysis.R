Analysis_flag = 1 #set to 1 to run the analysis, 0 otherwise

# Analysis running-time is >> 10 hours (analyzing 100 ICU patients and healthy individuals)
#Analysis---
if(Analysis_flag == 1)
  source("Run_ICU_vs_healthy_analysis.R") 

Plot_flag = 1  #set to 1 to plot the results, 0 otherwise

#Plot---
if(Plot_flag == 1)
    source("Plot_ICU_vs_healthy.R")

  