Analysis_flag = 1 #set to 1 to run the analysis, 0 otherwise

# Analysis running-time is approx. 0.5 hours 
  
  #Analysis---
  if(Analysis_flag == 1)
    source("Run_HomeMicrobiome_analysis.R")

Plot_flag = 1 #set to 1 to plot the results, 0 otherwise

  #Plot---
  if(Plot_flag == 1) 
   source("Plot_Home_Microbiome.R")

  
  
