Analysis_flag = 1 #set to 1 to run the analysis, 0 otherwise
  
  #Analysis---
  if(Analysis_flag == 1){
  
    ST_unknown_flag = 1  #Change to 1 to run SourceTracker 
    source("Run_unknown_source_simulation.R")
  }
    

Plot_flag = 1 #set to 1 to plot the results, 0 otherwise

  #Plot---
  if(Plot_flag == 1) 
   source("Plot_Unknown_source_simulation.R")

  
  
