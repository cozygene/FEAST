Analysis_flag = 1 #set to 1 to run the analysis, 0 otherwise

#Analysis---
if(Analysis_flag == 1){
    
   
   # Analysis running-time is approx. 0.5 hours with st_flag_seq_depth = 0
   # Analysis running-time is approx. 8 hours with st_flag_seq_depth = 0
   
   ST_simulations_flag = 1  #set to 1 to run SourceTracker in addition to FEAST
   source("Run_simulation_study.R") 
}
  

Plot_flag = 1 #set to 1 to plot the results, 0 otherwise

#Plot---
if(Plot_flag == 1)
  source("Plot_simulation_results.R")