Analysis_flag = 1 #set to 1 to run the analysis, 0 otherwise


#Analysis---
if(Analysis_flag == 1){

   # Analysis running-time is approx. 1.5 hours with st_flag_seq_depth = 0
   # Analysis running-time is approx. 10 hours with st_flag_seq_depth = 0
   
    st_flag_seq_depth = 1  #set to 1 to run SourceTracker in addition to FEAST 
    source("Run_seq_depth_simulation.R") 
}
  