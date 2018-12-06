

#!/bin/bash
set -ex 

mkdir -p ../results/allo_HSCT
mkdir -p ../results/Baby_time_series 
mkdir -p ../results/Simulation 
mkdir -p ../results/Seq_depth_simulation
mkdir -p ../results/Unknown_source 
mkdir -p ../results/ICU_vs_healthy
mkdir -p ../results/HomeMicrobiome


# time # time a script 
Rscript "Main_allo_HSCT-analysis.R" >> ../results/allo_HSCT/allo_HSCT_output.txt&
Rscript "Main_Baby_time_series.R" >> ../results/Baby_time_series/Baby_time_series_output.txt&
Rscript "Main_HomeMicrobiome_analysis.R"  >> ../results/HomeMicrobiome/HomeMicrobiome_output.txt&
Rscript "Main_simulations.R"  >> ../results/Simulation/Simulation_output.txt&
Rscript "Main_seq_depth_simulation.R"  >> ../results/Seq_depth_simulation/Seq_depth_simulation_output.txt&
Rscript "Main_unknown_source_analysis.R"  >> ../results/Unknown_source/Unknown_source_output.txt&
Rscript "Main_ICU_vs_healthy_analysis.R"  >> ../results/ICU_vs_healthy/ICU_vs_healthy_output.txt&



wait