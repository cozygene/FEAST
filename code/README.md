# FEAST


## Guidelines
There are seven analyses in the FEAST paper :

1. Simulation study (using realistic source distribution from the Earth Microbiome Project).
2. Sequencing depth simulation (using realistic source distribution from the Earth Microbiome Project).
3. Unknown source simulation (using realistic source distribution from Lax et al. 2014).
4. Baby time series (real data analysis using data from Backhed et al. 2015).
5. Home microbiome (real data analysis using data from Lax et al. 2014).
6. allo-HSCT (real data analysis using data from Taur et al. 2012).
7. ICU vs. Healthy adults (real data analysis using data from McDonald et al. 2016 and the American Gut Project).

### How to run the code?

To run all the analyses in FEAST, set 'run.sh' as main and click on 'Run'. 

### How to run each analysis separately?

For each analysis in the paper, you will find a 'Main...' file, which will 
enable you to re-run the analysis and/or plot the results.
In order to run each analysis:

1. Right click on the relevant analysis (e.g., 'Main_simulations.R') 
and choose 'Set as Main' from the Actions menu.
2. For some analyses, you can choose whether to run it from scratch
(Analysis_Flag == 1) and/or to plot the results (Plot_flag == 1). 
3. For some analyses, you can choose whether to compare to other methods i.e.,  
(Compare_Flag == 1). 
2. Click Run. 
3. The output and plot will be saved automatically.

### comments

1. In the ICU code, we used a smaller number of sinks (20) since the  original analysis is very lengthy (over 24 hours).

© 2018 Big Data and Genomics Lab at UCLA All Rights Reserved