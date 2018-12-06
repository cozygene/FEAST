FEAST
-----------------------

FEAST is a scalable method for investigating the compositional patterns  of complex microbial communities. By efficiently evaluating hundreds to thousands of potential source environments,  FEAST quantifies the fraction of these sources in a target microbial community (sink), as well as the contribution of an unknown, uncharacterized source (‘unknown’ source). FEAST is agnostic to the sequencing data type (i.e., 16s rRNA or metagenomic data).  FEAST is available in R. The input supported by FEAST is a file wherein the sources are rows, and the columns are OTU counts. The number of distinct OTUs should be equal for the both the sources and the sink.


There are seven analyses in the FEAST paper :

1. Simulation study (using realistic source distribution from the Earth Microbiome Project).
2. Sequencing depth simulation (using realistic source distribution from the Earth Microbiome Project).
3. Unknown source simulation (using realistic source distribution from Lax et al. 2014).
4. Baby time series (real data analysis using data from Backhed et al. 2015).
5. Home microbiome (real data analysis using data from Lax et al. 2014).
6. allo-HSCT (real data analysis using data from Taur et al. 2012).
7. ICU vs. Healthy adults (real data analysis using data from McDonald et al. 2016 and the American Gut Project).


© 2018 Big Data and Genomics Lab at UCLA All Rights Reserved


Software Requirements
-----------------------

- R version 3.4.4
- FEAST uses the following dependencies : "doParallel", "foreach", "mgcv", "reshape2", "ggplot2", "philentropy", "MCMCpack", "lsei", "Rcpp", "RcppArmadillo" and "cowplot".

Contact: liashenhav@gmail.com
