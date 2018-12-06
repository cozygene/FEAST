FEAST
-----------------------

FEAST is a scalable method for investigating the compositional patterns  of complex microbial communities. By efficiently evaluating hundreds to thousands of potential source environments,  FEAST quantifies the fraction of these sources in a target microbial community (sink), as well as the contribution of an unknown, uncharacterized source (‘unknown’ source). FEAST is available in R and agnostic to the sequencing data type (i.e., 16s rRNA or metagenomic data). The input supported by FEAST is either .txt or .csv file wherein the sources are rows, and the columns are the read counts. The number of distinct OTUs should be equal for the both the sources and the sink.

© 2018 Big Data and Genomics Lab at UCLA All Rights Reserved


Software Requirements
-----------------------

- R version 3.4.4
- FEAST uses the following dependencies : "doParallel", "foreach", "mgcv", "reshape2", "ggplot2", "philentropy", "MCMCpack", "lsei", "Rcpp", "RcppArmadillo" and "cowplot".

Contact: liashenhav@gmail.com
