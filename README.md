FEAST
-----------------------

FEAST is an efficient Expectation-Maximization-based method to identify the sources and proportions of contamination in marker-gene and functional metagenomic studies. FEAST is currently in beta and is available in both R and Python. The data supported by FEAST must come in a file wherein the sources are rows, and the columns are OTU counts. The number of distinct OTUs should be equal for the both the sources and the sink. You can use the current code to run your own simulations, or modify the simulated sink to a real sink. 


Software Requirements:
-----------------------

- R version 3.4.0 or Python 2.7


FEAST in R
-----------------------

FEAST uses the following dependencies : "doParallel", "foreach", "mgcv", "reshape2", "ggplot2", "philentropy", "MCMCpack", "lsei", "Rcpp", "RcppArmadillo" and "cowplot".  



FEAST in Python 2.7 (pyFEAST)
-----------------------

pyFEAST mainly uses numpy and achieves faster performance than its R counterpart. The main difference between the R and python implementations of FEAST is that the R version initializes the mixing proportions of the EM algorithm using a solution to constrained least squares, and pyFEAST uses a solution from non-negative least squares. To our current knowledge, this does not cause differences in EM estimations.



Contact: liashenhav@gmail.com