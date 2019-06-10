FEAST
-----------------------

FEAST is a scalable algorithm for investigating the compositional patterns  of complex microbial communities. 

Support
-----------------------

For support using FEAST, please email: liashenhav@gmail.com


Software Requirements and dependencies
-----------------------

FEAST is written R. In addition to R 3.4.4 (and higher), it has the following dependencies::

"doParallel", "foreach", "mgcv", "reshape2", "ggplot2", "philentropy", "MCMCpack", "lsei", "Rcpp", "RcppArmadillo" and "cowplot".

Input format
-----------------------
The input to FEAST are two tab-separated ASCII text files :

count table  - A matrix of samples by taxa with the sources and sink. The first row contains the sample headers (SampleID). The first column contains taxa ids. Then every consecutive column contains read counts for each sample. Note that this order must be respected (see example below).

Metadata -  The first row contains the headers (SampleID, Env, SourceSink). The first column contains gene ids. Then every consecutive pair of columns contains read counts for reference alleles and alternate alleles, respectively. Note that this order must be respected  (see example below).



Output format
-----------------------

The output is a vector of  contributions of the known and unknown sources (with the pre-defined source environments as headers). 

Usage instructions
---------------------------

FEAST will be available on Qiime II in July 2019. Until then you can easily run it on your computer in just a few easy steps which I will walk you through in the following lines. 

	1. Clone this repository ('FEAST') and save it on your computer.
	2. Save your input data (metadata and count table) in the directory 'Data_files'.
	3. Run the file 'FEAST_main' from 'FEAST_src' after inserting the following arguments as input:

path = the path in which you saved the directory 'FEAST_src' (string)
metadata_file =  your_metadata_file_name (string)
count_matrix =  your taxa count matrix (string)
num_sources <- number of sources (int)




Example
---------------------------

To run FEAST on example data (using multiple sinks) do:

	
	1. Clone this repository ('FEAST') and save it on your computer.
	2. Run the file 'FEAST_example' which takes the following arguments as input:
	path = the path in which you saved the directory 'FEAST' (string)
	

Input - 

Metadata:



Count matrix:

 

Output - 



© 2018 Big Data and Genomics Lab at UCLA All Rights Reserved
