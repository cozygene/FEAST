STENSL - Source Tracking with ENvironmental SeLection
-----------------------

Microbial source tracking analysis has emerged as a widespread technique for characterizing the properties of complex microbial communities. However, this analysis is currently limited to source environments sampled in a specific study. In order to expand the scope beyond one single study, and allow the ‘exploration’ of source environments using large databases and repositories such as the Earth Microbiome Project, a source selection procedure is required.  Such a procedure will allow differentiating between contributing environments and nuisance ones when the number of potential sources considered is high. Here, we introduce STENSL (microbial Source Tracking with ENvironment SeLection), a machine learning method that extends common microbial source tracking analysis by performing an unsupervised source selection and enabling sparse identification of latent source environments. By incorporating sparsity into the estimation of potential source environments, STENSL improves the accuracy of true source contribution, while significantly reducing the noise introduced by non-contributing ones. We therefore anticipate that source selection will augment microbial source tracking analyses, enabling exploration of multiple source environments from publicly available repositories while maintaining high accuracy of the statistical inference. 

For more details [url TBD].


Support
-----------------------

For support using this package, please email: ulzee@ucla.edu. This is our new beta version - your comments/insights would be greatly appreciated. 


Software Requirements and dependencies
-----------------------

*STENSL* is implemented in R (>= 3.4.4) and requires the following dependencies: **Rcpp**, **RcppArmadillo**, **vegan**, **dplyr**, **reshape2**, **gridExtra**, **ggplot2**, **ggthemes**, **glmnet**, **CVXR**. Please install and load them prior to trying to install *STENSL*. If you are using a mac and having installation issues with **Rcpp** and or **RcppArmadillo**, try installing homebrew or xcode then reinstalling **Rcpp** and **RcppArmadillo**. 


```
Packages <- c("Rcpp", "RcppArmadillo", "vegan", "dplyr", "reshape2", "gridExtra", "ggplot2", "ggthemes", "glmnet", "CVXR")
install.packages(Packages)
lapply(Packages, library, character.only = TRUE)

```


Installation
---------------------------

*STENSL* will be available on QIIME 2 very soon. Until then you can you can simply install *STENSL* using **devtools**: 
```
devtools::install_github("cozygene/FEAST@STENSL")
```

## Usage
As input, *STENSL* takes mandatory arguments:

- _C_ - An _m_ by _n_ count matrix, where _m_ is the number samples and _n_ is the number of taxa.
- _metadata_ - An _m_ by 3 table, where _m_ is the number of samples. The metadata table has three columns (i.e., 'Env', 'SourceSink', 'id'). The first column is a description of the sampled environment (e.g., human gut), the second column indicates if this sample is a source or a sink (can take the value 'Source' or 'Sink'). The third column is the Sink-Source id. When using multiple sinks, each tested with the same group of sources, only the rows with 'SourceSink' = Sink will get an id (between 1 - number of sinks in the data). In this scenario, the sources’ ids are blank. When using multiple sinks, each tested with a distinct group of sources, each combination of sink and its corresponding sources should get the same id (between 1 - number of sinks in the data). Note that these names must be respected.
- _EM_iterations_ - A numeric value indicating the number of EM iterations (default 1000).
- _COVERAGE_ - A numeric value indicating the rarefaction depth (default = minimal sequencing depth within each group of sink and its corresponding sources).

Value: 

*STENSL* returns a vector whose length is the number of sources (including an unknown source). Each element is the estimation proportion of contribution where all proportions sum to 1.

<!-- *STENSL* returns an S1 by S2 matrix P, where S1 is the number sinks and S2 is the number of sources (including an unknown source). Each row in matrix P sums to 1. Pij is the contribution of source j to sink i. If Pij == NA it indicates that source j was not used in the analysis of sink i. *FEAST* will save the file "demo_FEAST.txt" (a file containing matrix P) . -->




Demo
-----------------------
We provide a dataset for an example of STENSL usage. Download the demo files <a href="https://github.com/cozygene/FEAST/tree/STENSL/Data_files">here</a>.

STENSL is implemented as part of the FEAST package.
First load the **FEAST** packages into R which will the STENSL function:
```
library(FEAST)
```

Then, load the datasets:
```
meta = read.table('Data_files/metadata_example_stensl.txt', sep='\t', header=T)
otus = read.table('Data_files/otu_example_stensl.txt', sep='\t', header=T)
```

```
result <- STENSL(
	C=as.matrix(otus),
	metadata=meta,
	EM_iterations=MAX_ITERS,
	COVERAGE=COVERAGE_DEPTH,
	l.range=c(0.1,1,10)
)
```

*STENSL* uses the same inpurt format as *FEAST*.
Details of the input are described in https://github.com/cozygene/FEAST/blob/master/README.md

