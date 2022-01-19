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

First load the **FEAST** packages into R which will also import STENSL:
```
library(FEAST)
```

Then, load the datasets:
```
metadata <- Load_metadata(metadata_path = "~/FEAST/Data_files/metadata_example_multi.txt")
otus <- Load_CountMatrix(CountMatrix_path = "~/FEAST/Data_files/otu_example_multi.txt")
```
Run _FEAST_, saving the output with prefix "demo":

```
FEAST_output <- FEAST(C = otus, metadata = metadata, different_sources_flag = 1, dir_path = "~/FEAST/Data_files/",
                      outfile="demo")
```

_FEAST_ will then save the file
*demo_FEAST.txt* - A file containing an S1 by S2 matrix P, where S1 is the number sinks and S2 is the number of sources (including an unknown source). Each row in matrix P sums to 1.

Graphical representation: 

As input, *PlotSourceContribution* takes mandatory arguments:

- _SinkNames_ - A vector with the sink names to plot.
- _SourceNames_ - A vector with all the sources' names.
- _Same_sources_flag_ - A Boolean value indicating the source-sink plotting assignment. Same_sources_flag = 1 if the same sources are assigned to the pre-defined sink samples , otherwise = 0.
- _dir_path_ - A path to an output .png file.
- _mixing_proportions_ - A list of vectors, where entry i corresponds to the vector of source contributions (summing to 1) to sink i.
- _Plot_title_ -  Plot's title and output .png file's name.
- _N_ - Number of barplots in each output .png file.


```
PlotSourceContribution(SinkNames = rownames(FEAST_output)[c(5:8)],
                       SourceNames = colnames(FEAST_output), dir_path = "~/FEAST/Data_files/",
                       mixing_proportions = FEAST_output, Plot_title = "Test_",Same_sources_flag = 0, N = 4)
```



Input format
-----------------------
The input to *FEAST* is composed of two tab-delimited ASCII text files:

(1) count table - An m by n count matrix, where m is the number samples and n is the number of taxa. Row names are the sample ids ('SampleID'). Column names are the taxa ids. Every consecutive column contains read counts for each sample. Note that this order must be respected.


count matrix (first 4 rows and columns):

| | ERR525698 |ERR525693 | ERR525688| ERR525699|
| ------------- | ------------- |------------- |------------- |------------- |
| taxa_1  |  0 | 5 | 0|20 |
| taxa_2  |  15 | 5 | 0|0 |
| taxa_3  |  0 | 13 | 200|0 |
| taxa_4  |  4 | 5 | 0|0 |



(2) metadata - An m by 3 table, where m is the number of samples. The metadata table has three columns (i.e., 'Env', 'SourceSink', 'id'). The first column is a description of the sampled environment (e.g., human gut), the second column indicates if this sample is a source or a sink (can take the value 'Source' or 'Sink'). The third column is the Sink-Source id. When using multiple sinks, each tested with the same group of sources, only the rows with 'SourceSink' = Sink will get an id (between 1 - number of sinks in the data). In this scenario, the sources’ ids are blank. When using multiple sinks, each tested with a distinct group of sources, each combination of sink and its corresponding sources should get the same id (between 1 - number of sinks in the data). Note that these names must be respected.


*using multiple sinks, each tested with the same group of sources:

| SampleID | Env |SourceSink | id |
| ------------- | ------------- |------------- |-------------|
| ERR525698  |  infant gut 1 | Sink | 1
| ERR525693  |  infant gut 2 | Sink | 2 |
| ERR525688   |  Adult gut 1 | Source| NA |
| ERR525699  |  Adult gut 2 | Source | NA |
| ERR525697  |  Adult gut 3 | Source | NA |


*using multiple sinks, each tested with a different group of sources:

| SampleID | Env |SourceSink | id |
| ------------- | ------------- |------------- |-------------|
| ERR525698  |  infant gut 1 | Sink | 1
| ERR525688   |  Adult gut 1 | Source| 1 |
| ERR525691  |  Adult gut 2 | Source | 1 |
| ERR525699  |  infant gut 2 | Sink | 2 |
| ERR525697  |  Adult gut 3 | Source | 2 |
| ERR525696  |  Adult gut 4 | Source | 2 |


 

Output - 

| infant gut 2  |Adult gut 1 | Adult gut 2| Adult gut 3| Adult skin 1 |  Adult skin 2|  Adult skin 3| Soil 1 | Soil 2 | unknown|
| ------------- | ------------- |------------- |------------- |------------- |------------- |------------- |------------- |------------- |------------- |
|  5.108461e-01  |  9.584116e-23 | 4.980321e-12 | 2.623358e-02|5.043635e-13 | 8.213667e-59| 1.773058e-10 |  2.704118e-14 |  3.460067e-02 |  4.283196e-01 |




