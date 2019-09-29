
# FEAST QIIME2 Plugin


One critical challenge in analyzing microbiome communities is due to their composition; each of them is typically comprised of several source environments, including different contaminants as well as other microbial communities that interacted with the sampled habitat. To account for this structure, we developed FEAST (Fast Expectation-mAximization microbial Source Tracking), a ready-to-use scalable framework that can simultaneously estimate the contribution of thousands of potential source environments in a timely manner, thereby helping unravel the origins of complex microbial communities. Specifically, FEAST is quantifying the fraction, or proportion, of different microbial samples (sources) in a target microbial community (sink), by leveraging its structure and measuring the respective similarities between a sink community and potential source environments. For more details see Shenhav et al., Nature Methods 2019 (https://www.nature.com/articles/s41592-019-0431-x).

## Install q2-FEAST

If you have not already done so, activate your QIIME environment.

```shell
source activate qiime2-20xx.x
```
Next we will need to ensure some dependancies are installed.

```shell
conda install -c bioconda -c conda-forge -c r bioconductor-phyloseq r-devtools r-magrittr r-dplyr r-vgam r-tidyr r-vegan r-reshape2 r-rcpp r-rcpparmadillo r-gridextra r-ggplot2 r-ggthemes   
```

Now we will install FEAST and the q2-plugin.

```R
# the main FEAST package
> R
> devtools::install_github("cozygene/FEAST")
> quit()
```
```shell
# the QIIME2 plugin
pip install git+https://github.com/cozygene/FEAST.git
```

# Tutorial 

A QIIME2 tutorial is available [here](https://github.com/cozygene/FEAST/q2_feast/tutorials/DIABIMMUNE.md)

## Running q2-FEAST microbial tracking

```shell
$ qiime feast microbialtracking --help
Usage: qiime feast microbialtracking [OPTIONS]

  A major challenge of analyzing the compositional structure of microbiome
  data is identifying its potential origins. Here, we introduce fast
  expectation-maximization microbial source tracking (FEAST), a ready-to-use
  scalable framework that can simultaneously estimate the contribution of
  thousands of potential source environments in a timely manner, thereby
  helping unravel the origins of complex microbial communities. The
  information gained from FEAST may provide insight into quantifying
  contamination, tracking the formation of developing microbial communities,
  as well as distinguishing and characterizing bacteria-related health
  conditions.

Inputs:
  --i-table ARTIFACT FeatureTable[Frequency]
                       Feature table file containing sources and sinks for
                       source tracking.                             [required]
Parameters:
  --m-metadata-file METADATA...
    (multiple          Sample metadata file containing sources and sinks for
     arguments will    source tracking.
     be merged)                                                     [required]
  --p-environment-column TEXT
                       Sample metadata column with a description of the
                       sampled environment (e.g., human gut).       [required]
  --p-source-sink-column TEXT
                       Sample metadata column with labels for source or a
                       sink.All the sub-classes in this column must be in
                       either source-ids or sink-ids.               [required]
  --p-source-ids TEXT  Comma-separated list (without spaces) of class ids
                       contained in source-sink-column to be considered as
                       sources.                                     [required]
  --p-sink-ids TEXT    Comma-separated list (without spaces) of class ids
                       contained in source-sink-column to be considered as
                       sinks.                                       [required]
  --p-shared-id-column TEXT
                       Sample metadata column with the Sink-Source id. When
                       using multiple sinks, each tested with the same group
                       of sources                                   [required]
  --p-em-iterations INTEGER
                       A numeric value indicating the number of EM
                       iterations.                             [default: 1000]
  --p-different-sources / --p-no-different-sources
                       A Boolean value indicating the source-sink
                       assignment.Different-sources is True if different
                       sources are assignedto each sink, otherwise
                       different-sources should be False.      [default: True]
Outputs:
  --o-mixing-proportions ARTIFACT FeatureTable[RelativeFrequency]
                       The mixing proportions returned from FEAST. The mixing
                       proportions table is an S1 by S2 matrix P, where S1 is
                       the number sinks and S2 is the number of sources
                       (including an unknown source). Each row in matrix P
                       sums to 1. Pij is the contribution of source j to sink
                       i. If Pij == NA it indicates that source j was not used
                       in the analysis of sink i.                   [required]
Miscellaneous:
  --output-dir PATH    Output unspecified results to a directory
  --verbose / --quiet  Display verbose output to stdout and/or stderr during
                       execution of this action. Or silence output if
                       execution is successful (silence is golden).
  --citations          Show citations and exit.
  --help               Show this message and exit.
```
