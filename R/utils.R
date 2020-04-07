

#' Load_metadata
#'
#' Load metadata file - A tab-separated ASCII text file.
#' The row names of  metadata table are the unique sample ids.
#' The metadata file has three colunms (i.e., 'Env', 'SourceSink', 'id').
#' The first column is a description of the sampled environment (e.g., human gut), the second column indicates
#' if this sample is a source or a sink (can take the value 'Source' or 'Sink'). The fourth column is the Sink-Source id.
#' When using multiple sinks, each tested with the same group of sources, only the rows with 'SourceSink' = Sink will get an id
#' (between 1 - number of sinks in the data). In this scenatio, the sources ids are blank. When using multiple sinks,
#' each tested with a distinct group of sources, each combination of sink and its corresponding sources should get the
#' same id (between 1 - number of sinks in the data). Note that these names must be respected.
#'
#' @param metadata_path A path to a tab-separated ASCII text metadata file.
#' @export
Load_metadata <- function(metadata_path){

  metadata <- read.csv(metadata_path,h=T, sep = '\t', row.names = 1,check=FALSE)
  return(metadata)
}


#' Load count matrix
#'
#' Load counts file - A tab-separated ASCII text file.
#' A matrix \eqn{C} of samples by taxa with the sources and sinks.
#' The row names are the unique sample ids. The column names are the unique taxa ids.
#' \eqn{Cij} is the read counts of taxa j in sample i each sample. Note that this order must be respected
#'
#' @param CountMatrix_path A path to a tab-separated ASCII text counts file
#' @export
Load_CountMatrix <- function(CountMatrix_path){

  CountMatrix <- read.table(CountMatrix_path, header = T, comment = '', check = F, sep = '\t', row.names = 1)
  CountMatrix <- t(as.matrix(CountMatrix))
  return(CountMatrix)
}

#' Perform rarefaction analysis
#' @export
FEAST_rarefy <- function(x,maxdepth){


  if(is.null(maxdepth)) return(x)

  if(!is.element(class(x), c('matrix', 'data.frame','array')))
    x <- matrix(x,nrow=nrow(x))
  nr <- nrow(x)
  nc <- ncol(x)

  for(i in 1:nrow(x)){
    if(sum(x[i,]) > maxdepth){
      prev.warn <- options()$warn
      options(warn=-1)
      s <- sample(nc, size=maxdepth, prob=x[i,], replace=T)
      options(warn=prev.warn)
      x[i,] <- hist(s,breaks=seq(.5,nc+.5,1), plot=FALSE)$counts
    }
  }
  return(x)
}

#' Calculate Kullback–Leibler divergence
#' @export
kld <- function(p,q){
  nonzero <- p>0 & q>0
  return(sum(p[nonzero] * log2(p[nonzero]/q[nonzero])))
}

#' Calculate Jensen–Shannon divergence
#' @export
jsd <- function(p,q){
  m <- (p + q)/2
  return((kld(p,m) + kld(q,m))/2)
}

#' Calculate pairwise Jensen–Shannon divergence
#' @export
jsdmatrix <- function(x){
  d <- matrix(0,nrow=nrow(x),ncol=nrow(x))
  for(i in 1:(nrow(x)-1)){
    for(j in (i+1):nrow(x)){
      d[i,j] <- jsd(x[i,], x[j,])
      d[j,i] <- d[i,j]
    }
  }
  return(d)
}

