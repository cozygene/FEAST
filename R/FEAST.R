#' A scalable algorithm for quantifying the origins of complex microbial communities
#'
#' FEAST performs microbial source tracking. As input, it takes
#' a count matrix \eqn{C}  of size \eqn{m} by \eqn{n} and a metadata table, of size \eqn{m} by \eqn{3},
#' where \eqn{m} is the number samples and \eqn{n} is the number of taxa.
#' The rownames of both \eqn{C} and the metadata table are the unique sample ids.


#' @param C An \eqn{m} by \eqn{n} count matrix, where \eqn{m} is the number samples and \eqn{n} is the number of taxa.
#' @param metadata An \eqn{m} by \eqn{3} table, where \eqn{m} is the number samples.
#' The metadata table has three colunms (i.e., 'Env', 'SourceSink', 'id').
#' The first column is a description of the sampled environment (e.g., human gut), the second column indicates
#' if this sample is a source or a sink (can take the value 'Source' or 'Sink'). The fourth column is the Sink-Source id.
#' When using multiple sinks, each tested with the same group of sources, only the rows with 'SourceSink' = Sink will get an id
#' (between 1 - number of sinks in the data). In this scenatio, the sources ids are blank. When using multiple sinks,
#' each tested with a distinct group of sources, each combination of sink and its corresponding sources should get the
#' same id (between 1 - number of sinks in the data). Note that these names must be respected.
#' @param EM_iterations A numeric value indicating the number of EM iterations (default 1000).
#' @param COVERAGE A numeric value indicating the rarefaction depth (default = minimal sequencing depth within each group of sink
#' and its corresponding sources).
#' @param different_sources_flag A boolian value indicating the source-sink assignment.
#' @param dir_path A path to an output .txt file.
#' @param outfile the prefix for saving the output file.
#' different_sources_flag = 1 if different sources are assigned to each sink , otherwise = 0.
#' @return P - an \eqn{S1} by \eqn{S2} matrix, where \eqn{S1} is the number sinks and \eqn{S2}
#' is the number of sources (including an unknown source). Each row in matrix \eqn{P} sums to 1.
#' \eqn{Pij} is the contribution of source j to sink i.
#' If \eqn{Pij} == NA it indicateds that source j was not used in the analysis of sink i.
#'
#' @examples
#' \donttest{
#'
#' #Load metadata table
#'
#' #Load count matrix
#'
#' #Calculate the sources' contributions to
#' #each sink sample in the data
#'FEAST_output <- FEAST(count_matrix = otus, metadata = metadata,
#'                      different_sources_flag = 1)
#' }
#'
#' @export
FEAST <- function (C, metadata, EM_iterations = 1000, COVERAGE = NULL, 
                   different_sources_flag, dir_path, outfile) 
{
  
  # C = otus;
  # metadata = metadata_sub;
  # different_sources_flag = 1;
  # dir_path = "~/FEAST_OLD_04072020/Test_results";
  # outfile="demo"
  
  if (sum(colnames(metadata) == "SampleID") == 0) 
    stop("The metadata file must contain a 'SampleID' column matching samples' names in the count table.")
  if (sum(colnames(metadata) == "Env") == 0) 
    stop("The metadata file must contain an 'Env' column naming the source environment for each sample.")
  if (sum(colnames(metadata) == "SourceSink") == 0) 
    stop("The metadata file must contain a 'SourceSink' column indicating 'source' or 'sink' for each sample.")
  if (sum(colnames(metadata) == "id") == 0) 
    stop("The metadata file must contain an 'id' column matching the source environments for each sink sample.")
  sink_ids <- grep("sink", as.character(metadata$SourceSink), 
                   ignore.case = TRUE, value = FALSE)
  source_ids <- grep("source", as.character(metadata$SourceSink), 
                     ignore.case = TRUE, value = FALSE)
  metadata$SourceSink[sink_ids] = "Sink"
  metadata$SourceSink[source_ids] = "Source"
  if (sum(as.integer(C) != as.numeric(C)) > 0) {
    stop("Data must be integral. Consider using \"ceiling(datatable)\" or ceiling(1000*datatable) to convert floating-point data to integers.")
  }
  common.sample.ids <- intersect(unique(metadata$SampleID), rownames(C))
  C <- C[common.sample.ids, ]
  if(length(unique(metadata$SampleID)) == length(metadata$SampleID))
    metadata <- metadata[common.sample.ids, ]
  if (length(common.sample.ids) <= 1) {
    message <- paste(sprintf("Error: there are %d sample ids in common "), 
                     "between the metadata file and data table")
    stop(message)
  }
  
  all_sinks_sampleID <- unique(metadata$SampleID[which(metadata$SourceSink == "Sink")])
  all_sources_sampleID <- unique(metadata$SampleID[which(metadata$SourceSink == "Source")])
  num_sinks <- length(all_sinks_sampleID)
  num_sources <- length(all_sources_sampleID)
  Ids <- na.omit(unique(metadata$id))
  envs <- paste0(metadata$Env, "_", metadata$SampleID)
  envs_sources <- paste0(all_sources_sampleID, "_", unique(metadata$Env[metadata$SampleID %in% all_sources_sampleID]))
  envs_sink <- paste0(all_sinks_sampleID, "_", unique(metadata$Env[metadata$SampleID %in%  all_sinks_sampleID]))
  proportions_mat <- matrix(NA, ncol = (num_sources + 1), nrow = num_sinks)
  idx.unknown <- num_sources + 1
  for (it in 1:num_sinks) {
    if (it%%10 == 0 || it == num_sinks) 
      print(paste0("Calculating mixinig proportions for sink ", 
                   it))
    if (different_sources_flag == 1) {
      # train.ix <- which(metadata$SourceSink == "Source" & 
      #                     metadata$id == Ids[it])
      # test.ix <- which(metadata$SourceSink == "Sink" & 
      #                    metadata$id == Ids[it])
      
      tmp_train <- metadata$SampleID[which(metadata$SourceSink == "Source" &
                                             metadata$id == Ids[it])]
      train.ix <- which(rownames(C) %in% tmp_train)
      tmp_test <- metadata$SampleID[which(metadata$SourceSink == "Sink" &
                                            metadata$id == Ids[it])]
      test.ix <- which(rownames(C) %in% tmp_test)
      
    }
    else {
      train.ix <- which(metadata$SourceSink == "Source")
      test.ix <- which(metadata$SourceSink == "Sink" & 
                         metadata$id == Ids[it])
    }
    num_sources <- length(train.ix)
    if (is.null(COVERAGE)) 
      COVERAGE <- min(rowSums(C[c(train.ix, test.ix), ]))
    if (COVERAGE > 0) {
      if (length(train.ix) == 1) 
        sources <- as.matrix(FEAST_rarefy(t(as.matrix(C[train.ix, 
                                                        ])), COVERAGE))
      if (length(train.ix) > 1) 
        sources <- as.matrix(FEAST_rarefy(C[train.ix, 
                                            ], COVERAGE))
      sinks <- as.matrix(FEAST_rarefy(t(as.matrix(C[test.ix, 
                                                    ])), COVERAGE))
      FEAST_output <- Infer.SourceContribution(source = sources, 
                                               sinks = t(sinks), env = envs[train.ix], em_itr = EM_iterations, 
                                               COVERAGE = COVERAGE)
      idx.sources <- which(all_sources_sampleID %in% rownames(sources))
      proportions_mat[it, c(idx.sources, idx.unknown)] <- FEAST_output$data_prop[,1]
    }
    else {
      cat(sprintf("Error: the sequencing depth in one sink or source sample in Id %d is zero. No proportions were generated\n", 
                  it))
    }
  }
  proportions_mat <- data.frame(proportions_mat)
  colnames(proportions_mat) <- c(envs_sources, "Unknown")
  rownames(proportions_mat) <- envs_sink
  setwd(dir_path)
  write.table(proportions_mat, file = paste0(outfile, "_source_contributions_matrix.txt"), 
              sep = "\t")
  return(proportions_mat)
}
