rm(list = ls())
gc()

print("Change dir_path")
dir_path = paste("~/Dropbox/FEAST-master/")
setwd(paste0(dir_path, "code"))
source("src.R")

setwd(paste0(dir_path, "Data_files"))

# Load sample metadata
metadata <- read.csv('metadata_example.txt',h=T, sep = "\t", row.names = 1)

# Load OTU table
otus <- read.table('otu_example.txt', header = T, comment = '', check = F, sep = '\t')
otus <- t(as.matrix(otus))


# Extract only those samples in common between the two tables
common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
# Double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}



EM_iterations = 100
envs <- metadata$Env
Ids <- unique(metadata$id)
Proportions_est <- list()
num_sources <- 3

for(it in 1:length(Ids)){
  

  # Extract the source environments and source/sink indices
  
  train.ix <- which(metadata$SourceSink=='Source' & metadata$id == Ids[it])
  test.ix <- which(metadata$SourceSink=='Sink' & metadata$id == Ids[it])
  COVERAGE =  min(rowSums(otus[c(train.ix, test.ix),]))  #Can be adjusted by the user
  
  # Define sources and sinks
  
  sources <- as.matrix(rarefy(otus[train.ix,], COVERAGE))
  sinks <- as.matrix(rarefy(t(as.matrix(otus[test.ix,])), COVERAGE))
  
  
  print(paste("Number of OTUs in the sink sample = ",length(which(sinks > 0))))
  print(paste("Seq depth in the sources and sink samples = ",COVERAGE))
  print(paste("The sink is:", envs[test.ix]))
  
  # Estimate source proportions for each sink
  
  FEAST<-EM_results(source=sources, sinks = t(sinks), env = envs[train.ix], em_itr = EM_iterations, COVERAGE = COVERAGE)
  Proportions_est[[it]] <- FEAST$data_prop[,1]
  
  
  names(Proportions_est[[it]]) <- c(as.character(envs[train.ix]), "unknown")
  
  if(length(Proportions_est[[it]]) < num_sources +1){

    tmp = Proportions_est[[it]]
    Proportions_est[[it]][num_sources] = NA
    Proportions_est[[it]][num_sources+1] = tmp[num_sources]
  }
  
  print("Source mixing proportions")
  print(Proportions_est[[it]])
  
  

}


Proportions_est_data = matrix(unlist(Proportions_est), ncol = (num_sources+1), byrow = T)
colnames(Proportions_est_data) = c("Baby_B", "Baby_4M", "Baby_M", "unknown")

