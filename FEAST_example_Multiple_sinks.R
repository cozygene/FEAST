rm(list = ls())
gc()

print("Change directory path")
dir_path = paste("~/FEAST/")
setwd(paste0(dir_path, "FEAST_src"))
source("src.R")

#Set the arguments of your data
metadata_file = "metadata_example_multi.txt"
count_matrix = "otu_example_multi.txt"
num_sources <- 3 
EM_iterations = 1000 #default value

setwd(paste0(dir_path, "Data_files"))
# Load sample metadata
metadata <- read.csv(metadata_file,h=T, sep = "\t", row.names = 1)

# Load OTU table
otus <- read.table(count_matrix, header = T, comment = '', check = F, sep = '\t')
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




envs <- metadata$Env
Ids <- unique(metadata$id)
Proportions_est <- list()

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
  
  FEAST_output<-FEAST(source=sources, sinks = t(sinks), env = envs[train.ix], em_itr = EM_iterations, COVERAGE = COVERAGE)
  Proportions_est[[it]] <- FEAST_output$data_prop[,1]
  
  
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
colnames(Proportions_est_data) = c("Env_2", "Env_3", "Env_4", "unknown")