rm(list = ls())
gc()

print("Change directory path")
dir_path = paste("~/FEAST_src/")
setwd(paste0(dir_path, "code"))
source("src.R")

#Set the arguments of your data
metadata_file = "your_metadata_file_name"
count_matrix = "your_taxa_count_matrix"
num_sources <- number_of_sources 
EM_iterations = 1000 #default value

setwd(paste0(dir_path, "Data_files"))

# Load sample metadata
metadata <- read.csv(file = metadata_file ,header = T, sep = "\t", row.names = 1)

# Load OTU table
otus <- read.table(file = count_matrix, header = T, comment = '', check = F, sep = '\t')
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

# Extract the source environments and source/sink indices

train.ix <- which(metadata$SourceSink=='Source')
test.ix <- which(metadata$SourceSink=='Sink')
COVERAGE =  min(rowSums(otus[c(train.ix, test.ix),]))  #Can be adjusted by the user

# Define sources and sinks

sources <- as.matrix(rarefy(otus[train.ix,], COVERAGE))
sinks <- as.matrix(rarefy(t(as.matrix(otus[test.ix,])), COVERAGE))


print(paste("Number of OTUs in the sink sample = ",length(which(sinks > 0))))
print(paste("Seq depth in the sources and sink samples = ",COVERAGE))
print(paste("The sink is:", envs[test.ix]))

# Estimate source proportions for each sink

FEAST_output<-FEAST(source=sources, sinks = t(sinks), env = envs[train.ix], em_itr = EM_iterations, COVERAGE = COVERAGE)
Proportions_est <- FEAST_output$data_prop[,1]
names(Proportions_est) <- c(as.character(envs[train.ix]), "unknown")

print("Source mixing proportions")
Proportions_est

