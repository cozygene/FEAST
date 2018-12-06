library("Rcpp")
library("vegan")
library("dplyr")
library("doParallel")
library("foreach")
library("mgcv")
library("reshape2")
library("ggplot2")
library("cowplot")
library("RcppArmadillo")
cppFunction("arma::mat schur(arma::mat& a, arma::mat& b) 
            {return(a % b); }", depends="RcppArmadillo")

source("SourceTracker.R")
source("src.R")


dat<-read.csv("/data/dysbiosis_otu_table.txt", sep='\t',h=T,check=F,comment='', row.names = 1)
dat<-dat[which(rowSums(dat) > 8000),]
sink_mapping<-read.table("/data/sinks_true_mapping.txt", sep='\t', header=T, 
                         stringsAsFactors = F, comment.char = '')
source_mapping<-read.csv("/data/complete_mapping_dysbiosis.csv", sep=',', header = T,
                         stringsAsFactors = F)
source_mapping<-source_mapping[-which(source_mapping$SourceSink=='sink'),]
sinks_included <- sink_mapping$X.SampleID[which(sink_mapping$X.SampleID %in% rownames(dat))]
sinks<-dat[sinks_included,]



load(file = "/data/FEAST_healthy.RData")
load(file = "/data/FEAST_icu.RData")


#Healthy sources
sources_icu_mapping <- em_noise_results_icu$sources_names
sources_healthy_mapping <- em_noise_healthy$sources_names


Healthy_sources_icu <- sinks[which(rownames(sinks) %in% sources_icu_mapping),]
Healthy_sources_healthy <- sinks[which(rownames(sinks) %in% sources_healthy_mapping),]

Healthy_sources_icu <- as.matrix(Healthy_sources_icu)
Healthy_sources_healthy <- as.matrix(Healthy_sources_healthy)

sinks_icu_mapping <- em_noise_results_icu$sink_icu
sinks_healthy_mapping <- em_noise_healthy$sink_healthy

sinks_icu <- sinks[which(rownames(sinks) %in% sinks_icu_mapping),]
sinks_healthy <- sinks[which(rownames(sinks) %in% sinks_healthy_mapping),]

sinks_icu <- as.matrix(sinks_icu)
sinks_healthy <- as.matrix(sinks_healthy)
COVERAGE = min(rowSums(Healthy_sources_icu))



#sink = ICU ----
st_flag = 0
em_flag = 1
EM_ITERATIONS = 100
healthy_flag = 0
icu_flag = 1
N = dim(Healthy_sources_icu)[1]
Num_ind = 10
start_ind = 1

em_results_icu_feces_new  = matrix(NA, ncol = 2, nrow = N )
em_noise_results_icu_feces_new  = matrix(NA, ncol = 2, nrow = N )
ST_results_icu_feces_new = matrix(NA, ncol = 2, nrow = N )


i = 1
for(samp in sinks_icu_mapping[start_ind:Num_ind]){
  
  print(samp)
  
  if(length(which(rownames(sinks_icu) %in% samp)) > 0){

  
  sink<-as.matrix(rarefy(as.matrix(sinks_icu[which(rownames(sinks_icu) %in% samp),]), apply(Healthy_sources_icu, 1, sum)[1]))
  
  sources_data = Healthy_sources_icu
  sink_data = sink
  
  if(st_flag == 1){
    
    env = rep("AGP-data",N)
    
    st_type <- sourcetracker(sources_data,  env, rarefaction_depth = COVERAGE)
    results <- predict(st_type,sink, alpha1=alpha1, alpha2=alpha2) #change alpha
    
    if(healthy_flag == 1)
      ST_results_healthy[i,] = results$proportions
    
    
    if(icu_flag == 1)
      ST_results_icu_feces_new[i,] = results$proportions
  }
  
  if(em_flag == 1){
    
    em_res<-EM_results(source=sources_data, sinks = sink, env = rownames(sources_data), 
                       em_itr = EM_ITERATIONS, include_epsilon = T, COVERAGE = COVERAGE)
    print(em_res$data_prop)
    
    
    em_res$data_prop$Type = c(rep("AGP-data",N) , "Unknown")
    
    Type_uniq = unique(em_res$data_prop$Type)
    
    
    if(healthy_flag == 1){
      
      
      for(k in 1:length(Type_uniq)){
        
        em_results_healthy[i,k] = sum(em_res$data_prop$pred_em_all[em_res$data_prop$Type == Type_uniq[k]])
        em_noise_results_healthy[i,k] = sum(em_res$data_prop$pred_emnoise_all[em_res$data_prop$Type == Type_uniq[k]])
        
        
      }
      
      print(Type_uniq)
      print(em_results_healthy[i,])
      print(em_noise_results_healthy[i,])
      
      
    }
    
    if(icu_flag == 1){
      
      # em_results_icu_feces_new = c()
      # em_noise_results_icu_feces_new = c()
      
      
      for(k in 1:length(Type_uniq)){
        
        em_results_icu_feces_new[i,k] = sum(em_res$data_prop$pred_em_all[em_res$data_prop$Type == Type_uniq[k]])
        em_noise_results_icu_feces_new[i,k] = sum(em_res$data_prop$pred_emnoise_all[em_res$data_prop$Type == Type_uniq[k]])
        
        
      }
      
      
      print(Type_uniq)
      # print(em_results_icu_feces_new[i,])
      print("ICU")
      print(em_noise_results_icu_feces_new[i,])
      i = i+1
      
    }
   }
  }
}


st_flag = 0
em_flag = 1
EM_ITERATIONS = 100
healthy_flag = 1
icu_flag = 0
N = dim(Healthy_sources_healthy)[1]
Num_ind = 10
start_ind = 1


em_results_healthy = matrix(NA, ncol = 2, nrow = N )
em_noise_results_healthy = matrix(NA, ncol = 2, nrow = N )
ST_results_healthy = matrix(NA, ncol = 2, nrow = N )

i = 1
#sink = healthy ----
for(samp in sinks_healthy_mapping[start_ind:Num_ind]){
  
  print(samp)
  
  if(length(which(rownames(sinks_healthy) %in% samp)) > 0){
  
  sink<-as.matrix(rarefy(as.matrix(sinks_healthy[which(rownames(sinks_healthy) %in% samp),]), apply(Healthy_sources_healthy, 1, sum)[1]))
  
  sources_data = Healthy_sources_healthy
  sink_data = sink
  
  if(st_flag == 1){
    
    env = rep("AGP-data",N)
    
    st_type <- sourcetracker(sources_data,  env, rarefaction_depth = COVERAGE)
    results <- predict(st_type,sink, alpha1=alpha1, alpha2=alpha2) #change alpha
    
    if(healthy_flag == 1)
      ST_results_healthy[i,] = results$proportions
    
    
    if(icu_flag == 1)
      ST_results_icu_feces_new[i,] = results$proportions
  }
  
  if(em_flag == 1){
    
    em_res<-EM_results(source=sources_data, sinks = sink, env = rownames(sources_data), 
                       em_itr = EM_ITERATIONS, include_epsilon = T, COVERAGE = COVERAGE)
    print(em_res$data_prop)
    
    
    em_res$data_prop$Type = c(rep("AGP-data",N) , "Unknown")
    
    Type_uniq = unique(em_res$data_prop$Type)
    
    
    if(healthy_flag == 1){
      
      
      for(k in 1:length(Type_uniq)){
        
        em_results_healthy[i,k] = sum(em_res$data_prop$pred_em_all[em_res$data_prop$Type == Type_uniq[k]])
        em_noise_results_healthy[i,k] = sum(em_res$data_prop$pred_emnoise_all[em_res$data_prop$Type == Type_uniq[k]])
        
        
      }
      
      
      print(Type_uniq)
      # print(em_results_healthy[i,])
      print("healthy")
      print(em_noise_results_healthy[i,])
      i = i+1
      
    }
    
    if(icu_flag == 1){
      
      
      for(k in 1:length(Type_uniq)){
        
        em_results_icu_feces_new[i,k] = sum(em_res$data_prop$pred_em_all[em_res$data_prop$Type == Type_uniq[k]])
        em_noise_results_icu_feces_new[i,k] = sum(em_res$data_prop$pred_emnoise_all[em_res$data_prop$Type == Type_uniq[k]])
        
        
      }
      
      print(Type_uniq)
      print(em_results_icu_feces_new[i,])
      print(em_noise_results_icu_feces_new[i,])
      
      
    }
  }
 }
  
}


FEAST_icu = list(em_noise_results_icu_feces = em_noise_results_icu_feces_new,
                   sink_icu = sinks_icu_mapping,
                   sources_names = rownames(sources_data))
save(FEAST_icu, file = "/data/FEAST_ICU.RData")
  
FEAST_healthy = list(em_noise_results_healthy = em_results_healthy,
                      sink_healthy = sinks_healthy_mapping,
                      sources_names = rownames(sources_data))
save(FEAST_healthy, file = "/data/FEAST__Healthy.RData")
