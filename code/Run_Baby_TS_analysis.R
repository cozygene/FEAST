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

  
  #Baby time-series----
  # rm(list = ls())
  # gc()
  

  source("SourceTracker.R")
  source("src.R")
 
  unk=T; noise=T; bs=F; clsinit=F; em_itr=100; eps = T; include_epsilon = T

  # Load sample delivery mode metadata 
  Delivery_mode <- read.csv("../data/Backhed_metadata.csv")
  Delivery_mode <- Delivery_mode[,c(1,2,6)]
  
  
  Delivery_mode_id = c()
  for(j in 1:length(unique(Delivery_mode$patient_number))){
    
    Delivery_mode_id[j] = Delivery_mode$delivery_mode[Delivery_mode$patient_number %in% unique(Delivery_mode$patient_number)[j]][1]
  }
  
  Delivery_mode_id = as.character(Delivery_mode_id)
  
    
  # load sample metadata
  metadata <- read.table('../data/metadata_baby_ts.txt',sep='\t',h=T,check=F,comment='', row.names = 1)
  head(metadata)
  
  # load OTU table
  otus <- read.table('../data/otu_baby_ts.txt',sep='\t', header = T)
  otus <- t(as.matrix(otus))
  
  
  # extract only those samples in common between the two tables
  common.sample.ids <- intersect(rownames(metadata), rownames(otus))
  otus <- otus[common.sample.ids,]
  metadata <- metadata[common.sample.ids,]
  # double-check that the mapping file and otu table
  # had overlapping samples
  if(length(common.sample.ids) <= 1) {
    message <- paste(sprintf('Error: there are %d sample ids in common '),
                     'between the metadata file and data table')
    stop(message)
  }
  
  # extract the source environments and source/sink indices
  train.ix <- which(metadata$SourceSink=='Source')
  test.ix <- which(metadata$SourceSink=='Sink')
  if(is.element('Description',colnames(metadata))) desc <- metadata$Description
  
  
  alpha1 <- alpha2 <- 0.001
  cls_results_1 = c();
  cls_results_2 = c();
  st_results = c();
  
  em_rts<-c(); emnoise_rts =c();
  em_results<-c();
  emnoise_results<-c();
  
  st_runtimes <- c()
  cls_runtimes <- c()
  em_runtimes<-c()
  emnoise_runtimes<-c()
  
  envs <- metadata$Env
  num_sources <- 3
  
  sinks = otus[test.ix,]
  sink_samples = metadata[metadata$SourceSink=='Sink', 1]
  k = which(metadata$SourceSink == 'Sink' )
  id_uniq = unique(metadata$id)
  itr = length(k)
  

  CI = list()
  js_values = c()
  emnoise_results = c()
  em_results = c()
  js_values = c()

  
  labels <- c("4M", "Birth","Mother", "Unknown")
  
  for (it in 1:length(id_uniq)){
    
    metadata[k[it] ,1]
    
    sinks = otus[which(rownames(otus) == as.character(rownames(metadata)[which(metadata$Env == sink_samples[it])])),]
    train.ix <- which(metadata$SourceSink=='Source' & metadata$id == id_uniq[it])
    test.ix <- which(metadata$SourceSink=='Sink' & metadata$id == id_uniq[it])
    
    COVERAGE = min(apply(otus[c(train.ix,test.ix),], 1, sum))
    
    st <- sourcetracker(otus[train.ix,], envs[train.ix], rarefaction_depth = COVERAGE)
   
    print(paste("Number of OTUs in the sink sample = ",length(which(sinks > 0))))
    print(paste("Seq depth in the sink sample = ",sum(sinks)))
    print(paste("The sink is:", envs[test.ix]))
    
    source = st$sources[-dim(st$sources)[1],]
    totalsource<-source
    totalsource<-as.matrix(totalsource)
    sources <- split(totalsource, seq(nrow(totalsource)))
    sources<-lapply(sources, as.matrix)
    dists<-lapply(sources, function(x) x/(sum(colSums(x))))
    totaldist<-t(Reduce("cbind", dists))
    sinks<-matrix(sinks, nrow = 1, ncol = dim(totalsource)[2])
    
    num_sources = dim(source)[1]
    envs_simulation = c(1:(num_sources))
    
    source_old = source
    totalsource_old = totalsource
    
    source_old=lapply(source_old,t)
    source_old<- split(totalsource_old, seq(nrow(totalsource_old)))
    source_old<-lapply(source_old, as.matrix)
    
    #Creating the unknown source initial per mixing iteration
    if(include_epsilon == TRUE){
      
      ##Adding the initial value of the unknown source for FEAST
      source_2 = list()
      totalsource_2 = matrix(NA, ncol = dim(totalsource_old)[2], nrow = ( dim(totalsource_old)[1] + 1))
      
      for(j in 1:num_sources){
        
        source_2[[j]] = source_old[[j]]
        totalsource_2[j,] = totalsource_old[j,]
      }
      
      #create unknown for each sink i
      sinks_rarefy = rarefy(matrix(sinks, nrow = 1), maxdepth = apply(totalsource_old, 1, sum)[1]) #make
      
      unknown_source_1 = unknown_initialize(sources = totalsource[c(1:num_sources),], sink = as.numeric(sinks),
                                            n_sources = num_sources)
      
      unknown_source = unknown_source_1 + rpois(n = length(sinks), lambda = 1)
      
      unknown_source_rarefy = rarefy(matrix(unknown_source, nrow = 1), maxdepth = COVERAGE)
      source_2[[j+1]] = t(unknown_source_rarefy)
      totalsource_2[(j+1),] = t(unknown_source_rarefy)
      totalsource = totalsource_2
      
      source=lapply(source_2,t)
      source<- split(totalsource, seq(nrow(totalsource_2)))
      source<-lapply(source_2, as.matrix)
      envs_simulation <- c(1:(num_sources+1))
      
    }
    
    
    samps <- source
    samps<-lapply(samps, t)
    observed_samps <- samps
    observed_samps[[(num_sources + 1)]] = t(rep(0, dim(samps[[1]])[2]))

    
    #Calculate JSD value
    x <- totalsource[c(1:num_sources),]
    JSDMatrix <- jsdmatrix(x)
    JSDMatrix <- JSDMatrix/COVERAGE
    JS = mean(JSDMatrix[-which(JSDMatrix == 0)])
    js_values = append(js_values, JS)
    print(js_values)
    
    initalphs<-runif(num_sources+1, 0.0, 1.0)
    initalphs=initalphs/Reduce("+", initalphs)
    sink_em = as.matrix(sinks)
    em_start_time<-as.numeric(proc.time()[3])
    pred_em<-do_EM_basic(alphas=initalphs, sources=samps, sink=sink_em, iterations=em_itr)
    
    
    em_rts=c(em_rts, as.numeric(proc.time()[3])-em_start_time)
    emnoise_start_time<-as.numeric(proc.time()[3])
    tmp<-do_EM(alphas=initalphs, sources=samps, sink=sink_em, iterations=em_itr, observed=observed_samps)
    pred_emnoise = tmp$toret
    
    if(num_sources < 3){
      
      pred_emnoise = c(pred_emnoise[1:2], NA, pred_emnoise[3])
      pred_em = c(pred_em[1:2], NA, pred_em[3])
    }
    
    
    
    emnoise_rts=c(emnoise_rts, as.numeric(proc.time()[3])-emnoise_start_time)
    emnoise_results=c(emnoise_results, pred_emnoise)
    em_results<-c(em_results, pred_em)
    
    names(pred_emnoise) = labels
    names(pred_em) = labels
    
    sort(pred_em, decreasing = T)
    sort(pred_emnoise, decreasing = T)
    
    
    print(paste0("FEAST"))
    print(paste0(round(pred_em*100, 3), "%", sep = ""))
    
  }
  
  
  baby_prop_em = matrix(em_results, ncol = 4, nrow = 98, byrow = T)
  rownames(baby_prop_em) = unique(metadata$id)
  colnames(baby_prop_em) = c("Birth", "4Months", "Mother", "Unknown")

  
  print(t.test(baby_prop_em[,3][Delivery_mode_id == 2], baby_prop_em[,3][Delivery_mode_id == 1]))
  print(wilcox.test(baby_prop_em[,3][Delivery_mode_id == 2], baby_prop_em[,3][Delivery_mode_id == 1]))

  
  
  toplot<-data.frame(baby_prop_em)
  colnames(toplot) = c("Birth", "4Months", "Mother", "Unknown")
  toplot$Delivery_mode_id <- Delivery_mode_id
  
  toplot_Delivery_1 <- toplot[toplot$Delivery_mode_id == 1,]
  toplot_Delivery_1.melt <- melt(toplot_Delivery_1, variable.name= "Origin",
                                 value.name = "Proportion")
  
  toplot_Delivery_1_plot <- ggplot(toplot_Delivery_1.melt, aes(x= Origin, y = Proportion, col= Origin)) + 
    geom_boxplot() + theme_bw() + ggtitle("C-section")
  
  
  toplot_Delivery_2 <- toplot[toplot$Delivery_mode_id == 2,]
  toplot_Delivery_2.melt <- melt(toplot_Delivery_2, variable.name= "Origin",
                                 value.name = "Proportion")
  
  toplot_Delivery_2_plot <- ggplot(toplot_Delivery_2.melt, aes(x= Origin, y = Proportion, col= Origin)) +
    geom_boxplot() + theme_bw() + ggtitle("Vaginal delivery")
  
  Baby_ts_plot = plot_grid(toplot_Delivery_1_plot, toplot_Delivery_2_plot,  nrow = 1, align = 'h')
  
  ggsave(filename="../results/Baby_time_series.png", plot = Baby_ts_plot , dpi = 600, width = 8.75, height = 6.1, units = "in")
