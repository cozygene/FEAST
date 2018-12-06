library("Rcpp")
library("vegan")
library("dplyr")
library("doParallel")
library("foreach")
library("mgcv")
library("reshape2")
library("ggplot2")
library("cowplot")
library("lsei")
library("RcppArmadillo")
cppFunction("arma::mat schur(arma::mat& a, arma::mat& b) 
            {return(a % b); }", depends="RcppArmadillo")

  #Home microbiome ----
  # rm(list = ls())
  # gc()
  

  start_time = 1
  end_time = 14
  first_sink = 1
  last_sink = 1
  rarefy_sink_flag = -1
  locale = "local" 
  COVERAGE= 10000
  BODYPART=T
  HOUSEONLY = T
  filter_otu = F
  init = "init" #"init1"
  
  house_number = 1 # 1m #7b #7_b #7mb
  dataset = paste0("house",house_number,"_",COVERAGE,"_",init,"_nofilter")

  
  source("SourceTracker.R")
  source("src.R")
  
  EM_ITERATIONS<-100
  include_epsilon <- T

  all_study_metadata <- read.table(paste0("../data/house",house_number,"_study_metadata.txt"),h=T,sep="\t",check=F,stringsAsFactors = FALSE)
  all_study_otus <- read.table(paste0("../data/house",house_number,"_study_otus.txt"),h=T,sep="\t",row.names=1,check=F,stringsAsFactors = FALSE)
  
  CI = list()
  
  for(t in start_time:end_time){
    savestr<-paste0('_house',house_number,'_t', t, '_')
    if(BODYPART){
      savestr<-paste0('_bodypart_house',house_number,'_t', t, '_')
    }
    
    
    metadata <- all_study_metadata[grep(paste0("_",house_number,"_",t,"$"),all_study_metadata$Env),]
    otus <- all_study_otus[grep(paste0("_",house_number,"_",t,"$"),all_study_metadata$Env),]
    row.names(metadata)
    idx<-apply(metadata, 1, function(var) all(var != ''))
    metadata<-metadata[idx,]
    print("new assign otu")
    

    ####  Get the data
    if(HOUSEONLY==T){
      houseidx<-as.numeric(max(grep(pattern='2192', row.names(metadata))))
      metadata<-metadata[1:houseidx,]
    }
    
    if(BODYPART==T){
      people<-as.character(metadata[grep('person', metadata$Env), "Env"])
      people<-unique(unlist(lapply(strsplit(people, split='_'), function(l) l[[1]])))
      for(pers in people){
        pidx<-grep(pers, metadata$Env)
        for(partix in pidx){
          toset<-strsplit(as.character(metadata[partix, "Description"]), split = '[.]')[[1]][2]
          toset<-paste0(paste0(strsplit(as.character(metadata[partix, "Env"]), split='_')[[1]][1:2], collapse='', sep='_'), toset)
          metadata[partix, 'Env']<-toset
        }
      }
      
      people<-as.character(metadata[grep('cat', metadata$Env), "Env"])
      people<-unique(unlist(lapply(strsplit(people, split='_'), function(l) l[[1]])))
      for(pers in people){
        pidx<-grep(pers, metadata$Env)
        for(partix in pidx){
          toset<-strsplit(as.character(metadata[partix, "Description"]), split = '[.]')[[1]][2]
          toset<-paste0(paste0(strsplit(as.character(metadata[partix, "Env"]), split='_')[[1]][1:2], collapse='', sep='_'), toset)
          metadata[partix, 'Env']<-toset
        }
      }
    }
    
    

    common_names<-intersect(rownames(metadata), rownames(otus))
    metadata<-metadata[common_names,]
    otus<-otus[common_names,]
    sinks<-rownames(metadata)[which(metadata$SourceSink == "sink")]
    sinks<-otus[sinks,]
    
    
    sources<-rownames(metadata)[which(metadata$SourceSink == "source")]
    envs<-as.character(metadata[sources, 2])
    
    sources<-otus[sources,]
    print("Sequencing depth")
    print(rowSums(sources))
    
    #######################
    
    
    st<-sourcetracker(sources, envs, COVERAGE)
    sources<-st$sources[-which(rownames(st$sources) == 'Unknown'),]
    
    envs <- row.names(sources)
    
    js_values<-c()
    x <- sources[c(1:dim(sources)[1]),]
    JSDMatrix <- jsdmatrix(x)
    JSDMatrix <- JSDMatrix/apply(sources, 1, sum)[1]
    JS = mean(JSDMatrix[-which(JSDMatrix == 0)])
    js_values = append(js_values, JS)
    print(js_values)
    my_sinks <- c("Kitchen.Counter","Bathroom.Door.Knob","Front.Door.Knob","Kitchen.Light.Switch","Bedroom.Floor","Kitchen.Floor")
    my_sources <- c(rep("Hand",4),rep("Foot",2))
    
    original_sources <- sources
    original_sinks <- sinks
    
    
    for(i in first_sink:last_sink){
      print(paste("time",i))
      
      print(row.names(sinks)[i])
      
      if(grepl("norarefy",dataset) == FALSE){
        sink<-rarefy(t(sinks[i,]), apply(sources, 1, sum)[1])
      }

      ####  Initialize the unknown source using raw sources and sink
      unk_source<-unknown__initialize_1(sources, sink, n_sources=dim(sources)[1])
      if(grepl("norarefy",dataset) == FALSE){
        unk_source<-rarefy(t(unk_source), apply(sources, 1, sum)[1])
      }

      em_source<-rbind(sources, unk_source)
      CI_em_source <- em_source
      labels<-c(levels(st$envs), "Unknown")
      print(labels)
      obs_source<-st$sources
      initalphs<-runif(dim(em_source)[1], min = 0.0, max = 1.0)
      initalphs <- initalphs / Reduce('+', initalphs)
      ##  Run EM
      em_source<-split(em_source, seq(nrow(em_source)))
      obs_source<-split(obs_source, seq(nrow(obs_source)))
      em_source<-lapply(em_source,t)
      obs_source<-lapply(obs_source,t)
      print("Start FEAST")
      print(Sys.time())
      
      tmp<-do_EM(alphas=initalphs, sources=em_source, 
                 observed=obs_source, sink = t(sink), iterations = EM_ITERATIONS)
      print(Sys.time())
      print("End FEAST")
      print(tmp$toret)
      #save(results$proportions, file=paste0('sourcetracker',savestr ))
      names(tmp$toret)<-labels
      pred_emnoise = tmp$toret
      

      
      

    }
  }
  
  
