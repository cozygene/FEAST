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
            

  
  #Cancer allo-HSCT----- 
  # rm(list = ls())
  # gc()
  
  source("src.R")
  
  MSK_flag = 1
  if(MSK_flag == 1){
    
    metadata <- read.csv('../data/events.tsv',h=T, sep = "\t")
    otus <- read.table('../data/otuTable.tsv', header = F, comment = '', check = F)
    
    day = as.numeric(otus[2,-1])
    patiendId = as.numeric(otus[1,-1])
    
  }
  
  taxa_summary = as.character(otus[-c(1:2),1])
  metadata_2 = data.frame(patiendId, day)
  
  otus = otus[-c(1,2),-1]
  colnames(otus) = paste0(patiendId,"_", day)
  otus = t(otus)
 
  EM_ITERATIONS = 500
  
  day_after = c()
  ind_after = c()
  
  day_before = c()
  ind_before = c()
  
  Results_array_emnoise_BD = list()
  Results_array_em_BD = list()
  Results_array_st_BD = list()
  em_results_BD = list()
  
  Results_array_emnoise_AD = list()
  Results_array_em_AD = list()
  Results_array_st_AD = list()
  em_results_AD = list()
  
  k = 1
  for(l in 1:length(unique(patiendId))){
    
    print(l)
    
    
    tmp = metadata_2[metadata_2$patiendId %in% l,]
    print(tmp)
    
    sink <- t(as.matrix(otus[patiendId == l,]))

    max_rel_abun <- c()
    for(j in 1:dim(sink)[2]){
      
      max_rel_abun[j] = max(sink[,j])/sum(sink[,j])
    }
    
    
    Dom_flag = 1
    D = min(which(max_rel_abun >= 0.3)) - 1
    
    if(D == Inf){
      
      Dom_flag = 0
      D = length(max_rel_abun)
      next
      
    }
    
    
    if(length(D) == 1 & D > 1){
      
      
      #Before domination
      sink <- t(as.matrix(otus[day == tmp$day[D] & patiendId == l,]))
      
      if(D > 2)
        sources = as.matrix(otus[day < tmp$day[D] & patiendId == l,])
      if(D == 2)
        sources = as.vector(otus[day < tmp$day[D] & patiendId == l,])
      
      
      if(D > 2){
        
        COVERAGE = 10000
        totalsource<-sources
        totalsource<-as.matrix(totalsource)
        totalsource<- rarefy(x = totalsource, maxdepth = COVERAGE)
        
        sources <- split(totalsource, seq(nrow(totalsource)))
        sources<-lapply(sources, as.matrix)
        dists<-lapply(sources, function(x) x/(sum(colSums(x))))
        totaldist<-t(Reduce("cbind", dists))
        
        N = 5
        k = 1
        totalnew_source = matrix(NA, ncol = dim(totalsource)[2], nrow = N*dim(totalsource)[1])
        
        for(i in 1:dim(totalsource)[1]){
          
          
          for(j in 1:N){
            print(k)
            totalnew_source[k,] = as.numeric(t(rmultinom(n = 1, size = sum(totalsource[i,]), prob=totaldist[i,] )))
            k = k+1
            
          }
          
        }
        
        sources = totalnew_source
        rownames(sources) = c(1:dim(sources)[1])
        
        sources<-rarefy(sources, COVERAGE)
        sink<-as.matrix(rarefy(sink, COVERAGE))
        
        em_res<-EM_results(source=sources, sinks = t(sink), env= rownames(sources), em_itr = EM_ITERATIONS, COVERAGE = COVERAGE)
        print(em_res$data_prop)
        
        em_results_BD[[l]] = em_res$data_prop
        Results_array_emnoise_BD[[l]] = em_res$data_prop[,1]
        Results_array_em_BD[[l]] = em_res$data_prop[,2]
        
      }
      
      if(D == 2){
        
        sources_tmp <- as.numeric(t(rmultinom(n = 1, size = sum(sources) , prob= sources/sum(sources)   )))
        
        sources_tmp_2 <- matrix(NA, ncol = length(sources), nrow = 2)
        sources_tmp_2[1,] <- as.numeric(sources)
        sources_tmp_2[2,] <- sources_tmp
        
        sources <- as.matrix(sources_tmp_2)
        
        COVERAGE = 10000
        totalsource<-sources
        totalsource<-as.matrix(totalsource)
        totalsource<- rarefy(x = totalsource, maxdepth = COVERAGE)
        
        sources <- split(totalsource, seq(nrow(totalsource)))
        sources<-lapply(sources, as.matrix)
        dists<-lapply(sources, function(x) x/(sum(colSums(x))))
        totaldist<-t(Reduce("cbind", dists))
        
        N = 5
        totalnew_source = matrix(NA, ncol = dim(totalsource)[2], nrow = N*dim(totalsource)[1])
        k = 1
        for(i in 1:dim(totalsource)[1]){
          
          
          for(j in 1:N){
            print(k)
            totalnew_source[k,] = as.numeric(t(rmultinom(n = 1, size = sum(totalsource[i,]), prob=totaldist[i,])))
            k = k+1
            
          }
          
        }
        
        sources = totalnew_source
        rownames(sources) = c(1:dim(sources)[1])
        
        # COVERAGE = min(rowSums(sources))
        sources<-rarefy(sources, COVERAGE)
        sink<-as.matrix(rarefy(sink, COVERAGE))
        
        em_res<-EM_results(source=sources, sinks = t(sink), env= rownames(sources), em_itr = EM_ITERATIONS, COVERAGE = COVERAGE)
        print(em_res$data_prop)
        
        em_results_BD[[l]] = em_res$data_prop
        Results_array_emnoise_BD[[l]] = em_res$data_prop[,1]
        Results_array_em_BD[[l]] = em_res$data_prop[,2]
        
      }
      
      
      
      if(Dom_flag == 1){
        #Domination event
        sink <- t(as.matrix(otus[day == tmp$day[D+1] & patiendId == l,]))
        sources = as.matrix(otus[day < tmp$day[D+1] & patiendId == l,])
        # str(sources)
        
        COVERAGE = 10000
        totalsource<-sources
        totalsource<-as.matrix(totalsource)
        totalsource<- rarefy(x = totalsource, maxdepth = COVERAGE)
        
        sources <- split(totalsource, seq(nrow(totalsource)))
        sources<-lapply(sources, as.matrix)
        dists<-lapply(sources, function(x) x/(sum(colSums(x))))
        totaldist<-t(Reduce("cbind", dists))
        
        N = 5
        totalnew_source = matrix(NA, ncol = dim(totalsource)[2], nrow = N*dim(totalsource)[1])
        k = 1
        for(i in 1:dim(totalsource)[1]){
          
          
          for(j in 1:N){
            print(k)
            totalnew_source[k,] = as.numeric(t(rmultinom(n = 1, size = sum(totalsource[i,]), prob=totaldist[i,])))
            k = k+1
            
          }
          
        }
        
        sources = totalnew_source
        rownames(sources) = c(1:dim(sources)[1])
        COVERAGE = min(rowSums(sources))
        sources<-rarefy(sources, COVERAGE)
        sink<-as.matrix(rarefy(sink, COVERAGE))
        
        em_res<-EM_results(source=sources, sinks = t(sink), env= rownames(sources), em_itr = EM_ITERATIONS, COVERAGE = COVERAGE)
        print(em_res$data_prop)
        
        em_results_AD[[l]] = em_res$data_prop
        Results_array_emnoise_AD[[l]] = em_res$data_prop[,1]
        Results_array_em_AD[[l]] = em_res$data_prop[,2]
      }

      
    }
    
    
  }
  
  
  unknown_emnoise_BD = c()
  unknown_emnoise_AD = c()
  
  for(i in 1:length(em_results_BD)){
    
    if(length(em_results_BD[[i]]) > 0){
      
      tmp = dim(em_results_BD[[i]])[1]
      unknown_emnoise_BD[i] = em_results_BD[[i]][tmp,1]
      
    }
    
    if(length(em_results_AD[[i]]) > 0){
      
      tmp_2 = dim(em_results_AD[[i]])[1]
      unknown_emnoise_AD[i] = em_results_AD[[i]][tmp_2,1]
      
    }
    
    
  }
  
  print(wilcox.test(na.omit(unknown_emnoise_BD), na.omit(unknown_emnoise_AD)))
  print(t.test(na.omit(unknown_emnoise_BD), na.omit(unknown_emnoise_AD)))
  
  
  unknown_dist_BD_AD = data.frame(c(rep(0, length(na.omit(unknown_emnoise_BD))), rep(1, length(na.omit(unknown_emnoise_BD)))),
                                  c(na.omit(unknown_emnoise_BD), na.omit(unknown_emnoise_AD)) )
  names(unknown_dist_BD_AD) = c("BD_AD", "unknown_source_proportion")

  
  toplot_Before_D <- unknown_dist_BD_AD[unknown_dist_BD_AD$BD_AD == 0,]
  # toplot_Delivery_1.melt <- melt(toplot_Delivery_1, variable.name= "Origin",
  #                                value.name = "Proportion")
  
  
  toplot_Before_plot <- ggplot(toplot_Before_D, aes(x= BD_AD, y = unknown_source_proportion)) + 
    geom_boxplot(show.legend=F) + theme_bw()+ ylim(0,1) + ggtitle("Unknown source distribution before Domination")
  
  
  toplot_After_D <- unknown_dist_BD_AD[unknown_dist_BD_AD$BD_AD == 1,]
  # toplot_Delivery_1.melt <- melt(toplot_Delivery_1, variable.name= "Origin",
  #                                value.name = "Proportion")
  
  toplot_After_plot <- ggplot(toplot_After_D, aes(x= BD_AD, y = unknown_source_proportion)) + 
    geom_boxplot(show.legend=F) + theme_bw()+ ylim(0,1) + ggtitle("Unknown source distribution after Domination")
  
  allo_HSCT_plot = plot_grid(toplot_Before_plot, toplot_After_plot,  nrow = 1, align = 'h')
  
  ggsave(filename="../results/allo-HSCT-bacteria_domination.png", plot = allo_HSCT_plot , dpi = 600, width = 8.75, height = 6.1, units = "in")

  