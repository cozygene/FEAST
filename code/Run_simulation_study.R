
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


#Simulation----
# rm(list = ls())
# gc()

source("SourceTracker.R")
source("src.R")
source_data<-readRDS("../data/60KFILT10KSOURCE.rds")


COVERAGE=10000
include_epsilon<-T
n_sources<-20
mixing_iterations<- 30
st_mixing_iterations<- 30
epmean<-0.01
set_ep_prop<- abs(rnorm(n = mixing_iterations, sd = 0.05, mean = epmean) + 0.01) ##proportion of the unknown source
sample_iterations<- 10
em_iterations<-100 ###stop if m changes by 10^-6 or 100 iterations
n_cores<-8
cl<- makeCluster(n_cores)
registerDoParallel(cl)
epstring="NOEPS"
if(include_epsilon==T){
  epstring<-substr(paste0(epmean),start=nchar(paste0(epmean))-1, stop=nchar(paste0(epmean)))
}
unk=T; dat="60KFILT10KSOURCE"; itr=sample_iterations; noise=T;plot = 0;
unknown=""
eps=include_epsilon; bs=F; clsinit=F; em_itr=em_iterations; num_sources=n_sources
noisestring<-"nonoise"
if(noise==T){
  noisestring<-"noise"
}
bootstring<-"nobs" ##"no bootstrapping"
if(bs==T){
  bootstring<-"bs"
}
clsstring<-"noclsinit"
if(clsinit==T){
  clsstring<-"clsinit"
}

epsilon<-0.0 #initiate to 0 in case eps is False
alpha1 <- alpha2 <- 0.001 #alphas for training
js_values = c() #Jensen-Shannon divergence values vector
unknown<-""; if(unk==T){  unknown<-"" }
st_R2_averages = c() ;st_runtimes = c(); st_ep_R2<-c();
cls_R2_averages = c(); cls_runtimes = c(); cls_ep_R2<-c();
em_R2_averages = c(); em_runtimes = c(); em_ep_R2<-c();
emnoise_R2_averages =c(); emnoise_runtimes = c(); emnoise_ep_R2=c();
m_matrices=list()
st_m_matrices<-list(); st_sources_R2s<-list()
cls_m_matrices<-list(); cls_sources_R2s<-list()
em_m_matrices<-list(); em_sources_R2s<-list()
emnoise_m_matrices=list(); emnoise_sources_R2s=list()
num_sources = n_sources



t = 1
dat = "M3"
source_data<-read_pseudo_data(dat)
itr = 100

for (it in 1:itr){
  
  
  print("Creating the known sources")
  envs_simulation <- 1:num_sources
  totalsource<-source_data[sample(c(1:5),num_sources,T),]
  totalsource<-as.matrix(totalsource)
  totalsource<- rarefy(x = totalsource, maxdepth = COVERAGE)
  
  sources <- split(totalsource, seq(nrow(totalsource)))
  sources<-lapply(sources, as.matrix)
  dists<-lapply(sources, function(x) x/(sum(colSums(x))))
  totaldist<-t(Reduce("cbind", dists))
  
  
  print("Adding the unknown source")
  epsource<-rpois(n = length(totalsource[1,]), lambda = 10) #use a different source for the unknown 
  epdist<-epsource/(sum(epsource))
  
  #Adding the unknown source

  if(include_epsilon == TRUE){
    
    totalsource_new = matrix(NA, nrow = (dim(totalsource)[1] + 1), ncol = (dim(totalsource)[2]))
    colnames(totalsource_new) = colnames(totalsource)
    rownames(totalsource_new) = seq(from = 1, to = (dim(totalsource)[1] + 1), by = 1)
    
    for(j in 1:(dim(totalsource)[1])){
      
      totalsource_new[j,] = totalsource[j,]
    }
    totalsource_new[(j+1),] = as.numeric(epsource)
  }
  
  totalsource = totalsource_new
  
  list_ms<-list()
  ep_values<-rep(0, mixing_iterations)
  
  
  ######change_C(X = totalnew_source[1,], newcov = COVERAGE) dont uncomment
  
  ##Known sources
  
  print("Known sources")
  new_source = list()
  totalnew_source = matrix(NA, ncol = dim(totalsource)[2], nrow = length(sources))
  for(j in 1:length(sources)){
    new_source[[j]] = t(rmultinom(n = 1, size = sum(totalsource[j,]), prob=totaldist[j,]))
    totalnew_source[j,] = as.numeric(t(rmultinom(n = 1, size = sum(totalsource[j,]), prob=totaldist[j,])))}
  new_source=lapply(new_source,t)
  # totalnew_source=t(apply(totalnew_source, 1, function(x) change_C(COVERAGE, x)))  ##COVERAGE CHANGE
  totalnew_source <- rarefy(x = totalnew_source, maxdepth = COVERAGE)
  new_source<- split(totalnew_source, seq(nrow(totalnew_source)))
  new_source<-lapply(new_source, as.matrix)
  
  
  print("Creating sinks - 1")
  
  newdists<-lapply(new_source, function(x) x/(sum(colSums(x))))
  totalnew_dist<-t(Reduce("cbind", newdists))
  list_ms<-list()
  sinks<-matrix(nrow = mixing_iterations, ncol = dim(totalsource)[2])
  
  print("Creating sinks - 2")
  for(k in 1:mixing_iterations){
    
    list_ms[[k]] = create_m(num_sources = num_sources, n = 3, EPSILON = max(0, set_ep_prop[k]))
    sinks[k,] = t(as.matrix(list_ms[[k]])) %*% totalsource
    sinks[k,] = as.numeric(t(rmultinom(n=1, size=COVERAGE, prob=sinks[k,]/sum(sinks[k,]))))
  }
  
  print("ms_mat")
  #The sink is created using the known and unknown sources (if include_epsilon == TRUE )
  sinks=round(sinks)
  m_matrices=append(m_matrices, list(mapply(rbind, list_ms)))
  ms_mat=c()
  if(include_epsilon==F){
    ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources, byrow = T)
  }else{
    ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
  }
  
  print("Jensen-Shannon divergence")
  
  # Jensen-Shannon Divergence between P and Q - new version
  JSDMatrix <- jsdmatrix(totalnew_source[c(1:num_sources),]) #instead of x
  JSDMatrix <- JSDMatrix/COVERAGE
  JS = mean(JSDMatrix[-which(JSDMatrix == 0)])
  print(JS)
  
  
  
  if(length(js_values) == 0) {
    
    if(JS < 0.1){
      
      
      
      t = t+1
      js_values = append(js_values, JS) 
      
      
      
      new_source_old = new_source
      totalnew_source_old = totalnew_source
      
      ms_mat<-c()
      if(include_epsilon==F){
        ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      
      
      
      #############START FEAST###########
      
      em_rts<-c(); emnoise_rts =c();
      predict_values_em<-c()
      predict_eps_em<-c()
      em_results<-c()
      emnoise_results<-c()
      predict_values_emnoise<-c()
      predict_eps_emnoise<-c()
      
      print("Start FEAST")
      
      for(i in 1:mixing_iterations){
        
        #Creating the unknown source per mixing iteration
        if(include_epsilon == TRUE){
          
          ##Adding the initial value of the unknown source for CLS and EM
          new_source_2 = list()
          totalnew_source_2 = matrix(NA, ncol = dim(totalnew_source_old)[2], nrow = ( dim(totalnew_source_old)[1] + 1))
          
          for(j in 1:n_sources){
            
            new_source_2[[j]] = new_source_old[[j]]
            totalnew_source_2[j,] = totalnew_source_old[j,]
          }
          
          #create unknown for each sink i
          unknown_source = unknown_initialize(sources = totalnew_source_old[c(1:n_sources),], sink = as.numeric(sinks[i,]), 
                                              n_sources = n_sources) 
          
          new_source_2[[j+1]] = round(unknown_source)
          totalnew_source_2[(j+1),] = round(unknown_source)
          
          
          totalnew_source = totalnew_source_2
          new_source=lapply(new_source_2,t)
          totalnew_source <- rarefy(x = totalnew_source, maxdepth = COVERAGE)
          new_source<- split(totalnew_source_2, seq(nrow(totalnew_source_2)))
          new_source<-lapply(new_source, as.matrix)
          envs_simulation <- c(1:(num_sources+1))
          
        }
        
        samps <- list()
        samps <- new_source
        samps<-lapply(samps, t)
        
        
        observed_samps <- samps
        observed_samps[[(num_sources + 1)]] = t(rep(0,  dim(samps[[1]])[2]))
        
        
        if(eps==T) {initalphs<-runif(num_sources+1, 0.0, 1.0)
        }else {initalphs<-runif(num_sources, 0.0, 1.0)}
        initalphs=initalphs/Reduce("+", initalphs)
        m_s<-unlist(list_ms[[i]])
        sink=t(as.matrix(sinks[i,]))
        em_start_time<-as.numeric(proc.time()[3])
        pred_em<-do_EM_basic(alphas=initalphs, sources=samps, sink=sink, iterations=em_iterations)
        em_results<-c(em_results, pred_em)
        em_rts=c(em_rts, as.numeric(proc.time()[3])-em_start_time)
        emnoise_start_time<-as.numeric(proc.time()[3])
        tmp<-do_EM(alphas=initalphs, sources=samps, sink=sink, iterations=em_iterations, observed=observed_samps)
        pred_emnoise = tmp$toret
        emnoise_rts=c(emnoise_rts, as.numeric(proc.time()[3])-emnoise_start_time)
        emnoise_results=c(emnoise_results, pred_emnoise)
      }
      
      
      em_predictions=c()
      if(include_epsilon==F){
        em_predictions <- matrix(em_results, nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        em_predictions <- matrix(em_results, nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      em_runtimes<-c(em_runtimes, mean(em_rts))
      emnoise_runtimes<-c(emnoise_runtimes, mean(emnoise_rts))
      
      em_R2s<-c()
      for(j in 1:dim(em_predictions)[2]){
        em_R2s[j] = cor(em_predictions[,j], ms_mat[,j])^2
      }
      em_R2_averages<-c(em_R2_averages, mean(na.omit(em_R2s)))
      em_m_matrices<-append(em_m_matrices, list(mapply(rbind, em_predictions)))
      em_sources_R2s<-append(em_sources_R2s, list(as.matrix(t(em_R2s))))
      emnoise_predictions=c()
      if(include_epsilon==F){
        emnoise_predictions <- matrix(emnoise_results, nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        emnoise_predictions <- matrix(emnoise_results, nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      emnoise_R2s<-c()
      for(j in 1:dim(em_predictions)[2]){
        emnoise_R2s[j] = cor(emnoise_predictions[,j], ms_mat[,j])^2
      }
      emnoise_R2_averages<-c(emnoise_R2_averages, mean(na.omit(emnoise_R2s)))
      emnoise_m_matrices<-append(emnoise_m_matrices, list(mapply(rbind, emnoise_predictions)))
      emnoise_sources_R2s<-append(emnoise_sources_R2s, list(as.matrix(t(emnoise_R2s))))
      
      
    #   print(paste("FEAST_R2_averages = ", round(emnoise_R2_averages, 4)))
      print("End FEAST")
      
      #############END FEAST#################
      
      
      
      if(ST_simulations_flag == 1){
        
        print("Start ST")
        
        #############ST START############
        envs_simulation_st = LETTERS[c(1:num_sources)]
        st <- sourcetracker(totalnew_source_old, envs_simulation_st, rarefaction_depth=COVERAGE)
        
        result<- foreach(i=1:st_mixing_iterations, .combine='c', .export=ls(.GlobalEnv)) %dopar% {
          sink=t(as.matrix(sinks[i,]))
          st_start_time<-as.numeric(proc.time()[3])
          results <- predict(st,sink, alpha1=alpha1, alpha2=alpha2)
          predict <- c(results$proportions) #predicted proportions in sample
          c(list(predict),list(as.numeric(proc.time()[3])-st_start_time))
        }
        
        
        st_rts<-split(result, 1:2)[[2]]
        st_rts<-unlist(st_rts)
        st_runtimes<-c(st_runtimes, mean(st_rts))
        st_hold=split(result, 1:2)[[1]]
        st_hold=unlist(st_hold)
        #####st_predictions=split(result, 1:(num_sources+1))dont uncomment
        st_predictions <- matrix(st_hold, nrow = st_mixing_iterations, ncol = num_sources+1, byrow = T)
        st_R2s<-c()
        if(include_epsilon==F){
          st_predictions=st_predictions[,1:num_sources]
        }
        for(j in 1:dim(st_predictions)[2]){
          st_R2s[j] = cor(st_predictions[c(1:st_mixing_iterations),j], ms_mat[c(1:st_mixing_iterations),j])^2
        }
        st_R2_averages = append(st_R2_averages,mean(na.omit(st_R2s)))
        newstmat<-as.matrix(mapply(rbind, st_predictions))
        st_m_matrices<-append(st_m_matrices, list(newstmat))
        st_sources_R2s<-append(st_sources_R2s, list(as.matrix(t(st_R2s))))
        
        ############ST###########
        print("End ST")
        
      }
      
      
      break;
    }
    
    
  }
  
  
}

dat = "F4"
source_data<-read_pseudo_data(dat)
itr = 1000


for (it in 1:itr){
  
  
  print("Creating the known sources")
  envs_simulation <- 1:num_sources
  totalsource<-source_data[sample(1:dim(source_data)[1],num_sources,F),]
  totalsource<-as.matrix(totalsource)
  totalsource<- rarefy(x = totalsource, maxdepth = COVERAGE)
  
  sources <- split(totalsource, seq(nrow(totalsource)))
  sources<-lapply(sources, as.matrix)
  dists<-lapply(sources, function(x) x/(sum(colSums(x))))
  totaldist<-t(Reduce("cbind", dists))
  
  
  print("Adding the unknown source")
  epsource<-rpois(n = length(totalsource[1,]), lambda = 10) #use a different source for the unknown 
  epdist<-epsource/(sum(epsource))
  
  # epsource<-source[sample(1:dim(source)[1],1,F),]
  # epdist<-epsource/(sum(colSums(epsource)))
  
  
  #Adding the unknown source

  if(include_epsilon == TRUE){
    
    totalsource_new = matrix(NA, nrow = (dim(totalsource)[1] + 1), ncol = (dim(totalsource)[2]))
    colnames(totalsource_new) = colnames(totalsource)
    rownames(totalsource_new) = seq(from = 1, to = (dim(totalsource)[1] + 1), by = 1)
    
    for(j in 1:(dim(totalsource)[1])){
      
      totalsource_new[j,] = totalsource[j,]
    }
    totalsource_new[(j+1),] = as.numeric(epsource)
  }
  
  totalsource = totalsource_new
  
  list_ms<-list()
  ep_values<-rep(0, mixing_iterations)
  
  
  ######change_C(X = totalnew_source[1,], newcov = COVERAGE) dont uncomment
  
  ##Known sources
  
  print("Known sources")
  new_source = list()
  totalnew_source = matrix(NA, ncol = dim(totalsource)[2], nrow = length(sources))
  for(j in 1:length(sources)){
    new_source[[j]] = t(rmultinom(n = 1, size = sum(totalsource[j,]), prob=totaldist[j,]))
    totalnew_source[j,] = as.numeric(t(rmultinom(n = 1, size = sum(totalsource[j,]), prob=totaldist[j,])))}
  new_source=lapply(new_source,t)
  # totalnew_source=t(apply(totalnew_source, 1, function(x) change_C(COVERAGE, x)))  ##COVERAGE CHANGE
  totalnew_source <- rarefy(x = totalnew_source, maxdepth = COVERAGE)
  new_source<- split(totalnew_source, seq(nrow(totalnew_source)))
  new_source<-lapply(new_source, as.matrix)
  
  
  print("Creating sinks - 1")
  
  newdists<-lapply(new_source, function(x) x/(sum(colSums(x))))
  totalnew_dist<-t(Reduce("cbind", newdists))
  list_ms<-list()
  sinks<-matrix(nrow = mixing_iterations, ncol = dim(totalsource)[2])
  
  print("Creating sinks - 2")
  for(k in 1:mixing_iterations){
    
    list_ms[[k]] = create_m(num_sources = num_sources, n = 3, EPSILON = max(0, set_ep_prop[k]))
    sinks[k,] = t(as.matrix(list_ms[[k]])) %*% totalsource
    sinks[k,] = as.numeric(t(rmultinom(n=1, size=COVERAGE, prob=sinks[k,]/sum(sinks[k,]))))
  }
  
  print("ms_mat")
  #The sink is created using the known and unknown sources (if include_epsilon == TRUE )
  sinks=round(sinks)
  # sinks=t(apply(sinks, 1, function(x) change_C(COVERAGE, x)))
  # sinks = rarefy(x = sinks, maxdepth = COVERAGE)
  m_matrices=append(m_matrices, list(mapply(rbind, list_ms)))
  ms_mat=c()
  if(include_epsilon==F){
    ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources, byrow = T)
  }else{
    ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
  }
  
  print("Jensen-Shannon divergence")
  #Jensen-Shannon divergence calculation - OLD VERSION
  # weights<-rep(1/num_sources, num_sources)
  # JS = mult_JSD(as.matrix(totaldist), q=weights); print(cat('jsd: ', JS))
  # js_values = append(js_values, JS)
  
  # Jensen-Shannon Divergence between P and Q - new version
  #x <- totalsource[c(1:num_sources),]
  JSDMatrix <- jsdmatrix(totalnew_source[c(1:num_sources),]) #instead of x
  JSDMatrix <- JSDMatrix/COVERAGE
  ##not good because the coverage is a huge range...
  #JSDMatrix<-sapply(1:num_sources, function(z) JSDMatrix[z,]/sum(totalsource[z,]))
  JS = mean(JSDMatrix[-which(JSDMatrix == 0)])
  print(JS)
  ###save the JSD values
  # plot a heatmap of the corresponding JSD matrix
  # heatmap(JSDMatrix)
  
  
  
  if(length(js_values) == 0) {
    
    if(JS < 0.19){
      
      
      
      t = t+1
      js_values = append(js_values, JS) 
      
      
      #############CLS###########
      
      
      #predict_values_cls <- c(); predict_eps_cls<-c(); cls_runtime<-c()
      #cls_start_time<-as.numeric(proc.time()[3])
      #bootsource<-matrix(NA,ncol = dim(sources[[1]])[1])
      #if(bs==T){ ###if we want to use bootstrapping 
      #  envs_simulation<-c()
      #  for(i in 1:length(sources)){
      #    bootsource=rbind(bootsource,t(as.matrix(sources[[i]])))
      #    bootsource=rbind(bootsource, t(
      #      rmultinom(n = 99, size = sum(totalnew_source[i,]), prob=totalnew_dist[i,])))
      #    envs_simulation=c(envs_simulation, rep(i, 100))
      #  }
      #  totalnew_source<-bootsource[complete.cases(bootsource),]
      #  rownames(totalnew_source)<-make.names(rownames(totalnew_source), unique = T)
      #  envs_simulation=factor(envs_simulation)
      #}
      
      new_source_old = new_source
      totalnew_source_old = totalnew_source
      
      
      
      
      ms_mat<-c()
      if(include_epsilon==F){
        ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      
      #############CLS END#############
      
      
      #############EM###########
      
      em_rts<-c(); emnoise_rts =c();
      predict_values_em<-c()
      predict_eps_em<-c()
      em_results<-c()
      emnoise_results<-c()
      predict_values_emnoise<-c()
      predict_eps_emnoise<-c()
      
      print("Start FEAST")
      
      for(i in 1:mixing_iterations){
        
        #Creating the unknown source per mixing iteration
        if(include_epsilon == TRUE){
          
          ##Adding the initial value of the unknown source for CLS and EM
          new_source_2 = list()
          totalnew_source_2 = matrix(NA, ncol = dim(totalnew_source_old)[2], nrow = ( dim(totalnew_source_old)[1] + 1))
          
          for(j in 1:n_sources){
            
            new_source_2[[j]] = new_source_old[[j]]
            totalnew_source_2[j,] = totalnew_source_old[j,]
          }
          
          #create unknown for each sink i
          unknown_source = unknown_initialize(sources = totalnew_source_old[c(1:n_sources),], sink = as.numeric(sinks[i,]), 
                                              n_sources = n_sources) 
          # new_source_2[[j+1]] = unknown_source
          # totalnew_source_2[(j+1),] = unknown_source
          # plot(as.numeric(unknown_source))
          
          # tmp_unk = change_C(newcov = COVERAGE, X = unknown_source)
          # tmp_unk =  rarefy(x = unknown_source, maxdepth = COVERAGE)
          # unknown_source = tmp_unk
          # new_source_2[[j+1]] = t(round(unknown_source))
          # source_2[[j+1]] = round(unknown_source)
          new_source_2[[j+1]] = round(unknown_source)
          totalnew_source_2[(j+1),] = round(unknown_source)
          
          
          totalnew_source = totalnew_source_2
          
          new_source=lapply(new_source_2,t)
          # totalnew_source <- t(apply(totalnew_source, 1, function(x) change_C(COVERAGE, x)))  ##COVERAGE CHANGE
          totalnew_source <- rarefy(x = totalnew_source, maxdepth = COVERAGE)
          new_source<- split(totalnew_source_2, seq(nrow(totalnew_source_2)))
          new_source<-lapply(new_source, as.matrix)
          
          # View(data.frame(t(totalnew_source), sinks[i,]))
          
          envs_simulation <- c(1:(num_sources+1))
          
          # View(data.frame(t(totalnew_source), sinks[i,]) )
        }
        
        samps <- list()
        samps <- new_source
        samps<-lapply(samps, t)
        
        # num_sources = 5
        observed_samps <- samps
        observed_samps[[(num_sources + 1)]] = t(rep(0,  dim(samps[[1]])[2]))
        
        # observed_samps[[(num_sources + 1)]] = t(rpois(lambda = 2,n = dim(samps[[1]])[2]))
        
        # observed_samps<-lapply(observed_samps, t)
        # str(observed_samps)
        
        # View(data.frame(t(totalsource),  sinks[1,]))
        
        
        # num_sources = 5
        if(eps==T) {initalphs<-runif(num_sources+1, 0.0, 1.0)
        }else {initalphs<-runif(num_sources, 0.0, 1.0)}
        initalphs=initalphs/Reduce("+", initalphs)
        m_s<-unlist(list_ms[[i]])
        # if(clsinit==T){initalphs<-c(predict_values_cls[i], 1-predict_values_cls[i])}
        sink=t(as.matrix(sinks[i,]))
        em_start_time<-as.numeric(proc.time()[3])
        pred_em<-do_EM_basic(alphas=initalphs, sources=samps, sink=sink, iterations=em_iterations)
        em_results<-c(em_results, pred_em)
        em_rts=c(em_rts, as.numeric(proc.time()[3])-em_start_time)
        emnoise_start_time<-as.numeric(proc.time()[3])
        tmp<-do_EM(alphas=initalphs, sources=samps, sink=sink, iterations=em_iterations, observed=observed_samps)
        pred_emnoise = tmp$toret
        emnoise_rts=c(emnoise_rts, as.numeric(proc.time()[3])-emnoise_start_time)
        emnoise_results=c(emnoise_results, pred_emnoise)
      }
      
      
      em_predictions=c()
      if(include_epsilon==F){
        em_predictions <- matrix(em_results, nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        em_predictions <- matrix(em_results, nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      em_runtimes<-c(em_runtimes, mean(em_rts))
      emnoise_runtimes<-c(emnoise_runtimes, mean(emnoise_rts))
      
      em_R2s<-c()
      for(j in 1:dim(em_predictions)[2]){
        em_R2s[j] = cor(em_predictions[,j], ms_mat[,j])^2
      }
      em_R2_averages<-c(em_R2_averages, mean(na.omit(em_R2s)))
      em_m_matrices<-append(em_m_matrices, list(mapply(rbind, em_predictions)))
      em_sources_R2s<-append(em_sources_R2s, list(as.matrix(t(em_R2s))))
      emnoise_predictions=c()
      if(include_epsilon==F){
        emnoise_predictions <- matrix(emnoise_results, nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        emnoise_predictions <- matrix(emnoise_results, nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      emnoise_R2s<-c()
      for(j in 1:dim(em_predictions)[2]){
        emnoise_R2s[j] = cor(emnoise_predictions[,j], ms_mat[,j])^2
      }
      emnoise_R2_averages<-c(emnoise_R2_averages, mean(na.omit(emnoise_R2s)))
      emnoise_m_matrices<-append(emnoise_m_matrices, list(mapply(rbind, emnoise_predictions)))
      emnoise_sources_R2s<-append(emnoise_sources_R2s, list(as.matrix(t(emnoise_R2s))))
      
      
    #   print(paste("FEAST_R2_averages = ", round(emnoise_R2_averages, 4)))
      print("End FEAST")
      
      #############EM#################
      
      if(ST_simulations_flag == 1){
        
        print("Start ST")
        
        #############ST START############
        envs_simulation_st = LETTERS[c(1:num_sources)]
        st <- sourcetracker(totalnew_source_old, envs_simulation_st, rarefaction_depth=COVERAGE)
        
        result<- foreach(i=1:st_mixing_iterations, .combine='c', .export=ls(.GlobalEnv)) %dopar% {
          sink=t(as.matrix(sinks[i,]))
          st_start_time<-as.numeric(proc.time()[3])
          results <- predict(st,sink, alpha1=alpha1, alpha2=alpha2)
          predict <- c(results$proportions) #predicted proportions in sample
          c(list(predict),list(as.numeric(proc.time()[3])-st_start_time))
        }
        
        
        st_rts<-split(result, 1:2)[[2]]
        st_rts<-unlist(st_rts)
        st_runtimes<-c(st_runtimes, mean(st_rts))
        st_hold=split(result, 1:2)[[1]]
        st_hold=unlist(st_hold)
        #####st_predictions=split(result, 1:(num_sources+1))dont uncomment
        st_predictions <- matrix(st_hold, nrow = st_mixing_iterations, ncol = num_sources+1, byrow = T)
        st_R2s<-c()
        if(include_epsilon==F){
          st_predictions=st_predictions[,1:num_sources]
        }
        for(j in 1:dim(st_predictions)[2]){
          st_R2s[j] = cor(st_predictions[c(1:st_mixing_iterations),j], ms_mat[c(1:st_mixing_iterations),j])^2
        }
        st_R2_averages = append(st_R2_averages,mean(na.omit(st_R2s)))
        newstmat<-as.matrix(mapply(rbind, st_predictions))
        st_m_matrices<-append(st_m_matrices, list(newstmat))
        st_sources_R2s<-append(st_sources_R2s, list(as.matrix(t(st_R2s))))
        
        ############ST###########
        print("End ST")
        
      }
      
      
    }
    
  }
  
  if(length(js_values) > 0) {
    
    if(JS > (js_values[t - 1] + 0.07) & JS < (js_values[t - 1] + 0.1)){
      
      t = t+1
      js_values = append(js_values, JS) 
      
      
      new_source_old = new_source
      totalnew_source_old = totalnew_source
      
      ms_mat<-c()
      if(include_epsilon==F){
        ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      
      
      
      #############EM###########
      
      em_rts<-c(); emnoise_rts =c();
      predict_values_em<-c()
      predict_eps_em<-c()
      em_results<-c()
      emnoise_results<-c()
      predict_values_emnoise<-c()
      predict_eps_emnoise<-c()
      
      print("Start FEAST")
      
      for(i in 1:mixing_iterations){
        
        #Creating the unknown source per mixing iteration
        if(include_epsilon == TRUE){
          
          ##Adding the initial value of the unknown source for CLS and EM
          new_source_2 = list()
          totalnew_source_2 = matrix(NA, ncol = dim(totalnew_source_old)[2], nrow = ( dim(totalnew_source_old)[1] + 1))
          
          for(j in 1:n_sources){
            
            new_source_2[[j]] = new_source_old[[j]]
            totalnew_source_2[j,] = totalnew_source_old[j,]
          }
          
          #create unknown for each sink i
          unknown_source = unknown_initialize(sources = totalnew_source_old[c(1:n_sources),], sink = as.numeric(sinks[i,]), 
                                              n_sources = n_sources) 
          # new_source_2[[j+1]] = unknown_source
          # totalnew_source_2[(j+1),] = unknown_source
          # plot(as.numeric(unknown_source))
          
          # tmp_unk = change_C(newcov = COVERAGE, X = unknown_source)
          # tmp_unk =  rarefy(x = unknown_source, maxdepth = COVERAGE)
          # unknown_source = tmp_unk
          # new_source_2[[j+1]] = t(round(unknown_source))
          # source_2[[j+1]] = round(unknown_source)
          new_source_2[[j+1]] = round(unknown_source)
          totalnew_source_2[(j+1),] = round(unknown_source)
          
          
          totalnew_source = totalnew_source_2
          
          new_source=lapply(new_source_2,t)
          # totalnew_source <- t(apply(totalnew_source, 1, function(x) change_C(COVERAGE, x)))  ##COVERAGE CHANGE
          totalnew_source <- rarefy(x = totalnew_source, maxdepth = COVERAGE)
          new_source<- split(totalnew_source_2, seq(nrow(totalnew_source_2)))
          new_source<-lapply(new_source, as.matrix)
          
          # View(data.frame(t(totalnew_source), sinks[i,]))
          
          envs_simulation <- c(1:(num_sources+1))
          
          # View(data.frame(t(totalnew_source), sinks[i,]) )
        }
        
        samps <- list()
        samps <- new_source
        samps<-lapply(samps, t)
        
        # num_sources = 5
        observed_samps <- samps
        observed_samps[[(num_sources + 1)]] = t(rep(0,  dim(samps[[1]])[2]))
        
        # observed_samps[[(num_sources + 1)]] = t(rpois(lambda = 2,n = dim(samps[[1]])[2]))
        
        # observed_samps<-lapply(observed_samps, t)
        # str(observed_samps)
        
        # View(data.frame(t(totalsource),  sinks[1,]))
        
        
        # num_sources = 5
        if(eps==T) {initalphs<-runif(num_sources+1, 0.0, 1.0)
        }else {initalphs<-runif(num_sources, 0.0, 1.0)}
        initalphs=initalphs/Reduce("+", initalphs)
        m_s<-unlist(list_ms[[i]])
        # if(clsinit==T){initalphs<-c(predict_values_cls[i], 1-predict_values_cls[i])}
        sink=t(as.matrix(sinks[i,]))
        em_start_time<-as.numeric(proc.time()[3])
        pred_em<-do_EM_basic(alphas=initalphs, sources=samps, sink=sink, iterations=em_iterations)
        em_results<-c(em_results, pred_em)
        em_rts=c(em_rts, as.numeric(proc.time()[3])-em_start_time)
        emnoise_start_time<-as.numeric(proc.time()[3])
        tmp<-do_EM(alphas=initalphs, sources=samps, sink=sink, iterations=em_iterations, observed=observed_samps)
        pred_emnoise = tmp$toret
        emnoise_rts=c(emnoise_rts, as.numeric(proc.time()[3])-emnoise_start_time)
        emnoise_results=c(emnoise_results, pred_emnoise)
      }
      
      
      em_predictions=c()
      if(include_epsilon==F){
        em_predictions <- matrix(em_results, nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        em_predictions <- matrix(em_results, nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      em_runtimes<-c(em_runtimes, mean(em_rts))
      emnoise_runtimes<-c(emnoise_runtimes, mean(emnoise_rts))
      
      em_R2s<-c()
      for(j in 1:dim(em_predictions)[2]){
        em_R2s[j] = cor(em_predictions[,j], ms_mat[,j])^2
      }
      em_R2_averages<-c(em_R2_averages, mean(na.omit(em_R2s)))
      em_m_matrices<-append(em_m_matrices, list(mapply(rbind, em_predictions)))
      em_sources_R2s<-append(em_sources_R2s, list(as.matrix(t(em_R2s))))
      emnoise_predictions=c()
      if(include_epsilon==F){
        emnoise_predictions <- matrix(emnoise_results, nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        emnoise_predictions <- matrix(emnoise_results, nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      emnoise_R2s<-c()
      for(j in 1:dim(em_predictions)[2]){
        emnoise_R2s[j] = cor(emnoise_predictions[,j], ms_mat[,j])^2
      }
      emnoise_R2_averages<-c(emnoise_R2_averages, mean(na.omit(emnoise_R2s)))
      emnoise_m_matrices<-append(emnoise_m_matrices, list(mapply(rbind, emnoise_predictions)))
      emnoise_sources_R2s<-append(emnoise_sources_R2s, list(as.matrix(t(emnoise_R2s))))
      
      
    #   print(paste("FEAST_R2_averages = ", round(emnoise_R2_averages, 4)))
      print("End FEAST")
      
      #############EM#################
      
      if(ST_simulations_flag == 1){
        
        print("Start ST")
        
        #############ST START############
        envs_simulation_st = LETTERS[c(1:num_sources)]
        st <- sourcetracker(totalnew_source_old, envs_simulation_st, rarefaction_depth=COVERAGE)
        
        result<- foreach(i=1:st_mixing_iterations, .combine='c', .export=ls(.GlobalEnv)) %dopar% {
          sink=t(as.matrix(sinks[i,]))
          st_start_time<-as.numeric(proc.time()[3])
          results <- predict(st,sink, alpha1=alpha1, alpha2=alpha2)
          predict <- c(results$proportions) #predicted proportions in sample
          c(list(predict),list(as.numeric(proc.time()[3])-st_start_time))
        }
        
        
        st_rts<-split(result, 1:2)[[2]]
        st_rts<-unlist(st_rts)
        st_runtimes<-c(st_runtimes, mean(st_rts))
        st_hold=split(result, 1:2)[[1]]
        st_hold=unlist(st_hold)
        #####st_predictions=split(result, 1:(num_sources+1))dont uncomment
        st_predictions <- matrix(st_hold, nrow = st_mixing_iterations, ncol = num_sources+1, byrow = T)
        st_R2s<-c()
        if(include_epsilon==F){
          st_predictions=st_predictions[,1:num_sources]
        }
        for(j in 1:dim(st_predictions)[2]){
          st_R2s[j] = cor(st_predictions[c(1:st_mixing_iterations),j], ms_mat[c(1:st_mixing_iterations),j])^2
        }
        st_R2_averages = append(st_R2_averages,mean(na.omit(st_R2s)))
        newstmat<-as.matrix(mapply(rbind, st_predictions))
        st_m_matrices<-append(st_m_matrices, list(newstmat))
        st_sources_R2s<-append(st_sources_R2s, list(as.matrix(t(st_R2s))))
        
        ############ST###########
        
        print("End ST")
        
      }
    }
    
  }  
  
  
}

dat = "DB"
source_data<-read_pseudo_data(dat)
itr = 1000

for (it in 1:itr){
  
  
  print("Creating the known sources")
  envs_simulation <- 1:num_sources
  totalsource<-source_data[sample(1:dim(source_data)[1],num_sources,F),]
  totalsource<-as.matrix(totalsource)
  totalsource<- rarefy(x = totalsource, maxdepth = COVERAGE)
  
  sources <- split(totalsource, seq(nrow(totalsource)))
  sources<-lapply(sources, as.matrix)
  dists<-lapply(sources, function(x) x/(sum(colSums(x))))
  totaldist<-t(Reduce("cbind", dists))
  
  
  # print(str(totalsource))
  
  print("Adding the unknown source")
  epsource<-rpois(n = length(totalsource[1,]), lambda = 10) #use a different source for the unknown 
  epdist<-epsource/(sum(epsource))
  
  # epsource<-source[sample(1:dim(source)[1],1,F),]
  # epdist<-epsource/(sum(colSums(epsource)))
  
  
  #Adding the unknown source

  if(include_epsilon == TRUE){
    
    totalsource_new = matrix(NA, nrow = (dim(totalsource)[1] + 1), ncol = (dim(totalsource)[2]))
    colnames(totalsource_new) = colnames(totalsource)
    rownames(totalsource_new) = seq(from = 1, to = (dim(totalsource)[1] + 1), by = 1)
    
    for(j in 1:(dim(totalsource)[1])){
      
      totalsource_new[j,] = totalsource[j,]
    }
    totalsource_new[(j+1),] = as.numeric(epsource)
  }
  
  totalsource = totalsource_new
  
  list_ms<-list()
  ep_values<-rep(0, mixing_iterations)
  
  
  ######change_C(X = totalnew_source[1,], newcov = COVERAGE) dont uncomment
  
  ##Known sources
  
  print("Known sources")
  new_source = list()
  totalnew_source = matrix(NA, ncol = dim(totalsource)[2], nrow = length(sources))
  for(j in 1:length(sources)){
    new_source[[j]] = t(rmultinom(n = 1, size = sum(totalsource[j,]), prob=totaldist[j,]))
    totalnew_source[j,] = as.numeric(t(rmultinom(n = 1, size = sum(totalsource[j,]), prob=totaldist[j,])))}
  new_source=lapply(new_source,t)
  # totalnew_source=t(apply(totalnew_source, 1, function(x) change_C(COVERAGE, x)))  ##COVERAGE CHANGE
  totalnew_source <- rarefy(x = totalnew_source, maxdepth = COVERAGE)
  new_source<- split(totalnew_source, seq(nrow(totalnew_source)))
  new_source<-lapply(new_source, as.matrix)
  
  
  print("Creating sinks - 1")
  
  newdists<-lapply(new_source, function(x) x/(sum(colSums(x))))
  totalnew_dist<-t(Reduce("cbind", newdists))
  list_ms<-list()
  sinks<-matrix(nrow = mixing_iterations, ncol = dim(totalsource)[2])
  
  print("Creating sinks - 2")
  for(k in 1:mixing_iterations){
    
    list_ms[[k]] = create_m(num_sources = num_sources, n = 3, EPSILON = max(0, set_ep_prop[k]))
    sinks[k,] = t(as.matrix(list_ms[[k]])) %*% totalsource
    sinks[k,] = as.numeric(t(rmultinom(n=1, size=COVERAGE, prob=sinks[k,]/sum(sinks[k,]))))
  }
  
  print("ms_mat")
  #The sink is created using the known and unknown sources (if include_epsilon == TRUE )
  sinks=round(sinks)
  # sinks=t(apply(sinks, 1, function(x) change_C(COVERAGE, x)))
  # sinks = rarefy(x = sinks, maxdepth = COVERAGE)
  m_matrices=append(m_matrices, list(mapply(rbind, list_ms)))
  ms_mat=c()
  if(include_epsilon==F){
    ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources, byrow = T)
  }else{
    ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
  }
  
  print("Jensen-Shannon divergence")
  #Jensen-Shannon divergence calculation - OLD VERSION
  # weights<-rep(1/num_sources, num_sources)
  # JS = mult_JSD(as.matrix(totaldist), q=weights); print(cat('jsd: ', JS))
  # js_values = append(js_values, JS)
  
  # Jensen-Shannon Divergence between P and Q - new version
  #x <- totalsource[c(1:num_sources),]
  JSDMatrix <- jsdmatrix(totalnew_source[c(1:num_sources),]) #instead of x
  JSDMatrix <- JSDMatrix/COVERAGE
  ##not good because the coverage is a huge range...
  #JSDMatrix<-sapply(1:num_sources, function(z) JSDMatrix[z,]/sum(totalsource[z,]))
  JS = mean(JSDMatrix[-which(JSDMatrix == 0)])
  print(JS)
  ###save the JSD values
  # plot a heatmap of the corresponding JSD matrix
  # heatmap(JSDMatrix)
  
  
  
  if(length(js_values) == 0) {
    
    if(JS < 0.19){
      
      
      
      t = t+1
      js_values = append(js_values, JS) 
      
      
      #############CLS###########
      
      
      #predict_values_cls <- c(); predict_eps_cls<-c(); cls_runtime<-c()
      #cls_start_time<-as.numeric(proc.time()[3])
      #bootsource<-matrix(NA,ncol = dim(sources[[1]])[1])
      #if(bs==T){ ###if we want to use bootstrapping 
      #  envs_simulation<-c()
      #  for(i in 1:length(sources)){
      #    bootsource=rbind(bootsource,t(as.matrix(sources[[i]])))
      #    bootsource=rbind(bootsource, t(
      #      rmultinom(n = 99, size = sum(totalnew_source[i,]), prob=totalnew_dist[i,])))
      #    envs_simulation=c(envs_simulation, rep(i, 100))
      #  }
      #  totalnew_source<-bootsource[complete.cases(bootsource),]
      #  rownames(totalnew_source)<-make.names(rownames(totalnew_source), unique = T)
      #  envs_simulation=factor(envs_simulation)
      #}
      
      new_source_old = new_source
      totalnew_source_old = totalnew_source
      
      
      
      
      ms_mat<-c()
      if(include_epsilon==F){
        ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      
      #############CLS END#############
      
      
      #############EM###########
      
      em_rts<-c(); emnoise_rts =c();
      predict_values_em<-c()
      predict_eps_em<-c()
      em_results<-c()
      emnoise_results<-c()
      predict_values_emnoise<-c()
      predict_eps_emnoise<-c()
      
      print("Start FEAST")
      
      for(i in 1:mixing_iterations){
        
        #Creating the unknown source per mixing iteration
        if(include_epsilon == TRUE){
          
          ##Adding the initial value of the unknown source for CLS and EM
          new_source_2 = list()
          totalnew_source_2 = matrix(NA, ncol = dim(totalnew_source_old)[2], nrow = ( dim(totalnew_source_old)[1] + 1))
          
          for(j in 1:n_sources){
            
            new_source_2[[j]] = new_source_old[[j]]
            totalnew_source_2[j,] = totalnew_source_old[j,]
          }
          
          #create unknown for each sink i
          unknown_source = unknown_initialize(sources = totalnew_source_old[c(1:n_sources),], sink = as.numeric(sinks[i,]), 
                                              n_sources = n_sources) 
          # new_source_2[[j+1]] = unknown_source
          # totalnew_source_2[(j+1),] = unknown_source
          # plot(as.numeric(unknown_source))
          
          # tmp_unk = change_C(newcov = COVERAGE, X = unknown_source)
          # tmp_unk =  rarefy(x = unknown_source, maxdepth = COVERAGE)
          # unknown_source = tmp_unk
          # new_source_2[[j+1]] = t(round(unknown_source))
          # source_2[[j+1]] = round(unknown_source)
          new_source_2[[j+1]] = round(unknown_source)
          totalnew_source_2[(j+1),] = round(unknown_source)
          
          
          totalnew_source = totalnew_source_2
          
          new_source=lapply(new_source_2,t)
          # totalnew_source <- t(apply(totalnew_source, 1, function(x) change_C(COVERAGE, x)))  ##COVERAGE CHANGE
          totalnew_source <- rarefy(x = totalnew_source, maxdepth = COVERAGE)
          new_source<- split(totalnew_source_2, seq(nrow(totalnew_source_2)))
          new_source<-lapply(new_source, as.matrix)
          
          # View(data.frame(t(totalnew_source), sinks[i,]))
          
          envs_simulation <- c(1:(num_sources+1))
          
          # View(data.frame(t(totalnew_source), sinks[i,]) )
        }
        
        samps <- list()
        samps <- new_source
        samps<-lapply(samps, t)
        
        # num_sources = 5
        observed_samps <- samps
        observed_samps[[(num_sources + 1)]] = t(rep(0,  dim(samps[[1]])[2]))
        
        # observed_samps[[(num_sources + 1)]] = t(rpois(lambda = 2,n = dim(samps[[1]])[2]))
        
        # observed_samps<-lapply(observed_samps, t)
        # str(observed_samps)
        
        # View(data.frame(t(totalsource),  sinks[1,]))
        
        
        # num_sources = 5
        if(eps==T) {initalphs<-runif(num_sources+1, 0.0, 1.0)
        }else {initalphs<-runif(num_sources, 0.0, 1.0)}
        initalphs=initalphs/Reduce("+", initalphs)
        m_s<-unlist(list_ms[[i]])
        # if(clsinit==T){initalphs<-c(predict_values_cls[i], 1-predict_values_cls[i])}
        sink=t(as.matrix(sinks[i,]))
        em_start_time<-as.numeric(proc.time()[3])
        pred_em<-do_EM_basic(alphas=initalphs, sources=samps, sink=sink, iterations=em_iterations)
        em_results<-c(em_results, pred_em)
        em_rts=c(em_rts, as.numeric(proc.time()[3])-em_start_time)
        emnoise_start_time<-as.numeric(proc.time()[3])
        tmp<-do_EM(alphas=initalphs, sources=samps, sink=sink, iterations=em_iterations, observed=observed_samps)
        pred_emnoise = tmp$toret
        emnoise_rts=c(emnoise_rts, as.numeric(proc.time()[3])-emnoise_start_time)
        emnoise_results=c(emnoise_results, pred_emnoise)
      }
      
      
      em_predictions=c()
      if(include_epsilon==F){
        em_predictions <- matrix(em_results, nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        em_predictions <- matrix(em_results, nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      em_runtimes<-c(em_runtimes, mean(em_rts))
      emnoise_runtimes<-c(emnoise_runtimes, mean(emnoise_rts))
      
      em_R2s<-c()
      for(j in 1:dim(em_predictions)[2]){
        em_R2s[j] = cor(em_predictions[,j], ms_mat[,j])^2
      }
      em_R2_averages<-c(em_R2_averages, mean(na.omit(em_R2s)))
      em_m_matrices<-append(em_m_matrices, list(mapply(rbind, em_predictions)))
      em_sources_R2s<-append(em_sources_R2s, list(as.matrix(t(em_R2s))))
      emnoise_predictions=c()
      if(include_epsilon==F){
        emnoise_predictions <- matrix(emnoise_results, nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        emnoise_predictions <- matrix(emnoise_results, nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      emnoise_R2s<-c()
      for(j in 1:dim(em_predictions)[2]){
        emnoise_R2s[j] = cor(emnoise_predictions[,j], ms_mat[,j])^2
      }
      emnoise_R2_averages<-c(emnoise_R2_averages, mean(na.omit(emnoise_R2s)))
      emnoise_m_matrices<-append(emnoise_m_matrices, list(mapply(rbind, emnoise_predictions)))
      emnoise_sources_R2s<-append(emnoise_sources_R2s, list(as.matrix(t(emnoise_R2s))))
      
      
    #   print(paste("FEAST_R2_averages = ", round(emnoise_R2_averages, 4)))
      print("End FEAST")
      
      #############EM#################
      
      if(ST_simulations_flag == 1){
        
        print("Start ST")
        
        #############ST START############
        envs_simulation_st = LETTERS[c(1:num_sources)]
        st <- sourcetracker(totalnew_source_old, envs_simulation_st, rarefaction_depth=COVERAGE)
        
        result<- foreach(i=1:st_mixing_iterations, .combine='c', .export=ls(.GlobalEnv)) %dopar% {
          sink=t(as.matrix(sinks[i,]))
          st_start_time<-as.numeric(proc.time()[3])
          results <- predict(st,sink, alpha1=alpha1, alpha2=alpha2)
          predict <- c(results$proportions) #predicted proportions in sample
          c(list(predict),list(as.numeric(proc.time()[3])-st_start_time))
        }
        
        
        st_rts<-split(result, 1:2)[[2]]
        st_rts<-unlist(st_rts)
        st_runtimes<-c(st_runtimes, mean(st_rts))
        st_hold=split(result, 1:2)[[1]]
        st_hold=unlist(st_hold)
        #####st_predictions=split(result, 1:(num_sources+1))dont uncomment
        st_predictions <- matrix(st_hold, nrow = st_mixing_iterations, ncol = num_sources+1, byrow = T)
        st_R2s<-c()
        if(include_epsilon==F){
          st_predictions=st_predictions[,1:num_sources]
        }
        for(j in 1:dim(st_predictions)[2]){
          st_R2s[j] = cor(st_predictions[c(1:st_mixing_iterations),j], ms_mat[c(1:st_mixing_iterations),j])^2
        }
        st_R2_averages = append(st_R2_averages,mean(na.omit(st_R2s)))
        newstmat<-as.matrix(mapply(rbind, st_predictions))
        st_m_matrices<-append(st_m_matrices, list(newstmat))
        st_sources_R2s<-append(st_sources_R2s, list(as.matrix(t(st_R2s))))
        
        ############ST###########
        print("End ST")
        
      }
    }
    
  }
  
  if(length(js_values) > 0) {
    
    if(JS > (js_values[t - 1] + 0.07) & JS < (js_values[t - 1] + 0.1)){
      
      
      
      t = t+1
      js_values = append(js_values, JS) 
      
      
      #############CLS###########
      
      
      #predict_values_cls <- c(); predict_eps_cls<-c(); cls_runtime<-c()
      #cls_start_time<-as.numeric(proc.time()[3])
      #bootsource<-matrix(NA,ncol = dim(sources[[1]])[1])
      #if(bs==T){ ###if we want to use bootstrapping 
      #  envs_simulation<-c()
      #  for(i in 1:length(sources)){
      #    bootsource=rbind(bootsource,t(as.matrix(sources[[i]])))
      #    bootsource=rbind(bootsource, t(
      #      rmultinom(n = 99, size = sum(totalnew_source[i,]), prob=totalnew_dist[i,])))
      #    envs_simulation=c(envs_simulation, rep(i, 100))
      #  }
      #  totalnew_source<-bootsource[complete.cases(bootsource),]
      #  rownames(totalnew_source)<-make.names(rownames(totalnew_source), unique = T)
      #  envs_simulation=factor(envs_simulation)
      #}
      
      new_source_old = new_source
      totalnew_source_old = totalnew_source
      
      
      
      
      ms_mat<-c()
      if(include_epsilon==F){
        ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      
      #############CLS END#############
      
      
      #############EM###########
      
      em_rts<-c(); emnoise_rts =c();
      predict_values_em<-c()
      predict_eps_em<-c()
      em_results<-c()
      emnoise_results<-c()
      predict_values_emnoise<-c()
      predict_eps_emnoise<-c()
      
      print("Start FEAST")
      
      for(i in 1:mixing_iterations){
        
        #Creating the unknown source per mixing iteration
        if(include_epsilon == TRUE){
          
          ##Adding the initial value of the unknown source for CLS and EM
          new_source_2 = list()
          totalnew_source_2 = matrix(NA, ncol = dim(totalnew_source_old)[2], nrow = ( dim(totalnew_source_old)[1] + 1))
          
          for(j in 1:n_sources){
            
            new_source_2[[j]] = new_source_old[[j]]
            totalnew_source_2[j,] = totalnew_source_old[j,]
          }
          
          #create unknown for each sink i
          unknown_source = unknown_initialize(sources = totalnew_source_old[c(1:n_sources),], sink = as.numeric(sinks[i,]), 
                                              n_sources = n_sources) 
          # new_source_2[[j+1]] = unknown_source
          # totalnew_source_2[(j+1),] = unknown_source
          # plot(as.numeric(unknown_source))
          
          # tmp_unk = change_C(newcov = COVERAGE, X = unknown_source)
          # tmp_unk =  rarefy(x = unknown_source, maxdepth = COVERAGE)
          # unknown_source = tmp_unk
          # new_source_2[[j+1]] = t(round(unknown_source))
          # source_2[[j+1]] = round(unknown_source)
          new_source_2[[j+1]] = round(unknown_source)
          totalnew_source_2[(j+1),] = round(unknown_source)
          
          
          totalnew_source = totalnew_source_2
          
          new_source=lapply(new_source_2,t)
          # totalnew_source <- t(apply(totalnew_source, 1, function(x) change_C(COVERAGE, x)))  ##COVERAGE CHANGE
          totalnew_source <- rarefy(x = totalnew_source, maxdepth = COVERAGE)
          new_source<- split(totalnew_source_2, seq(nrow(totalnew_source_2)))
          new_source<-lapply(new_source, as.matrix)
          
          # View(data.frame(t(totalnew_source), sinks[i,]))
          
          envs_simulation <- c(1:(num_sources+1))
          
          # View(data.frame(t(totalnew_source), sinks[i,]) )
        }
        
        samps <- list()
        samps <- new_source
        samps<-lapply(samps, t)
        
        # num_sources = 5
        observed_samps <- samps
        observed_samps[[(num_sources + 1)]] = t(rep(0,  dim(samps[[1]])[2]))
        
        # observed_samps[[(num_sources + 1)]] = t(rpois(lambda = 2,n = dim(samps[[1]])[2]))
        
        # observed_samps<-lapply(observed_samps, t)
        # str(observed_samps)
        
        # View(data.frame(t(totalsource),  sinks[1,]))
        
        
        # num_sources = 5
        if(eps==T) {initalphs<-runif(num_sources+1, 0.0, 1.0)
        }else {initalphs<-runif(num_sources, 0.0, 1.0)}
        initalphs=initalphs/Reduce("+", initalphs)
        m_s<-unlist(list_ms[[i]])
        # if(clsinit==T){initalphs<-c(predict_values_cls[i], 1-predict_values_cls[i])}
        sink=t(as.matrix(sinks[i,]))
        em_start_time<-as.numeric(proc.time()[3])
        pred_em<-do_EM_basic(alphas=initalphs, sources=samps, sink=sink, iterations=em_iterations)
        em_results<-c(em_results, pred_em)
        em_rts=c(em_rts, as.numeric(proc.time()[3])-em_start_time)
        emnoise_start_time<-as.numeric(proc.time()[3])
        tmp<-do_EM(alphas=initalphs, sources=samps, sink=sink, iterations=em_iterations, observed=observed_samps)
        pred_emnoise = tmp$toret
        emnoise_rts=c(emnoise_rts, as.numeric(proc.time()[3])-emnoise_start_time)
        emnoise_results=c(emnoise_results, pred_emnoise)
      }
      
      
      em_predictions=c()
      if(include_epsilon==F){
        em_predictions <- matrix(em_results, nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        em_predictions <- matrix(em_results, nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      em_runtimes<-c(em_runtimes, mean(em_rts))
      emnoise_runtimes<-c(emnoise_runtimes, mean(emnoise_rts))
      
      em_R2s<-c()
      for(j in 1:dim(em_predictions)[2]){
        em_R2s[j] = cor(em_predictions[,j], ms_mat[,j])^2
      }
      em_R2_averages<-c(em_R2_averages, mean(na.omit(em_R2s)))
      em_m_matrices<-append(em_m_matrices, list(mapply(rbind, em_predictions)))
      em_sources_R2s<-append(em_sources_R2s, list(as.matrix(t(em_R2s))))
      emnoise_predictions=c()
      if(include_epsilon==F){
        emnoise_predictions <- matrix(emnoise_results, nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        emnoise_predictions <- matrix(emnoise_results, nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      emnoise_R2s<-c()
      for(j in 1:dim(em_predictions)[2]){
        emnoise_R2s[j] = cor(emnoise_predictions[,j], ms_mat[,j])^2
      }
      emnoise_R2_averages<-c(emnoise_R2_averages, mean(na.omit(emnoise_R2s)))
      emnoise_m_matrices<-append(emnoise_m_matrices, list(mapply(rbind, emnoise_predictions)))
      emnoise_sources_R2s<-append(emnoise_sources_R2s, list(as.matrix(t(emnoise_R2s))))
      
      
    #   print(paste("FEAST_R2_averages = ", round(emnoise_R2_averages, 4)))
      print("End FEAST")
      
      #############EM#################
      
      if(ST_simulations_flag == 1){
        
        print("Start ST")
        
        #############ST START############
        envs_simulation_st = LETTERS[c(1:num_sources)]
        st <- sourcetracker(totalnew_source_old, envs_simulation_st, rarefaction_depth=COVERAGE)
        
        result<- foreach(i=1:st_mixing_iterations, .combine='c', .export=ls(.GlobalEnv)) %dopar% {
          sink=t(as.matrix(sinks[i,]))
          st_start_time<-as.numeric(proc.time()[3])
          results <- predict(st,sink, alpha1=alpha1, alpha2=alpha2)
          predict <- c(results$proportions) #predicted proportions in sample
          c(list(predict),list(as.numeric(proc.time()[3])-st_start_time))
        }
        
        
        st_rts<-split(result, 1:2)[[2]]
        st_rts<-unlist(st_rts)
        st_runtimes<-c(st_runtimes, mean(st_rts))
        st_hold=split(result, 1:2)[[1]]
        st_hold=unlist(st_hold)
        #####st_predictions=split(result, 1:(num_sources+1))dont uncomment
        st_predictions <- matrix(st_hold, nrow = st_mixing_iterations, ncol = num_sources+1, byrow = T)
        st_R2s<-c()
        if(include_epsilon==F){
          st_predictions=st_predictions[,1:num_sources]
        }
        for(j in 1:dim(st_predictions)[2]){
          st_R2s[j] = cor(st_predictions[c(1:st_mixing_iterations),j], ms_mat[c(1:st_mixing_iterations),j])^2
        }
        st_R2_averages = append(st_R2_averages,mean(na.omit(st_R2s)))
        newstmat<-as.matrix(mapply(rbind, st_predictions))
        st_m_matrices<-append(st_m_matrices, list(newstmat))
        st_sources_R2s<-append(st_sources_R2s, list(as.matrix(t(st_R2s))))
        
        ############ST###########
        print("End ST")
        
      }
    }
    
  }  
  
  
}

###EMP
source_data = readRDS("../data/60KFILT10KSOURCE.rds")
sum_sources = apply(source_data, 1, sum)
ind <- which(sum_sources >=10000)
source_data_new = source_data[ind,]

itr = 1000
start_EMP_flag = 1

for (it in 1:itr){
  envs_simulation <- 1:num_sources
  totalsource<-source_data[sample(1:dim(source_data)[1],num_sources,F),]
  totalsource<-as.matrix(totalsource)
  totalsource<- rarefy(x = totalsource, maxdepth = COVERAGE)
  
  sources <- split(totalsource, seq(nrow(totalsource)))
  sources<-lapply(sources, as.matrix)
  dists<-lapply(sources, function(x) x/(sum(colSums(x))))
  totaldist<-t(Reduce("cbind", dists))
  
  
  # print(str(totalsource))
  
  epsource<-rpois(n = length(totalsource[1,]), lambda = 10) #use a different source for the unknown 
  epdist<-epsource/(sum(epsource))
  
  # epsource<-source[sample(1:dim(source)[1],1,F),]
  # epdist<-epsource/(sum(colSums(epsource)))
  
  
  #Adding the unknown source
  if(include_epsilon == TRUE){
    
    totalsource_new = matrix(NA, nrow = (dim(totalsource)[1] + 1), ncol = (dim(totalsource)[2]))
    colnames(totalsource_new) = colnames(totalsource)
    rownames(totalsource_new) = seq(from = 1, to = (dim(totalsource)[1] + 1), by = 1)
    
    for(j in 1:(dim(totalsource)[1])){
      
      totalsource_new[j,] = totalsource[j,]
    }
    totalsource_new[(j+1),] = as.numeric(epsource)
  }
  
  totalsource = totalsource_new
  
  list_ms<-list()
  ep_values<-rep(0, mixing_iterations)
  
  
  ######change_C(X = totalnew_source[1,], newcov = COVERAGE) dont uncomment
  
  ##Known sources
  new_source = list()
  totalnew_source = matrix(NA, ncol = dim(totalsource)[2], nrow = length(sources))
  for(j in 1:length(sources)){
    new_source[[j]] = t(rmultinom(n = 1, size = sum(totalsource[j,]), prob=totaldist[j,]))
    totalnew_source[j,] = as.numeric(t(rmultinom(n = 1, size = sum(totalsource[j,]), prob=totaldist[j,])))}
  new_source=lapply(new_source,t)
  # totalnew_source=t(apply(totalnew_source, 1, function(x) change_C(COVERAGE, x)))  ##COVERAGE CHANGE
  totalnew_source <- rarefy(x = totalnew_source, maxdepth = COVERAGE)
  new_source<- split(totalnew_source, seq(nrow(totalnew_source)))
  new_source<-lapply(new_source, as.matrix)
  
  
  newdists<-lapply(new_source, function(x) x/(sum(colSums(x))))
  totalnew_dist<-t(Reduce("cbind", newdists))
  list_ms<-list()
  sinks<-matrix(nrow = mixing_iterations, ncol = dim(totalsource)[2])
  for(k in 1:mixing_iterations){
    list_ms[[k]] = create_m(num_sources = num_sources, n = 3, EPSILON = max(0, set_ep_prop[k]))
    sinks[k,] = t(as.matrix(list_ms[[k]])) %*% totalsource
    sinks[k,] = as.numeric(t(rmultinom(n=1, size=COVERAGE, prob=sinks[k,]/sum(sinks[k,]))))
  }
  
  
  #The sink is created using the known and unknown sources (if include_epsilon == TRUE )
  sinks=round(sinks)
  # sinks=t(apply(sinks, 1, function(x) change_C(COVERAGE, x)))
  # sinks = rarefy(x = sinks, maxdepth = COVERAGE)
  m_matrices=append(m_matrices, list(mapply(rbind, list_ms)))
  ms_mat=c()
  if(include_epsilon==F){
    ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources, byrow = T)
  }else{
    ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
  }
  
  
  #Jensen-Shannon divergence calculation - OLD VERSION
  # weights<-rep(1/num_sources, num_sources)
  # JS = mult_JSD(as.matrix(totaldist), q=weights); print(cat('jsd: ', JS))
  # js_values = append(js_values, JS)
  
  # Jensen-Shannon Divergence between P and Q - new version
  #x <- totalsource[c(1:num_sources),]
  JSDMatrix <- jsdmatrix(totalnew_source[c(1:num_sources),]) #instead of x
  JSDMatrix <- JSDMatrix/COVERAGE
  ##not good because the coverage is a huge range...
  #JSDMatrix<-sapply(1:num_sources, function(z) JSDMatrix[z,]/sum(totalsource[z,]))
  JS = mean(JSDMatrix[-which(JSDMatrix == 0)])
  print(JS)
  ###save the JSD values
  # plot a heatmap of the corresponding JSD matrix
  # heatmap(JSDMatrix)
  
  if(start_EMP_flag == 1) {
    
    if(JS < 0.7){
      
      
      start_EMP_flag = 0
      t = t+1
      js_values = append(js_values, JS) 
      
      
      new_source_old = new_source
      totalnew_source_old = totalnew_source
      
      
      
      
      ms_mat<-c()
      if(include_epsilon==F){
        ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      
      #############CLS END#############
      
      
      #############EM###########
      
      em_rts<-c(); emnoise_rts =c();
      predict_values_em<-c()
      predict_eps_em<-c()
      em_results<-c()
      emnoise_results<-c()
      predict_values_emnoise<-c()
      predict_eps_emnoise<-c()
      
      
      print("Start FEAST")
      
      for(i in 1:mixing_iterations){
        
        #Creating the unknown source per mixing iteration
        if(include_epsilon == TRUE){
          
          ##Adding the initial value of the unknown source for CLS and EM
          new_source_2 = list()
          totalnew_source_2 = matrix(NA, ncol = dim(totalnew_source_old)[2], nrow = ( dim(totalnew_source_old)[1] + 1))
          
          for(j in 1:n_sources){
            
            new_source_2[[j]] = new_source_old[[j]]
            totalnew_source_2[j,] = totalnew_source_old[j,]
          }
          
          #create unknown for each sink i
          unknown_source = unknown_initialize(sources = totalnew_source_old[c(1:n_sources),], sink = as.numeric(sinks[i,]), 
                                              n_sources = n_sources) 
          # new_source_2[[j+1]] = unknown_source
          # totalnew_source_2[(j+1),] = unknown_source
          # plot(as.numeric(unknown_source))
          
          # tmp_unk = change_C(newcov = COVERAGE, X = unknown_source)
          # tmp_unk =  rarefy(x = unknown_source, maxdepth = COVERAGE)
          # unknown_source = tmp_unk
          # new_source_2[[j+1]] = t(round(unknown_source))
          # source_2[[j+1]] = round(unknown_source)
          new_source_2[[j+1]] = round(unknown_source)
          totalnew_source_2[(j+1),] = round(unknown_source)
          
          
          totalnew_source = totalnew_source_2
          
          new_source=lapply(new_source_2,t)
          # totalnew_source <- t(apply(totalnew_source, 1, function(x) change_C(COVERAGE, x)))  ##COVERAGE CHANGE
          totalnew_source <- rarefy(x = totalnew_source, maxdepth = COVERAGE)
          new_source<- split(totalnew_source_2, seq(nrow(totalnew_source_2)))
          new_source<-lapply(new_source, as.matrix)
          
          # View(data.frame(t(totalnew_source), sinks[i,]))
          
          envs_simulation <- c(1:(num_sources+1))
          
          # View(data.frame(t(totalnew_source), sinks[i,]) )
        }
        
        samps <- list()
        samps <- new_source
        samps<-lapply(samps, t)
        
        # num_sources = 5
        observed_samps <- samps
        observed_samps[[(num_sources + 1)]] = t(rep(0,  dim(samps[[1]])[2]))
        
        # observed_samps[[(num_sources + 1)]] = t(rpois(lambda = 2,n = dim(samps[[1]])[2]))
        
        # observed_samps<-lapply(observed_samps, t)
        # str(observed_samps)
        
        # View(data.frame(t(totalsource),  sinks[1,]))
        
        
        # num_sources = 5
        if(eps==T) {initalphs<-runif(num_sources+1, 0.0, 1.0)
        }else {initalphs<-runif(num_sources, 0.0, 1.0)}
        initalphs=initalphs/Reduce("+", initalphs)
        m_s<-unlist(list_ms[[i]])
        # if(clsinit==T){initalphs<-c(predict_values_cls[i], 1-predict_values_cls[i])}
        sink=t(as.matrix(sinks[i,]))
        em_start_time<-as.numeric(proc.time()[3])
        pred_em<-do_EM_basic(alphas=initalphs, sources=samps, sink=sink, iterations=em_iterations)
        em_results<-c(em_results, pred_em)
        em_rts=c(em_rts, as.numeric(proc.time()[3])-em_start_time)
        emnoise_start_time<-as.numeric(proc.time()[3])
        tmp<-do_EM(alphas=initalphs, sources=samps, sink=sink, iterations=em_iterations, observed=observed_samps)
        pred_emnoise = tmp$toret
        emnoise_rts=c(emnoise_rts, as.numeric(proc.time()[3])-emnoise_start_time)
        emnoise_results=c(emnoise_results, pred_emnoise)
      }
      
      
      em_predictions=c()
      if(include_epsilon==F){
        em_predictions <- matrix(em_results, nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        em_predictions <- matrix(em_results, nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      em_runtimes<-c(em_runtimes, mean(em_rts))
      emnoise_runtimes<-c(emnoise_runtimes, mean(emnoise_rts))
      
      em_R2s<-c()
      for(j in 1:dim(em_predictions)[2]){
        em_R2s[j] = cor(em_predictions[,j], ms_mat[,j])^2
      }
      em_R2_averages<-c(em_R2_averages, mean(na.omit(em_R2s)))
      em_m_matrices<-append(em_m_matrices, list(mapply(rbind, em_predictions)))
      em_sources_R2s<-append(em_sources_R2s, list(as.matrix(t(em_R2s))))
      emnoise_predictions=c()
      if(include_epsilon==F){
        emnoise_predictions <- matrix(emnoise_results, nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        emnoise_predictions <- matrix(emnoise_results, nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      emnoise_R2s<-c()
      for(j in 1:dim(em_predictions)[2]){
        emnoise_R2s[j] = cor(emnoise_predictions[,j], ms_mat[,j])^2
      }
      emnoise_R2_averages<-c(emnoise_R2_averages, mean(na.omit(emnoise_R2s)))
      emnoise_m_matrices<-append(emnoise_m_matrices, list(mapply(rbind, emnoise_predictions)))
      emnoise_sources_R2s<-append(emnoise_sources_R2s, list(as.matrix(t(emnoise_R2s))))
      
      
    #   print(paste("FEAST_R2_averages = ", round(emnoise_R2_averages, 4)))
      print("End FEAST")
      
      #############EM#################
      
      if(ST_simulations_flag == 1){
        
        print("Start ST")
        
        #############ST START############
        envs_simulation_st = LETTERS[c(1:num_sources)]
        st <- sourcetracker(totalnew_source_old, envs_simulation_st, rarefaction_depth=COVERAGE)
        
        result<- foreach(i=1:st_mixing_iterations, .combine='c', .export=ls(.GlobalEnv)) %dopar% {
          sink=t(as.matrix(sinks[i,]))
          st_start_time<-as.numeric(proc.time()[3])
          results <- predict(st,sink, alpha1=alpha1, alpha2=alpha2)
          predict <- c(results$proportions) #predicted proportions in sample
          c(list(predict),list(as.numeric(proc.time()[3])-st_start_time))
        }
        
        
        st_rts<-split(result, 1:2)[[2]]
        st_rts<-unlist(st_rts)
        st_runtimes<-c(st_runtimes, mean(st_rts))
        st_hold=split(result, 1:2)[[1]]
        st_hold=unlist(st_hold)
        #####st_predictions=split(result, 1:(num_sources+1))dont uncomment
        st_predictions <- matrix(st_hold, nrow = st_mixing_iterations, ncol = num_sources+1, byrow = T)
        st_R2s<-c()
        if(include_epsilon==F){
          st_predictions=st_predictions[,1:num_sources]
        }
        for(j in 1:dim(st_predictions)[2]){
          st_R2s[j] = cor(st_predictions[c(1:st_mixing_iterations),j], ms_mat[c(1:st_mixing_iterations),j])^2
        }
        st_R2_averages = append(st_R2_averages,mean(na.omit(st_R2s)))
        newstmat<-as.matrix(mapply(rbind, st_predictions))
        st_m_matrices<-append(st_m_matrices, list(newstmat))
        st_sources_R2s<-append(st_sources_R2s, list(as.matrix(t(st_R2s))))
        
        ############ST###########
        print("End ST")
        
      }
    }
    
  }
  
  if(length(js_values) > 0) {
    
    if(JS > (js_values[t - 1] + 0.07) & JS < (js_values[t - 1] + 0.1)){
      
      
      t = t+1
      js_values = append(js_values, JS) 
      
      new_source_old = new_source
      totalnew_source_old = totalnew_source
      
      
      ms_mat<-c()
      if(include_epsilon==F){
        ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      
      #############CLS END#############
      
      
      #############EM###########
      
      em_rts<-c(); emnoise_rts =c();
      predict_values_em<-c()
      predict_eps_em<-c()
      em_results<-c()
      emnoise_results<-c()
      predict_values_emnoise<-c()
      predict_eps_emnoise<-c()
      
      print("Start FEAST")
      
      for(i in 1:mixing_iterations){
        
        #Creating the unknown source per mixing iteration
        if(include_epsilon == TRUE){
          
          ##Adding the initial value of the unknown source for CLS and EM
          new_source_2 = list()
          totalnew_source_2 = matrix(NA, ncol = dim(totalnew_source_old)[2], nrow = ( dim(totalnew_source_old)[1] + 1))
          
          for(j in 1:n_sources){
            
            new_source_2[[j]] = new_source_old[[j]]
            totalnew_source_2[j,] = totalnew_source_old[j,]
          }
          
          #create unknown for each sink i
          unknown_source = unknown_initialize(sources = totalnew_source_old[c(1:n_sources),], sink = as.numeric(sinks[i,]), 
                                              n_sources = n_sources) 
          # new_source_2[[j+1]] = unknown_source
          # totalnew_source_2[(j+1),] = unknown_source
          # plot(as.numeric(unknown_source))
          
          # tmp_unk = change_C(newcov = COVERAGE, X = unknown_source)
          # tmp_unk =  rarefy(x = unknown_source, maxdepth = COVERAGE)
          # unknown_source = tmp_unk
          # new_source_2[[j+1]] = t(round(unknown_source))
          # source_2[[j+1]] = round(unknown_source)
          new_source_2[[j+1]] = round(unknown_source)
          totalnew_source_2[(j+1),] = round(unknown_source)
          
          
          totalnew_source = totalnew_source_2
          
          new_source=lapply(new_source_2,t)
          # totalnew_source <- t(apply(totalnew_source, 1, function(x) change_C(COVERAGE, x)))  ##COVERAGE CHANGE
          totalnew_source <- rarefy(x = totalnew_source, maxdepth = COVERAGE)
          new_source<- split(totalnew_source_2, seq(nrow(totalnew_source_2)))
          new_source<-lapply(new_source, as.matrix)
          
          # View(data.frame(t(totalnew_source), sinks[i,]))
          
          envs_simulation <- c(1:(num_sources+1))
          
          # View(data.frame(t(totalnew_source), sinks[i,]) )
        }
        
        samps <- list()
        samps <- new_source
        samps<-lapply(samps, t)
        
        # num_sources = 5
        observed_samps <- samps
        observed_samps[[(num_sources + 1)]] = t(rep(0,  dim(samps[[1]])[2]))
        
        # observed_samps[[(num_sources + 1)]] = t(rpois(lambda = 2,n = dim(samps[[1]])[2]))
        
        # observed_samps<-lapply(observed_samps, t)
        # str(observed_samps)
        
        # View(data.frame(t(totalsource),  sinks[1,]))
        
        
        # num_sources = 5
        if(eps==T) {initalphs<-runif(num_sources+1, 0.0, 1.0)
        }else {initalphs<-runif(num_sources, 0.0, 1.0)}
        initalphs=initalphs/Reduce("+", initalphs)
        m_s<-unlist(list_ms[[i]])
        # if(clsinit==T){initalphs<-c(predict_values_cls[i], 1-predict_values_cls[i])}
        sink=t(as.matrix(sinks[i,]))
        em_start_time<-as.numeric(proc.time()[3])
        pred_em<-do_EM_basic(alphas=initalphs, sources=samps, sink=sink, iterations=em_iterations)
        em_results<-c(em_results, pred_em)
        em_rts=c(em_rts, as.numeric(proc.time()[3])-em_start_time)
        emnoise_start_time<-as.numeric(proc.time()[3])
        tmp<-do_EM(alphas=initalphs, sources=samps, sink=sink, iterations=em_iterations, observed=observed_samps)
        pred_emnoise = tmp$toret
        emnoise_rts=c(emnoise_rts, as.numeric(proc.time()[3])-emnoise_start_time)
        emnoise_results=c(emnoise_results, pred_emnoise)
      }
      
      
      em_predictions=c()
      if(include_epsilon==F){
        em_predictions <- matrix(em_results, nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        em_predictions <- matrix(em_results, nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      em_runtimes<-c(em_runtimes, mean(em_rts))
      emnoise_runtimes<-c(emnoise_runtimes, mean(emnoise_rts))
      
      em_R2s<-c()
      for(j in 1:dim(em_predictions)[2]){
        em_R2s[j] = cor(em_predictions[,j], ms_mat[,j])^2
      }
      em_R2_averages<-c(em_R2_averages, mean(na.omit(em_R2s)))
      em_m_matrices<-append(em_m_matrices, list(mapply(rbind, em_predictions)))
      em_sources_R2s<-append(em_sources_R2s, list(as.matrix(t(em_R2s))))
      emnoise_predictions=c()
      if(include_epsilon==F){
        emnoise_predictions <- matrix(emnoise_results, nrow = mixing_iterations, ncol = num_sources, byrow = T)
      }else{
        emnoise_predictions <- matrix(emnoise_results, nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
      }
      emnoise_R2s<-c()
      for(j in 1:dim(em_predictions)[2]){
        emnoise_R2s[j] = cor(emnoise_predictions[,j], ms_mat[,j])^2
      }
      emnoise_R2_averages<-c(emnoise_R2_averages, mean(na.omit(emnoise_R2s)))
      emnoise_m_matrices<-append(emnoise_m_matrices, list(mapply(rbind, emnoise_predictions)))
      emnoise_sources_R2s<-append(emnoise_sources_R2s, list(as.matrix(t(emnoise_R2s))))
      
      
     
    #   print(paste("FEAST_R2_averages = ", round(emnoise_R2_averages, 4)))
      print("End FEAST")
      
      #############EM#################
      
      if(ST_simulations_flag == 1){
        
        print("Start ST")
        
        #############ST START############
        envs_simulation_st = LETTERS[c(1:num_sources)]
        st <- sourcetracker(totalnew_source_old, envs_simulation_st, rarefaction_depth=COVERAGE)
        
        result<- foreach(i=1:st_mixing_iterations, .combine='c', .export=ls(.GlobalEnv)) %dopar% {
          sink=t(as.matrix(sinks[i,]))
          st_start_time<-as.numeric(proc.time()[3])
          results <- predict(st,sink, alpha1=alpha1, alpha2=alpha2)
          predict <- c(results$proportions) #predicted proportions in sample
          c(list(predict),list(as.numeric(proc.time()[3])-st_start_time))
        }
        
        
        st_rts<-split(result, 1:2)[[2]]
        st_rts<-unlist(st_rts)
        st_runtimes<-c(st_runtimes, mean(st_rts))
        st_hold=split(result, 1:2)[[1]]
        st_hold=unlist(st_hold)
        #####st_predictions=split(result, 1:(num_sources+1))dont uncomment
        st_predictions <- matrix(st_hold, nrow = st_mixing_iterations, ncol = num_sources+1, byrow = T)
        st_R2s<-c()
        if(include_epsilon==F){
          st_predictions=st_predictions[,1:num_sources]
        }
        for(j in 1:dim(st_predictions)[2]){
          st_R2s[j] = cor(st_predictions[c(1:st_mixing_iterations),j], ms_mat[c(1:st_mixing_iterations),j])^2
        }
        st_R2_averages = append(st_R2_averages,mean(na.omit(st_R2s)))
        newstmat<-as.matrix(mapply(rbind, st_predictions))
        st_m_matrices<-append(st_m_matrices, list(newstmat))
        st_sources_R2s<-append(st_sources_R2s, list(as.matrix(t(st_R2s))))
        
        ############ST###########
        print("End ST")
        
      }
      
    }
    
  }  
  
  
}



noisestring<-"nonoise"
if(noise==T){
  noisestring<-"noise"
}
bootstring<-"nobs" ##"no bootstrapping"
if(bs==T){
  bootstring<-"bs"
}
clsstring<-"noclsinit"
if(clsinit==T){
  clsstring<-"clsinit"
}
epstring<-"noeps"
cov=paste0("FULLCOV")


emnoise_sd = c()
st_sd = c()
for(j in 1:length(emnoise_sources_R2s)){
  
  emnoise_sd[j] = sd(emnoise_sources_R2s[[j]])/sqrt(length(emnoise_sources_R2s[[j]]))
  if(ST_simulations_flag ==1)
    st_sd[j] = sd(st_sources_R2s[[j]])/sqrt(length(st_sources_R2s[[j]]))
  
}



if(ST_simulations_flag == 1){
  

  run_tables<-rbind(st_runtimes,emnoise_runtimes, js_values)
  rownames(run_tables) <- c("SourceTracker_run_times","FEAST_run_times", "jsd_values")
  print("Run time table")
  print(run_tables)
  write.table(run_tables, file = "runtimes.txt", sep = '\t', quote = F)
  
  plot_table_simulation<-data.frame(cbind(st_R2_averages, emnoise_R2_averages,js_values))
  
  
  if(plot == 1){
  
  comp_plot<-ggplot(plot_table_simulation, aes(js_values)) +
    geom_errorbar(aes(ymin=st_R2_averages-st_sd, ymax=st_R2_averages+st_sd),
                  width=.01) +
    geom_point(aes(y=st_R2_averages, colour="SourceTracker"), size = 1, shape = 16) +
    geom_line(aes(y=st_R2_averages, colour="SourceTracker")) +
    geom_point(aes(y=emnoise_R2_averages, colour="FEAST"), size = 1, shape = 16) +
    geom_errorbar(aes(ymin=emnoise_R2_averages-emnoise_sd, 
                      ymax=emnoise_R2_averages+emnoise_sd), 
                  width=.01) +
    geom_line(aes(y=emnoise_R2_averages, colour="FEAST")) +  theme_bw() +ylim (0,1) + xlim(0,1)
    comp_plot<- comp_plot + labs(colour="Method", y="R2 of proportion estimates", 
                               x="Jensen shannon divergence")
  ggsave(filename="../results/Simulation_plot.png", plot = comp_plot , dpi = 600, width = 8.75, height = 6.1, units = "in")

  }
}

if(ST_simulations_flag == 0){
  
  
  run_tables<-rbind(emnoise_runtimes, js_values)
  rownames(run_tables) <- c("FEAST_run_times", "jsd_values")
  print("Run time table")
  print(run_tables)
  write.table(run_tables, file="runtimes.txt", sep = '\t', quote = F)
  
  
  plot_table_simulation<-data.frame(cbind(emnoise_R2_averages,js_values))
  
  if(plot == 1){
  
  comp_plot<-ggplot(plot_table_simulation, aes(js_values)) +
    geom_point(aes(y=emnoise_R2_averages, colour="FEAST"), size = 1, shape = 16) +
    geom_line(aes(y=emnoise_R2_averages, colour="FEAST")) +  theme_bw() +ylim (0,1) + xlim(0,1)
  comp_plot<- comp_plot + labs(colour="Method", y="R2 of proportion estimates", 
                               x="Jensen shannon divergence")
                               
  ggsave(filename="../results/Simulation_plot.png", plot = comp_plot , dpi = 600, width = 8.75, height = 6.1, units = "in")

  }
  
}


