
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
  
  #Sequencing depth----
  # rm(list = ls())
  # gc()

source("SourceTracker.R")
source("src.R")
source_data<-readRDS("../data/60KFILT10KSOURCE.rds")
COVERAGE=c(seq(100,900,100), seq(1000,10000,1000))

include_epsilon<-T
n_sources<-20
mixing_iterations<- 30 #should be 100
st_mixing_iterations<- 30 #should be 100
epmean<-0.01
set_ep_prop<- abs(rnorm(n = mixing_iterations, sd = 0.01, mean = epmean) + 0.01) ##proportion of the unknown source
sample_iterations<- 10
m_iterations<-100 ###stop if m changes by 10^-6 or 1000 iterations
n_cores<-8
em_iterations <- 100
cl<- makeCluster(n_cores)#, outfile="logtest.txt")
registerDoParallel(cl)



epstring="NOEPS"
if(include_epsilon==T){
    epstring<-substr(paste0(epmean),start=nchar(paste0(epmean))-1, stop=nchar(paste0(epmean)))
}
unk=T; dat="60KFILT10KSOURCE"; itr=length(COVERAGE); noise=T;
unknown=""
eps=include_epsilon; bs=F; clsinit=T; em_itr=em_iterations; num_sources=n_sources
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


sum_sources = apply(source_data, 1, sum)
ind <- which(sum_sources >=10000)
source_data_new = source_data[ind,]


for (it in 1:itr){
  envs_simulation <- 1:num_sources
  #####adding support for multiple sources#####
  totalsource<-source_data_new[sample(1:dim(source_data_new)[1],num_sources,F),]
  totalsource<-as.matrix(totalsource)
  apply(totalsource, 1, sum)
  totalsource<- rarefy(x = totalsource, maxdepth = COVERAGE[it])
  apply(totalsource, 1, sum)
  
  sources <- split(totalsource, seq(nrow(totalsource)))
  sources<-lapply(sources, as.matrix)
  dists<-lapply(sources, function(x) x/(sum(colSums(x))))
  totaldist<-t(Reduce("cbind", dists))
  
  # print(str(totalsource))
  
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
  
  
  ##Known sources
  new_source = list()
  totalnew_source = matrix(NA, ncol = dim(totalsource)[2], nrow = length(sources))
  for(j in 1:length(sources)){
    new_source[[j]] = t(rmultinom(n = 1, size = sum(totalsource[j,]), prob=totaldist[j,]))
    totalnew_source[j,] = as.numeric(t(rmultinom(n = 1, size = sum(totalsource[j,]), prob=totaldist[j,])))}
  new_source=lapply(new_source,t)
  totalnew_source=t(apply(totalnew_source, 1, function(x) change_C(COVERAGE[it], x)))  ##COVERAGE[it] CHANGE
  # apply(totalnew_source, 1, sum)
  new_source<- split(totalnew_source, seq(nrow(totalnew_source)))
  new_source<-lapply(new_source, as.matrix)
  
  
  newdists<-lapply(new_source, function(x) x/(sum(colSums(x))))
  totalnew_dist<-t(Reduce("cbind", newdists))
  list_ms<-list()
  sinks<-matrix(nrow = mixing_iterations, ncol = dim(totalsource)[2])
  for(k in 1:mixing_iterations){
    list_ms[[k]] = create_m(num_sources = num_sources, n = 3, EPSILON = max(0, set_ep_prop[k]))
    sinks[k,] = t(as.matrix(list_ms[[k]])) %*% totalsource
    sinks[k,] = as.numeric(t(rmultinom(n=1, size=COVERAGE[it], prob=sinks[k,]/sum(sinks[k,]))))
  }
  
  
  #The sink is created using the known and unknown sources (if include_epsilon == TRUE )
  sinks=round(sinks)
  sinks=t(apply(sinks, 1, function(x) change_C(COVERAGE[it], x)))
  m_matrices=append(m_matrices, list(mapply(rbind, list_ms)))
  ms_mat=c()
  if(include_epsilon==F){
    ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources, byrow = T)
  }else{
    ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
  }
  
  

  
  # Jensen-Shannon Divergence between P and Q - new version
  JSDMatrix <- jsdmatrix(totalsource[c(1:num_sources),]) #instead of x
  JSDMatrix <- JSDMatrix/COVERAGE[it]
  JS = mean(JSDMatrix[-which(JSDMatrix == 0)])
  js_values = append(js_values, JS) 
 

  print(paste0( "COVERAGE = ",COVERAGE[it]) )   
  
  #############CLS###########
  new_source_old = new_source
  totalnew_source_old = totalnew_source
  
  predict_values_cls <- c(); predict_eps_cls<-c(); cls_runtime<-c()
  cls_start_time<-as.numeric(proc.time()[3])
  bootsource<-matrix(NA,ncol = dim(sources[[1]])[1])
  if(bs==T){ ###if we want to use bootstrapping 
    envs_simulation<-c()
    for(i in 1:length(sources)){
      bootsource=rbind(bootsource,t(as.matrix(sources[[i]])))
      bootsource=rbind(bootsource, t(
        rmultinom(n = 99, size = sum(totalnew_source[i,]), prob=totalnew_dist[i,])))
      envs_simulation=c(envs_simulation, rep(i, 100))
    }
    totalnew_source<-bootsource[complete.cases(bootsource),]
    rownames(totalnew_source)<-make.names(rownames(totalnew_source), unique = T)
    envs_simulation=factor(envs_simulation)
  }
  
  
  
  cls_result<-foreach(i=1:mixing_iterations, .combine='c', .export=ls(.GlobalEnv), .packages = "lsei") %dopar% {
    
    #Creating the unknown source per mixing iteration
    if(include_epsilon == TRUE){
      
      ##Adding the initial value of the unknown source for CLS and EM
      new_source_2 = list()
      totalnew_source_2 = matrix(NA, ncol = dim(totalnew_source_old)[2], nrow = ( dim(totalnew_source_old)[1] + 1))
      
      for(j in 1:length(sources)){
        
        new_source_2[[j]] = new_source_old[[j]]
        totalnew_source_2[j,] = totalnew_source_old[j,]
      }
      
      #create unknown for each sink i
      unknown_source = unknown_initialize(sources = totalnew_source_old[c(1:n_sources),], sink = as.numeric(sinks[i,]), 
                                          n_sources = n_sources) 
      # new_source_2[[j+1]] = unknown_source
      # totalnew_source_2[(j+1),] = unknown_source
      # plot(as.numeric(unknown_source))
      
      tmp_unk = rarefy(maxdepth = COVERAGE[it], x = as.matrix(unknown_source))
      # tmp_unk = change_C(newcov = COVERAGE[it], X = unknown_source)
      unknown_source = tmp_unk
      new_source_2[[j+1]] = t(round(unknown_source))
      # source_2[[j+1]] = round(unknown_source)
      totalnew_source_2[(j+1),] = round(unknown_source)
      
      
      totalnew_source = totalnew_source_2
      
      new_source=lapply(new_source_2,t)
      # totalsource=t(apply(totalsource, 1, function(x) change_C(COVERAGE[it], x)))  ##COVERAGE[it] CHANGE
      new_source<- split(totalnew_source_2, seq(nrow(totalnew_source_2)))
      new_source<-lapply(new_source, as.matrix)
      
      
      envs_simulation <- c(1:(num_sources+1))
      
      # View(data.frame(t(totalnew_source), sinks[i,]) )
    }
    
    # i= 1
    X <- source_process_nounknown(totalnew_source, envs_simulation, rarefaction_depth = COVERAGE[it])
    X <- t(t(X)/rowSums(t(X)))
    y<-t(as.matrix(sinks[i,]))
    #relative abundance
    y <- y/rowSums(y)
    X_mat = as.matrix(X)
    y_mat = as.matrix(y)
    cls_start_time<-as.numeric(proc.time()[3])
    test = lsei(a = X_mat, b = t(y_mat), e =  diag(length(envs_simulation)), f = rep(0, length(envs_simulation)))
    #test = test/sum(test)
    
    # test = lsei(a = X_mat, b = as.numeric(y_mat), e =  diag(length(envs_simulation)), 
    #             f = rep(0, length(envs_simulation)), c = t(rep(1, length(envs_simulation))),d = 1)
    test[test < 0] = 0
    
    cls_run_time<- as.numeric(proc.time()[3]) - cls_start_time
    predict_cls = test
    subtractsrc<-0.0
    for(j in 1:num_sources){subtractsrc=subtractsrc + c(predict_cls)[j]*X[,j] }
    subtractsrc=as.matrix(subtractsrc)
    c = y - t(subtractsrc) #prediction error vector
    # if(include_epsilon==T){
    #   predict_cls<-c(predict_cls, 1-sum(predict_cls))
    # }
    predict_cls=predict_cls/sum(predict_cls)
    c(list(predict_cls), list(c),list(cls_run_time)) #c(predcls[1], predcls[3], list(c))
  }
  
  cls_rts<-split(cls_result, 1:3)[[3]]
  cls_rts=unlist(cls_rts)
  hold<-split(cls_result,1:3)[[1]] ##so c's as [[2]], mlists[[1]]
  hold<-unlist(hold)
  
  if(include_epsilon==F){
    cls_predictions <- matrix(hold, nrow = mixing_iterations, ncol = num_sources, byrow = T)
  }else{
    cls_predictions <- matrix(hold, nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
  }
  ms_mat<-c()
  if(include_epsilon==F){
    ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources, byrow = T)
  }else{
    ms_mat <-  matrix(unlist(list_ms), nrow = mixing_iterations, ncol = num_sources+1, byrow = T)
  }
  cls_R2s<-c()
  for(j in 1:dim(cls_predictions)[2]){
    cls_R2s[j] = cor(cls_predictions[,j], ms_mat[,j])^2
  }
  cls_m_matrices<-append(cls_m_matrices, list(mapply(rbind, cls_predictions)))
  cls_sources_R2s<-append(cls_sources_R2s, list(as.matrix(t(cls_R2s))))
  cls_R2_averages <- append(cls_R2_averages, mean(na.omit(cls_R2s)))
  cls_runtimes<-c(cls_runtimes, mean(cls_rts))
  
#   print(paste("cls_R2_averages = ", round(cls_R2_averages, 4)))
  
  
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
    
    print(i)
    #Creating the unknown source per mixing iteration
    if(include_epsilon == TRUE){
      
      ##Adding the initial value of the unknown source for CLS and EM
      new_source_2 = list()
      totalnew_source_2 = matrix(NA, ncol = dim(totalnew_source_old)[2], nrow = ( dim(totalnew_source_old)[1] + 1))
      
      for(j in 1:length(sources)){
        
        new_source_2[[j]] = new_source_old[[j]]
        totalnew_source_2[j,] = totalnew_source_old[j,]
      }
      
      #create unknown for each sink i
      unknown_source = unknown_initialize(sources = totalnew_source_old[c(1:n_sources),], sink = as.numeric(sinks[i,]), 
                                          n_sources = n_sources) 
      # new_source_2[[j+1]] = unknown_source
      # totalnew_source_2[(j+1),] = unknown_source
      # plot(as.numeric(unknown_source))
      
      tmp_unk = rarefy(maxdepth = COVERAGE[it], x = as.matrix(unknown_source))
      unknown_source = tmp_unk
      new_source_2[[j+1]] = t(round(unknown_source))
      # source_2[[j+1]] = round(unknown_source)
      totalnew_source_2[(j+1),] = round(unknown_source)
      
      
      totalnew_source = totalnew_source_2
      
      new_source<- lapply(new_source_2,t)
      # totalnew_source <- t(apply(totalsource, 1, function(x) change_C(COVERAGE[it], x)))  ##COVERAGE[it] CHANGE
      new_source<- split(totalnew_source_2, seq(nrow(totalnew_source_2)))
      new_source<-lapply(new_source, as.matrix)
      
      
      envs_simulation <- c(1:(num_sources+1))
      
      # View(data.frame(t(totalnew_source), sinks[i,]) )
    }
    
    samps <- list()
    samps <- new_source
    samps<-lapply(samps, t)
    
    # num_sources = 5
    observed_samps <- samps
    observed_samps[[(num_sources + 1)]] = t(rep(0,  dim(samps[[1]])[2]))
    # observed_samps<-lapply(observed_samps, t)
    # str(observed_samps)
    
    # View(data.frame(t(totalsource),  sinks[1,]))
    
    
    # num_sources = 5
    if(eps==T) {initalphs<-runif(num_sources+1, 0.0, 1.0)
    }else {initalphs<-runif(num_sources, 0.0, 1.0)}
    initalphs=initalphs/Reduce("+", initalphs)
    m_s<-unlist(list_ms[[i]])
    if(clsinit==T){initalphs<-as.numeric(cls_predictions[i,])}
    sink=t(as.matrix(sinks[i,]))
    em_start_time<-as.numeric(proc.time()[3])
    pred_em<-do_EM_basic(alphas=initalphs, sources=samps, sink=sink, iterations=em_itr)
    em_results<-c(em_results, pred_em)
    em_rts=c(em_rts, as.numeric(proc.time()[3])-em_start_time)
    emnoise_start_time<-as.numeric(proc.time()[3])
    tmp<-do_EM(alphas=initalphs, sources=samps, sink=sink, iterations=em_itr, observed=observed_samps)
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
  for(j in 1:dim(emnoise_predictions)[2]){
    emnoise_R2s[j] = cor(emnoise_predictions[,j], ms_mat[,j])^2
  }
  emnoise_R2_averages<-c(emnoise_R2_averages, mean(na.omit(emnoise_R2s)))
  emnoise_m_matrices<-append(emnoise_m_matrices, list(mapply(rbind, emnoise_predictions)))
  emnoise_sources_R2s<-append(emnoise_sources_R2s, list(as.matrix(t(emnoise_R2s))))
  
  print("End FEAST")
  print(paste("FEAST_R2_averages = ", round(emnoise_R2_averages, 4)))
#   print(paste("em_R2_averages = ", round(em_R2_averages, 4)))
  
  
  #############EM#################
  

  
if(st_flag_seq_depth == 1){
  
  print("Start ST")
  
  #############ST START############
  
  envs_simulation_st = LETTERS[c(1:num_sources)]
  
  st <- sourcetracker(totalnew_source_old, envs_simulation_st, rarefaction_depth=COVERAGE[it])

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
  ######st_predictions=split(result, 1:(num_sources+1))dont uncomment
  st_predictions <- matrix(st_hold, nrow = st_mixing_iterations, ncol = num_sources+1, byrow = T)
  st_R2s<-c()
  if(include_epsilon==F){
    st_predictions=st_predictions
  }
  for(j in 1:dim(st_predictions)[2]){
    st_R2s[j] = cor(st_predictions[c(1:st_mixing_iterations),j], ms_mat[c(1:st_mixing_iterations),j])^2
  }
  st_R2_averages = append(st_R2_averages,mean(na.omit(st_R2s)))
  newstmat<-as.matrix(mapply(rbind, st_predictions))
  st_m_matrices<-append(st_m_matrices, list(newstmat))
  st_sources_R2s<-append(st_sources_R2s, list(as.matrix(t(st_R2s))))

  print(paste("st_R2_averages = ", round(st_R2_averages, 4)))
  
  
  print("End ST")
  #############END ST###########
    
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
  cov=paste0("COVERAGE[it]")
  

  
  
if(st_flag_seq_depth == 1){
  
  #R2_table<-rbind(st_R2_averages,cls_R2_averages, em_R2_averages, emnoise_R2_averages, js_values)
  R2_table<-rbind(st_R2_averages, emnoise_R2_averages, COVERAGE)
  rownames(R2_table) <- c("SourceTracker_R2","FEAST_R2", "COVERAGE")
  print("R2 table")
  print(R2_table)
  
  run_tables<-rbind(st_runtimes,emnoise_runtimes, COVERAGE)
  rownames(run_tables) <- c("SourceTracker_run_times","FEAST_run_times", "COVERAGE")
  print("Run time table")
  print(run_tables)
  write.table(R2_table, file = "R2values.txt", sep = '\t', quote = F)
  write.table(run_tables, file = "runtimes.txt", sep = '\t', quote = F)
  
  plot_table<-data.frame(cbind(st_R2_averages, emnoise_R2_averages,COVERAGE))
  
  comp_plot<-(ggplot(plot_table, aes(COVERAGE)) +
                geom_line(aes(y=st_R2_averages, colour="ST_R2")) +
                # geom_line(aes(y=cls_R2_averages, colour="CLS_R2")) + 
                geom_line(aes(y=emnoise_R2_averages, colour="FEAST")) + theme_bw())
  comp_plot<- comp_plot + labs(colour="Method", y="Average R2 (Source proportions)", x="JSD Value")
  
}

if(st_flag_seq_depth == 0){
  
  R2_table<-rbind(emnoise_R2_averages, COVERAGE)
  rownames(R2_table) <- c("FEAST_R2", "COVERAGE")
  print("R2 table")
  print(R2_table)
  
  run_tables<-rbind(emnoise_runtimes, COVERAGE)
  rownames(run_tables) <- c("FEAST_run_times", "COVERAGE")
  print("Run time table")
  print(run_tables)
  write.table(R2_table, file="R2values.txt", sep = '\t', quote = F)
  write.table(run_tables, file="runtimes.txt", sep = '\t', quote = F)
                          
  plot_table<-data.frame(cbind(emnoise_R2_averages,COVERAGE))
                          
  comp_plot<-(ggplot(plot_table, aes(COVERAGE)) +
              geom_line(aes(y=emnoise_R2_averages, colour="FEAST")) + theme_bw()) + 
              labs(colour="Method", y="Average R2 (Source proportions)", x="JS Value")
                          
}

ggsave(filename="../results/Seq_depth_plot.png", plot = comp_plot , dpi = 600, width = 8.75, height = 6.1, units = "in")
  
  
#   stopCluster(cl)

  
