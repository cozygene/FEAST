x<-c("vegan", "dplyr", "ggrepel", "doParallel", "foreach" ,"mgcv", "reshape2", "ggplot2", "Rcpp", "RcppArmadillo")
lapply(x, require, character.only = TRUE)
library("vegan")
library("dplyr")
library("doParallel")
library("foreach")
library("mgcv")
library("reshape2")
library("ggplot2")
library("cowplot")
library("Rcpp")
library("RcppArmadillo")
cppFunction("arma::mat schur(arma::mat& a, arma::mat& b) 
            {return(a % b); }", depends="RcppArmadillo")


"change_C"<-function(newcov, X){

  
  X=t(as.matrix(X))
  idx = 1:dim(X)[2]
  
  if(sum(X) > newcov){
    
    while(sum(X) > newcov){
      greaterone = X > 1
      samps = 20
      if(samps > length(X[greaterone]))
        samps = length(X[greaterone])
      changeidx = sample(idx[greaterone], samps, replace = F)
      X[changeidx] = X[changeidx] - 1
    }
    
  }

  if(sum(X) < newcov){
    
    while(sum(X) < newcov){
      greaterone = X > 1
      samps = 100
      if(samps > length(X[greaterone]))
        samps = length(X[greaterone])
      changeidx = sample(idx[greaterone], samps, replace = F)
      X[changeidx] = X[changeidx] + 1
    }
    
  }
  
  return(X)
}

rarefy <- function(x,maxdepth){
 
  
  if(is.null(maxdepth)) return(x)
  
  if(!is.element(class(x), c('matrix', 'data.frame','array')))
    x <- matrix(x,nrow=nrow(x))
  nr <- nrow(x)
  nc <- ncol(x)
  
  for(i in 1:nrow(x)){
    if(sum(x[i,]) > maxdepth){
      prev.warn <- options()$warn
      options(warn=-1)
      s <- sample(nc, size=maxdepth, prob=x[i,], replace=T)
      options(warn=prev.warn)
      x[i,] <- hist(s,breaks=seq(.5,nc+.5,1), plot=FALSE)$counts
    }
  }
  return(x)
}

"jsdmatrix" <- function(x){
  d <- matrix(0,nrow=nrow(x),ncol=nrow(x))
  for(i in 1:(nrow(x)-1)){
    for(j in (i+1):nrow(x)){
      d[i,j] <- jsd(x[i,], x[j,])
      d[j,i] <- d[i,j]
    }
  }
  return(d)
}

"jsd" <- function(p,q){
  m <- (p + q)/2
  return((kld(p,m) + kld(q,m))/2)
}


"h"<-function(x) {y <- x[x > 0]; -sum(y * log(y))}; 
"mult_JSD" <- function(p,q) {h(q %*% p) - q %*% apply(p, 1, h)}

"retrands"<-function(V){
  toret<-unlist(lapply(c(V), function(x) runif(1, x+1e-12, x+1e-09)))
  return(toret)
}

"getR2"<-function(x,y){
  return((cor(x,y))^2)
}

"E"<-function(alphas, sources){
  nums<-(sapply(1:length(alphas), function(n) Reduce("+", crossprod(as.numeric(alphas[n]),as.numeric(sources[[n]])))))
  denom<-(Reduce("+", nums))
  return(nums/denom)
}

"A"<-function(alph, XO, raos){
  tmp<-crossprod(alph, XO/raos)
  tmp<-rapply(list(tmp), f=function(x) ifelse(is.nan(x),0,x), how="replace" )
  tmp<-Reduce("+",unlist(tmp))
  return(tmp)
}

"M"<-function(alphas, sources, sink, observed){
    
  newalphs<-c()
  rel_sink <-sink/sum(sink)
  
  if(sum(sources[[1]]) > 1){
    
    sources <-lapply(sources, function(x) x/(sum(colSums(x))))
  }
  
  
  LOs<-lapply(sources, schur, b=rel_sink)
  BOs<-t(mapply(crossprod, x=sources, y=alphas))
  BOs<-split(BOs, seq(nrow(BOs)))
  BOs<-lapply(BOs, as.matrix)
  BOs<-lapply(BOs, t)
  num_list <- list()
  source_new <- list()
  
  
  for(i in 1:length(sources)){
    num <- c()
    denom <- c()
    num<-crossprod(alphas[i], (LOs[[i]]/(Reduce("+", BOs))))
    num<-rapply(list(num), f=function(x) ifelse(is.nan(x),0,x), how="replace" ) #replace na with zero
    # num_list[[i]]<- num[[1]][1,] + sources[[i]][1,]
    num_list[[i]]<- num[[1]][1,] + observed[[i]][1,]
    
    denom <- Reduce("+",unlist(num_list[[i]]))
    source_new[[i]] <- num_list[[i]]/denom
    source_new[[i]][is.na(source_new[[i]])] = 0
  }
  
  sources = source_new
  
  newalphs<-c()
  sources<-lapply(sources, t)
  XOs<-lapply(sources,schur, b=rel_sink)
  AOs<-t(mapply(crossprod, x=sources, y=alphas))
  AOs<-split(AOs, seq(nrow(AOs)))
  AOs<-lapply(AOs, as.matrix)
  AOs<-lapply(AOs, t)
  newAs<-c()
  for(i in 1:length(sources)){
    newA<-crossprod(alphas[i], (XOs[[i]]/(Reduce("+", AOs))))
    newA<-rapply(list(newA), f=function(x) ifelse(is.nan(x),0,x), how="replace" )
    newA<-Reduce("+",unlist(newA))
    newAs<-c(newAs, newA)
  }
  tot<-sum(newAs)
  Results <- list (new_alpha = newAs/(tot), new_sources = sources)
  return(Results)
}

"do_EM"<-function(alphas, sources, observed, sink, iterations){
    
  curalphas<-alphas
  newalphas<-alphas
  m_guesses<-c(alphas[1])
  for(itr in 1:iterations){
    
    curalphas<-E(newalphas, sources)
    tmp <- M(alphas = curalphas, sources = sources, sink = sink, observed = observed)
    newalphas <- tmp$new_alpha
    sources <- tmp$new_sources  
    m_guesses<-c(m_guesses, newalphas[1])
    if(abs(m_guesses[length(m_guesses)]-m_guesses[length(m_guesses)-1])<=10^-6) break
    
  }        
  toret<-c(newalphas)
  results <- list(toret = toret, sources = sources)
  
  return(results)
}

"M_basic"<-function(alphas, sources, sink){
  newalphs<-c()
  XOs<-lapply(sources,schur, b=sink)
  AOs<-t(mapply(crossprod, x=sources, y=alphas))
  AOs<-split(AOs, seq(nrow(AOs)))
  AOs<-lapply(AOs, as.matrix)
  AOs<-lapply(AOs, t)
  newAs<-c()
  for(i in 1:length(sources)){
    newA<-crossprod(alphas[i], (XOs[[i]]/(Reduce("+", AOs))))
    newA<-rapply(list(newA), f=function(x) ifelse(is.nan(x),0,x), how="replace" )
    newA<-Reduce("+",unlist(newA))
    newAs<-c(newAs, newA)
  }
  tot<-sum(newAs)
  return(newAs/(tot))
}

"do_EM_basic"<-function(alphas, sources, sink, iterations){
  curalphas<-alphas
  newalphas<-alphas
  m_guesses<-c(alphas[1])
  for(itr in 1:iterations){
    curalphas<-E(newalphas, sources)
    newalphas<-M_basic(curalphas, sources, sink)
    m_guesses<-c(m_guesses, newalphas[1])
    if(abs(m_guesses[length(m_guesses)]-m_guesses[length(m_guesses)-1])<=10^-6) break
 
  }        
  toret<-c(newalphas)
  return(toret)
}

"source_process_nounknown" <- function(train, envs, rarefaction_depth=1000){
  
  train <- as.matrix(train)
  
  # enforce integer data
  if(sum(as.integer(train) != as.numeric(train)) > 0){
    stop('Data must be integral. Consider using "ceiling(datatable)" or ceiling(1000*datatable) to convert floating-point data to integers.')
  }
  envs <- factor(envs)
  train.envs <- sort(unique(levels(envs)))
  
  # rarefy samples above maxdepth if requested
  if(!is.null(rarefaction_depth) && rarefaction_depth > 0) train <- rarefy(train, rarefaction_depth)
  
  # get source environment counts
  # sources is nenvs X ntaxa
  X <- t(sapply(split(data.frame(train), envs), colSums)) 
  
  rownames(X) <- c(train.envs)
  X <- t(as.matrix(X))
  
  return(X)
}

"read_pseudo_data"<-function(dataset){
  path_to_data<-"../data/"
  if(dataset=="DA"){
    df<-read.table(paste0(path_to_data,"DA_99_T_d10000_date_nan.txt"), fill = NA)
    return(df[complete.cases(df),])
  }else if(dataset=="DB"){
    df<-read.table(paste0(path_to_data,"DB_99_T_d10000_date_nan.txt"), fill = NA)
    return(df[complete.cases(df),])
  }else if (dataset=="F4"){
    df<-read.table(paste0(path_to_data,"F4_99_T_d10000_date_nan.txt"), fill = NA)
    return(df[complete.cases(df),])
  }else{
    df<-read.table(paste0(path_to_data,"M3_99_T_d10000_date_nan.txt"), fill = NA)
    return(df[complete.cases(df),])}
}

create_m <- function(num_sources, n, EPSILON){
  
  
  if( n == 1 ){
    
    index = sample(c(1:num_sources), 1)
    m_1 = runif(min = 0.6, max = 0.9, n = 1)
    resid = 1-m_1
    other_ms = resid/(num_sources-1)
    m = rep(NA, num_sources)
    m[index] = c(m_1)
    m[is.na(m)] = other_ms
    
  }
  
  
  if( n == 2 ){
    
    index = sample(c(1:num_sources), 2)
    m_1 = runif(min = 0.1, max = 0.2, n = 1)
    m_2 = runif(min = 0.4, max = 0.5, n = 1)
    resid = 1-(m_1+m_2)
    other_ms = resid/(num_sources-2)
    m = rep(NA, num_sources)
    m[index] = c(m_1, m_2)
    m[is.na(m)] = other_ms
    
  }
  
  
  if( n == 3 ){
    
    index = sample(c(1:num_sources), 3)
    m_1 = runif(min = 0.1, max = 0.5, n = 1)
    m_2 = runif(min = 0.2, max = 0.25, n = 1)
    m_3 = runif(min = 0.1, max = 0.15, n = 1)
    resid = 1-(m_1+m_2+m_3)
    other_ms = runif(min = 0.001, max = resid/(num_sources-3), n = (num_sources-3))
    m = rep(NA, num_sources)
    m[index] = c(m_1, m_2, m_3)
    m[is.na(m)] = other_ms
    m = m/sum(m)
    
  }
  subsum = 0
  idx = 1:length(m)
  
  while ((subsum+0.001) < EPSILON){
    tosub = EPSILON - subsum
    tosub = tosub / (num_sources+1)
    mask = m > tosub
    m[mask] = m[mask] - tosub
    subsum = subsum + length(m[mask]) * tosub
    
  }
  m = c(m,(EPSILON))
  
  # sum(m)
  return(m)
  
}


unknown_initialize <- function(sources, sink, n_sources){
  
  
  unknown_source = rep(0, length(sink))
  sum_sources = apply(sources, 2, sum)
  
  unknown_source = c()
  
  for(j in 1:length(sum_sources)){
    
    unknown_source[j] = max(sink[j]-sum_sources[j], 0)
    
  }

  
  return(unknown_source)
  
}


unknown_initialize_1 <- function(sources, sink, n_sources){
  
  
  unknown_source = rep(0, length(sink))
  
  #zero all the OTUs with at least 1 known source
  sources_sum = apply(sources, 2 ,sum)

  
  
  # unknown_source = rep(0, length(sink))
  
  unknown_source = c()
  
  for(j in 1:length(sources_sum)){
    
    unknown_source[j] = max(sink[j]-sources_sum[j], 0)
    
  }
  
  #Select the cor OTUs
  ind_cor = list()
  ind_known_source_abun = c()
  ind_cor_all = which(sources[1,] > 0)
  
  counter = matrix(0, ncol = dim(sources)[2], nrow =  dim(sources)[1])
  
  
  for(j in 1:n_sources){
    
    ind_cor[[j]] = which(sources[j,] > 0)
    
    for(k in 1:length(sources[j,])){
      
      if(sources[j,k] > 0){
        
        counter[j,k] = counter[j,k]+1
      }
      
      
    }
    
  }
  
  OTU_present_absent = apply(counter, 2, sum)
  ind_cor_all = which(OTU_present_absent >= round(n_sources*0.8))
  
  if(length(ind_cor_all) > 1){
    
    cor_abundance = round(apply(sources[,ind_cor_all], 2, median)/2) #take the min abundnace of the 'cor'
    unknown_source[ind_cor_all] = cor_abundance
    
  }
  
  
  #keep the sink abundance where there is no known source
  ind_no_known_source_abun = which(sources_sum == 0)
  
  for(j in 1:length(ind_no_known_source_abun)){
    
    # unknown_source[ind_no_known_source_abun[j]] = max(runif(n = 1, min = 1, max = 100), sink[ind_no_known_source_abun[j]])
    unknown_source[ind_no_known_source_abun[j]] = max((sink[ind_no_known_source_abun[j]] - rpois(n = 1, lambda = 0.5)), 0)
    
  }
  
  
  
  return(unknown_source)
  
}

unknown__initialize_1 <- function(sources, sink, n_sources){
  
  # sources = totalsource[c(1:n_sources),]
  # sink = as.numeric(sinks[1,])
  
  unknown_source = rep(0, length(sink))
  
  #zero all the OTUs with at least 1 known source
  sources_sum = apply(sources, 2 ,sum)
  ind_known_source_abun = which(sources_sum > 0)
  unknown_source[ind_known_source_abun] = 0
  
  
  #Select the cor OTUs
  ind_cor = list()
  ind_known_source_abun = c()
  ind_cor_all = which(sources[1,] > 0)
  
  counter = matrix(0, ncol = dim(sources)[2], nrow =  dim(sources)[1])
  
  
  for(j in 1:n_sources){
    
    ind_cor[[j]] = which(sources[j,] > 0)
    
    for(k in 1:length(sources[j,])){
      
      if(sources[j,k] > 0){
        
        counter[j,k] = counter[j,k]+1
      }
      
      
    }
    
  }
  
  OTU_present_absent = apply(counter, 2, sum)
  ind_cor_all = which(OTU_present_absent >= round(n_sources*0.8))
  
  if(length(ind_cor_all) > 1){
    
    cor_abundance = apply(sources[,ind_cor_all], 2, median) #take the median abundnace of the 'cor'
    unknown_source[ind_cor_all] = cor_abundance
    
  }
  
  
  
  #keep the sink abundance where there is no known source
  ind_no_known_source_abun = which(sources_sum == 0)
  
  for(j in 1:length(ind_no_known_source_abun)){
    
    # unknown_source[ind_no_known_source_abun[j]] = max(runif(n = 1, min = 1, max = 100), sink[ind_no_known_source_abun[j]])
    unknown_source[ind_no_known_source_abun[j]] = max( round(sink[ind_no_known_source_abun[j]]+ rnorm(n = length(sink[ind_no_known_source_abun[j]]))), 0)
    
  }
  
  
  
  return(unknown_source)
  
}



EM_results<- function(source = sources_data, sinks = sinks, em_itr = 1000, env = rownames(sources_data), include_epsilon = T, COVERAGE,
                      unknown_initialize = 0){
  
    
  tmp = source
  test_zeros = apply(tmp, 1, sum)
  ind_to_use = as.numeric(which(test_zeros > 0))
  ind_zero = as.numeric(which(test_zeros == 0))
  
  source = tmp[ind_to_use,]
  sinks = sinks
  
  
  
  #####adding support for multiple sources#####
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
  
  #Creating the unknown source per mixing iteration
  if(include_epsilon == TRUE){
    
    ##Adding the initial value of the unknown source for CLS and EM
    source_2 = list()
    totalsource_2 = matrix(NA, ncol = dim(totalsource_old)[2], nrow = ( dim(totalsource_old)[1] + 1))
    
    for(j in 1:num_sources){
      
      source_2[[j]] = source_old[[j]]
      totalsource_2[j,] = totalsource_old[j,]
    }
    
    #create unknown for each sink i
    
    sinks_rarefy = rarefy(matrix(sinks, nrow = 1), maxdepth = apply(totalsource_old, 1, sum)[1]) #make
    
    
    # unknown_source = unknown_initialize(sources = totalsource[c(1:num_sources),], sink = as.numeric(sinks),
    #                                     n_sources = num_sources)
    if(unknown_initialize == 1)
      unknown_source_1 = unknown_initialize_1(sources = totalsource[c(1:num_sources),], sink = as.numeric(sinks),
                                          n_sources = num_sources)
    
    
    if(unknown_initialize == 0)
      unknown_source_1 = unknown_initialize(sources = totalsource[c(1:num_sources),], sink = as.numeric(sinks),
                                          n_sources = num_sources)

      unknown_source = unknown_source_1 + rpois(n = length(sinks), lambda = 0.5)
    # unknown_source = rpois(n = length(sinks_rarefy), lambda = 1)
    
    # unknown_source_rarefy = rarefy(matrix(unknown_source, nrow = 1), maxdepth = apply(totalsource_old, 1, sum)[1])
    unknown_source_rarefy = rarefy(matrix(unknown_source, nrow = 1), maxdepth = COVERAGE)
    source_2[[j+1]] = t(unknown_source_rarefy)
    totalsource_2[(j+1),] = t(unknown_source_rarefy)
    totalsource = totalsource_2
    
    source=lapply(source_2,t)
    # totalsource <- rarefy(x = totalsource, maxdepth = COVERAGE)  
    source<- split(totalsource, seq(nrow(totalsource_2)))
    source<-lapply(source_2, as.matrix)
    
    envs_simulation <- c(1:(num_sources+1))
    
  }
  
  
  samps <- source
  samps<-lapply(samps, t)
  
  observed_samps <- samps
  observed_samps[[(num_sources + 1)]] = t(rep(0, dim(samps[[1]])[2]))
  
  
  initalphs<-runif(num_sources+1, 0.0, 1.0)
  initalphs=initalphs/Reduce("+", initalphs)
  sink_em = as.matrix(sinks)
  pred_em<-do_EM_basic(alphas=initalphs, sources=samps, sink=sink_em, iterations=em_itr)
  
  tmp<-do_EM(alphas=initalphs, sources=samps, sink=sink_em, iterations=em_itr, observed=observed_samps)
  pred_emnoise = tmp$toret
  
  k = 1
  pred_emnoise_all = c()
  pred_em_all = c()
  
  for(j in 1:length(env)){
    
    if(j %in% ind_to_use){
      
      pred_emnoise_all[j] = pred_emnoise[k]
      pred_em_all[j] = pred_em[k]
      k = k+1
      
    }
    
    else{
      
      pred_emnoise_all[j] = 0
      pred_em_all[j] = 0
    }
    
  }
  
  pred_emnoise_all[j+1] = pred_emnoise[k]
  pred_em_all[j+1] = pred_em[k]
  
  
  
  names(pred_emnoise_all) = c(env,"unknown")
  names(pred_em_all) = c(env,"unknown")
  
  
  Results = list(unknown_source = unknown_source, unknown_source_rarefy = unknown_source_rarefy, 
                 data_prop = data.frame(pred_emnoise_all,pred_em_all))
  return(Results)
  
}


create_CI_EM <- function(totalsource_orig, sinksCI, pred_emnoise, num_sources, bootstrapping_iterations, include_epsilon = T, n_cores = 8,EM_ITERATIONS){
  
  em_itr=EM_ITERATIONS
  
  
  
  n_cores <- n_cores
  cl<- makeCluster(n_cores)
  registerDoSNOW(cl)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  em_bootstrapping <- c()
  predict_cls <- c();
  
  
  pb <- txtProgressBar(max = bootstrapping_iterations, style = 3)
  starttime<-Sys.time()
  em_bootstrapping<-foreach(k=1:bootstrapping_iterations, .combine='c', .export=ls(.GlobalEnv), .options.snow = opts,
                            .packages=c("Rcpp","RcppArmadillo")) %dopar% {
                              
                              Rcpp::cppFunction("arma::mat schur(arma::mat& a, arma::mat& b)
                                                {return(a % b); }", depends="RcppArmadillo")
                              
                              subsample = unique(sample(c(1:dim(totalsource_orig)[2]), replace = T))
                              
                              totalsource= totalsource_orig[,subsample]
                              source<- split(totalsource, seq(nrow(totalsource)))
                              #tst<-lapply(source, t)
                              #source<-lapply(source, as.matrix)
                              #samps <- source
                              samps<-lapply(source, t)
                              
                              observed_samps <- samps
                              observed_samps[[(num_sources + 1)]] = t(rep(0, dim(samps[[1]])[1]))
                              
                              initalphs<-runif(num_sources+1, 0.0, 1.0)
                              initalphs=initalphs/Reduce("+", initalphs)
                              sink_final=t(as.matrix(sinksCI[,subsample]))
                              
                              emnoise_start_time<-as.numeric(proc.time()[3])
                              tmp<-do_EM(alphas=initalphs, sources=samps, sink=sink_final, iterations=em_itr, observed=observed_samps)
                              pred_emnoise = tmp$toret
                              emnoise_rts = as.numeric(proc.time()[3])-emnoise_start_time
                              c(list(pred_emnoise), list(emnoise_rts))
                              
                              
                              
                            }
  endttime<-Sys.time()
  
  close(pb)
  stopCluster(cl)
  
  cls_rts<-split(em_bootstrapping, 1:2)[[2]]
  cls_rts=unlist(cls_rts)
  hold<-split(em_bootstrapping,1:2)[[1]] ##so c's as [[2]], mlists[[1]]
  hold<-unlist(hold)
  
  if(include_epsilon==F){
    cls_predictions <- matrix(hold, nrow = bootstrapping_iterations, ncol = num_sources, byrow = T)
  }else{
    cls_predictions <- matrix(hold, nrow = bootstrapping_iterations, ncol = num_sources+1, byrow = T)
  }
  
  CI_L = c()
  CI_H = c()
  for(j in 1:(num_sources+1)){
    
    #bias correction
    Z_0 = qnorm (length(which(cls_predictions[,j] < pred_emnoise[j]))/bootstrapping_iterations)
    
    if(Z_0 == Inf || Z_0 == -Inf){
      
      CI_L[j] = quantile(cls_predictions[,j], probs = 0.025)
      CI_H[j] = quantile(cls_predictions[,j], probs = 0.975)
    }
    
    else{
      #acceleration a
      a = sum( (mean(cls_predictions[,j]) - cls_predictions[,j])^3) / 6*(( sum( (mean(cls_predictions[,j]) - 
                                                                                   cls_predictions[,j])^2) )^(2/3))
      Z_alpha_h = qnorm(0.975)
      Z_alpha_l = qnorm(0.025)
      
      alpha_new_h = Z_0 + (Z_0 + Z_alpha_h)/(1-a*(Z_0 + Z_alpha_h))
      alpha_new_l = Z_0 + (Z_0 + Z_alpha_l)/(1-a*(Z_0 + Z_alpha_l))
      
      pnorm(alpha_new_h)
      pnorm(alpha_new_l)
      
      
      CI_L[j] = quantile(x = cls_predictions[,j], probs = pnorm(alpha_new_l))
      CI_H[j] = quantile(x = cls_predictions[,j], probs = pnorm(alpha_new_h))
      
    }
    
    CI = data.frame(CI_L, CI_H)
  }
  
  return(CI)
}

create_CI_EM_2 <- function(totalsource_orig, sinksCI, pred_emnoise, num_sources, bootstrapping_iterations, include_epsilon = T, n_cores = 8,EM_ITERATIONS){
  
  
  # totalsource_orig = em_source;
  # sinksCI = sink;
  # pred_emnoise = tmp$toret
  # bootstrapping_iterations = 10
  # num_sources=12
  em_itr=EM_ITERATIONS
  
  
  n_cores <- n_cores
  cl<- makeCluster(n_cores)
  registerDoSNOW(cl)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  em_bootstrapping <- c()
  predict_cls <- c();
  
  
  pb <- txtProgressBar(max = bootstrapping_iterations, style = 3)
  starttime<-Sys.time()
  em_bootstrapping<-foreach(k=1:bootstrapping_iterations, .combine='c', .export=ls(.GlobalEnv), .options.snow = opts,
                            .packages=c("Rcpp","RcppArmadillo")) %dopar% {
                              
                              Rcpp::cppFunction("arma::mat schur(arma::mat& a, arma::mat& b)
                                                {return(a % b); }", depends="RcppArmadillo")
                              
                              subsample = unique(sample(c(1:dim(totalsource_orig)[2]), replace = T))
                              
                              totalsource= totalsource_orig[,subsample]
                              source<- split(totalsource, seq(nrow(totalsource)))
                              #tst<-lapply(source, t)
                              #source<-lapply(source, as.matrix)
                              #samps <- source
                              samps<-lapply(source, t)
                              
                              observed_samps <- samps
                              observed_samps[[(num_sources + 1)]] = t(rep(0, dim(samps[[1]])[1]))
                              
                              initalphs<-runif(num_sources+1, 0.0, 1.0)
                              initalphs=initalphs/Reduce("+", initalphs)
                              sink_final=t(as.matrix(sinksCI[,subsample]))
                              
                              emnoise_start_time<-as.numeric(proc.time()[3])
                              tmp<-do_EM(alphas=initalphs, sources=samps, sink=sink_final, iterations=em_itr, observed=observed_samps)
                              pred_emnoise = tmp$toret
                              emnoise_rts = as.numeric(proc.time()[3])-emnoise_start_time
                              c(list(pred_emnoise), list(emnoise_rts))
                              
                              
                              
                            }
  endttime<-Sys.time()
  
  close(pb)
  stopCluster(cl)
  
  cls_rts<-split(em_bootstrapping, 1:2)[[2]]
  cls_rts=unlist(cls_rts)
  hold<-split(em_bootstrapping,1:2)[[1]] ##so c's as [[2]], mlists[[1]]
  hold<-unlist(hold)
  
  if(include_epsilon==F){
    cls_predictions <- matrix(hold, nrow = bootstrapping_iterations, ncol = num_sources, byrow = T)
  }else{
    cls_predictions <- matrix(hold, nrow = bootstrapping_iterations, ncol = num_sources+1, byrow = T)
  }
  
  CI_L = c()
  CI_H = c()
  for(j in 1:c(num_sources+1)){
    
    delta =  cls_predictions[,j] - pred_emnoise[j]
    
    # CI_H[j] = pred_emnoise[j] - quantile(x = delta, probs = 0.025)
    # CI_L[j] = pred_emnoise[j] - quantile(x = delta, probs = 0.975)
    # 
    CI_L[j] = quantile(cls_predictions[,j], probs = 0.025)
    CI_H[j] = quantile(cls_predictions[,j], probs = 0.975)
    
    
    
    CI = data.frame(CI_L, CI_H)
  }
  
  return(CI)
}

create_CI_CLS <- function(totalnew_source, sinks, pred_emnoise, bootstrapping_iterations){
  
  
  n_cores<-8
  cl<- makeCluster(n_cores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = bootstrapping_iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cls_bootstrapping <- c()
  predict_cls <- c(); 
  
  cls_bootstrapping<-foreach(k=1:bootstrapping_iterations, .combine='c', .export=ls(.GlobalEnv), .options.snow = opts,
                             .packages = "lsei") %dopar% {
                               
                               subsample = sample(c(1:dim(totalnew_source)[2]), replace = T)
                               X <- source_process_nounknown(totalnew_source[,subsample], envs_simulation, rarefaction_depth = COVERAGE)
                               X <- t(t(X)/rowSums(t(X)))
                               y<-t(as.matrix(sinks[1,subsample]))
                               y <- y/rowSums(y)
                               X_mat = as.matrix(X)
                               y_mat = as.matrix(y)
                               
                               test = lsei(a = X_mat, b = t(y_mat), e =  diag(length(envs_simulation)), 
                                           f = rep(0, length(envs_simulation)))
                               test[test < 0] = 0
                               predict_cls = test
                               predict_cls=predict_cls/sum(predict_cls)
                               c(list(predict_cls)) 
                             }
  
  close(pb)
  stopCluster(cl)
  
  
  hold<-split(cls_bootstrapping,1:1)[[1]] ##so c's as [[2]], mlists[[1]]
  hold<-unlist(hold)
  
  if(include_epsilon==F){
    cls_predictions <- matrix(hold, nrow = bootstrapping_iterations, ncol = num_sources, byrow = T)
  }else{
    cls_predictions <- matrix(hold, nrow = bootstrapping_iterations, ncol = num_sources+1, byrow = T)
  }
  
  CI_L = c()
  CI_H = c()
  for(j in 1:n_sources){
    
    #bias correction
    Z_0 = qnorm (length(which(cls_predictions[,j] < pred_emnoise[j]))/bootstrapping_iterations)
    
    #acceleration a
    a = sum( (mean(cls_predictions[,j]) - cls_predictions[,j])^3) / 6*(( sum( (mean(cls_predictions[,j]) - 
                                                                                 cls_predictions[,j])^2) )^(2/3))
    Z_alpha_h = qnorm(0.975)
    Z_alpha_l = qnorm(0.025)
    
    alpha_new_h = Z_0 + (Z_0 + Z_alpha_h)/(1-a*(Z_0 + Z_alpha_h))
    alpha_new_l = Z_0 + (Z_0 + Z_alpha_l)/(1-a*(Z_0 + Z_alpha_l))
    
    pnorm(alpha_new_h)
    pnorm(alpha_new_l)
    
    
    CI_L[j] = quantile(x = cls_predictions[,j], probs = pnorm(alpha_new_l))
    CI_H[j] = quantile(x = cls_predictions[,j], probs = pnorm(alpha_new_h))
    
    CI = data.frame(CI_L, CI_H)
  }
  return(CI)
}

plot_pie_chart <- function(source_prop, source_label, plot_name, CI){
  
  
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=5, face="bold")
    )
  
  for(j in 1:length(source_label[6:10])){
    
    if(j <= (num_sources+1)){
      
      source_label[j] = paste(source_label[j], 
                              round(source_prop[j], 4), "                           ")
    }
    
    source_label[j+num_sources+1] = paste(source_label[j+num_sources+1], " ", round(source_prop[j+num_sources+1], 4)," [", 
                                          round(CI$CI_L[j], 4), ", ",
                                          round(CI$CI_H[j], 4), "]", sep = "")
  }
  
  
  dat = data.frame(count=source_prop[1:5]*100, category=source_label[1:5])
  dat$fraction = round(dat$count / sum(dat$count), 3)
  
  dat$ymax = cumsum(dat$fraction)
  dat$ymin = c(0, head(dat$ymax, n=-1))
  
  dat$category <- factor(dat$category, levels = unique(dat$category))
  
  
  p1 = ggplot(dat, aes(fill=category, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
    geom_rect(color='blue') +
    coord_polar(theta="y") +
    xlim(c(1, 4)) 
  
  
  edu1<-p1 + scale_fill_brewer("Source", palette="Set1") + blank_theme +
    theme(axis.text.x=element_blank()) + 
    ggtitle("") +
    theme(panel.grid=element_blank()) +
    theme(axis.text=element_blank()) +
    theme(axis.ticks=element_blank()) +
    # theme(legend.title = element_text(size=16, face="bold")) +
    # theme(legend.text = element_text(size = 14, face = "bold")) +  
    # theme(legend.position=c(1, 0)) +
    geom_label(aes(label=paste(fraction*100,"%"),x=3.5,y=(ymin+ymax)/2),inherit.aes = TRUE, 
               show.legend = FALSE)
  
  
  #####Dat2 
  dat = data.frame(count= source_prop[6:10]*100, category=source_label[6:10])
  dat$fraction = round(dat$count / sum(dat$count), 3)
  
  dat$ymax = cumsum(dat$fraction)
  dat$ymin = c(0, head(dat$ymax, n=-1))
  
  dat$category <- factor(dat$category, levels = unique(dat$category))
  
  
  p2 = ggplot(dat, aes(fill=category, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
    geom_rect(color='blue') +
    coord_polar(theta="y") +
    xlim(c(1, 4)) 
  
  
  edu2<-p2 + scale_fill_brewer("Source", palette="Set1") + blank_theme +
    theme(axis.text.x=element_blank()) + ggtitle("") +
    theme(panel.grid=element_blank()) +
    theme(axis.text=element_blank()) +
    theme(axis.ticks=element_blank()) +
    # theme(legend.title = element_text(size=16, face="bold")) +
    # theme(legend.text = element_text(size = 14, face = "bold")) +
    # theme(legend.position=c(1, 0)) +
    geom_label(aes(label=paste(fraction*100,"%"),x=3.5,y=(ymin+ymax)/2),inherit.aes = TRUE, 
               show.legend = FALSE)
  plot_grid(edu1, edu2, labels = c(paste('ST', plot_name), paste('EM', plot_name)),  nrow = 2, align = 'v')
  
  ggsave(filename = paste0(plot_name, ".jpeg"), 
         plot = plot_grid(edu1, edu2, labels = c(paste('ST', plot_name), paste('EM', plot_name)),  nrow = 2, align = 'v'))
  
}

plot_pie_chart_new <- function(source_prop, source_label, plot_name, CI){
  
  
  
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=3, face="bold")
    )
  
  
  for(j in 1:length(source_label[6:10])){
    
    if(j <= (num_sources+1)){
      
      source_label[j] = paste(source_label[j],
                              round(source_prop[j], 3), "                        ")
    }
    
    source_label[j+num_sources+1] = paste(source_label[j+num_sources+1], " ", 
                                          round(source_prop[j+num_sources+1], 3)," [",
                                          round(CI$CI_L[j], 3), ", ",
                                          round(CI$CI_H[j], 3), "]", sep = "")
  }
  
  
  
  dat = data.frame(count=source_prop[1:5]*100, category=source_label[1:5])
  dat$fraction = round(dat$count / sum(dat$count), 3)
  dat$fraction_known = c(dat$fraction[1:(length(dat$fraction) - 1)]/sum(dat$fraction[1:(length(dat$fraction) - 1)]), 0)
  
  dat$ymax = cumsum(dat$fraction)
  dat$ymin = c(0, head(dat$ymax, n=-1))
  
  dat$ymax_known = cumsum(dat$fraction_known)
  dat$ymin_known = c(0, head(dat$ymax_known, n=-1))
  
  dat$category <- factor(dat$category, levels = unique(dat$category))
  
  
  p1 = ggplot(dat, aes(fill=category, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
    geom_rect(color='black') +
    coord_polar(theta="y") +
    xlim(c(1, 4)) 
  
  
  edu1<-p1 + scale_fill_brewer("Source", palette="Set1") + blank_theme +
    theme(axis.text.x=element_blank()) + 
    ggtitle("") +
    theme(panel.grid=element_blank()) +
    theme(axis.text=element_blank()) +
    theme(axis.ticks=element_blank()) +
    geom_label_repel(aes(label=paste(fraction*100,"%"),
                         x=3.5,y=(ymin+ymax)/2),inherit.aes = TRUE, 
                     show.legend = FALSE) 
  
  
  dat_new = dat[c(1:4),]
  dat_new$fraction_known = round(dat_new$fraction_known, 2)
  
  
  p1.1 = ggplot(dat_new, aes(fill=category, ymax=ymax_known, ymin=ymin_known, xmax=4, xmin=3)) +
    geom_rect(color='black') +
    coord_polar(theta="y") +
    xlim(c(1, 4)) 
  
  
  edu1.1<-p1.1 + scale_fill_brewer("Source", palette="Set1") + blank_theme +
    theme(axis.text.x=element_blank()) + 
    ggtitle("") +
    theme(panel.grid=element_blank()) +
    theme(axis.text=element_blank()) +
    theme(axis.ticks=element_blank()) +
    # theme(legend.title = element_text(size=16, face="bold")) +
    # theme(legend.text = element_text(size = 14, face = "bold")) +  
    # theme(legend.position=c(1, 0)) +
    # geom_label(aes(label=paste(fraction_known*100,"%"),x=3.5,y=(ymin_known+ymax_known)/2),inherit.aes = TRUE, 
    #            show.legend = FALSE)  + theme(legend.position="none") 
    geom_label_repel(aes(label=paste(fraction_known*100,"%"),
                         x=3.5,y=(ymin_known+ymax_known)/2),inherit.aes = TRUE, 
                     show.legend = FALSE) + theme(legend.position="none") 
  
  
  #####Dat2 
  
  dat = data.frame(count=source_prop[6:10]*100, category=source_label[6:10])
  dat$fraction = round(dat$count / sum(dat$count), 3)
  dat$fraction_known = c(dat$fraction[1:(length(dat$fraction) - 1)]/sum(dat$fraction[1:(length(dat$fraction) - 1)]), 0)
  
  dat$ymax = cumsum(dat$fraction)
  dat$ymin = c(0, head(dat$ymax, n=-1))
  
  dat$ymax_known = cumsum(dat$fraction_known)
  dat$ymin_known = c(0, head(dat$ymax_known, n=-1))
  
  dat$category <- factor(dat$category, levels = unique(dat$category))
  
  dat_new = dat[c(1:4),]
  dat_new$fraction_known = round(dat_new$fraction_known, 2)
  
  p2 = ggplot(dat, aes(fill=category, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
    geom_rect(color='black') +
    coord_polar(theta="y") +
    xlim(c(1, 4)) 
  
  
  edu2<-p2 + scale_fill_brewer("Source", palette="Set1") + blank_theme +
    theme(axis.text.x=element_blank()) + ggtitle("") +
    theme(panel.grid=element_blank()) +
    theme(axis.text=element_blank()) +
    theme(axis.ticks=element_blank()) +
    # theme(legend.title = element_text(size=16, face="bold")) +
    # theme(legend.text = element_text(size = 14, face = "bold")) +
    # theme(legend.position=c(1, 0)) +
    # geom_label(aes(label=paste(fraction*100,"%"),x=3.5,y=(ymin+ymax)/2), label.padding = unit(0.2, "lines"),
    #            inherit.aes = TRUE,
    #            position = "identity",
    #            show.legend = FALSE)
    geom_label_repel(aes(label=paste(fraction*100,"%"),
                         x=3.5,y=(ymin+ymax)/2),inherit.aes = TRUE, 
                     show.legend = FALSE) 
  
  
  p2.1 = ggplot(dat_new, aes(fill=category, ymax=ymax_known, ymin=ymin_known, xmax=4, xmin=3)) +
    geom_rect(color='black') +
    coord_polar(theta="y") +
    xlim(c(1, 4)) 
  
  
  edu2.1<-p2.1 + scale_fill_brewer("Source", palette="Set1") + blank_theme +
    theme(axis.text.x=element_blank()) + ggtitle("") +
    theme(panel.grid=element_blank()) +
    theme(axis.text=element_blank()) +
    theme(axis.ticks=element_blank()) +
    # theme(legend.title = element_text(size=16, face="bold")) +
    # theme(legend.text = element_text(size = 14, face = "bold")) +
    # theme(legend.position=c(1, 0)) +
    # geom_label(aes(label=paste(fraction_known*100,"%"),
    #                x=3.5,y=(ymin_known+ymax_known)/2),inherit.aes = TRUE, 
    #            show.legend = FALSE) + theme(legend.position="none")
    geom_label_repel(aes(label=paste(fraction_known*100,"%"),
                         x=3.5,y=(ymin_known+ymax_known)/2),inherit.aes = TRUE, 
                     show.legend = FALSE) + theme(legend.position="none") 
  
  
  plot_grid(edu1.1, edu1,edu2.1, edu2,  nrow = 2, align = 'h')
  
  
  ggsave(filename = paste0(plot_name, ".jpeg"), 
         plot =   plot_grid(edu1.1, edu1, edu2.1, edu2,  nrow = 2, align = 'h'))
  
}



