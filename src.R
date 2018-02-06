library(ggplot2)
library(cowplot)


"change_C"<-function(newcov, X){
  
 # newcov = COVERAGE; X = totalnew_source[1,]
  
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
  
  # x = totalsource
  # maxdepth = COVERAGE
  
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

"kld" <- function(p,q){
  nonzero <- p>0 & q>0
  return(sum(p[nonzero] * log2(p[nonzero]/q[nonzero])))    
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
  
  # alphas = curalphas; sources = sources; sink = sink; observed = observed
  
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
  #sink<-as.matrix(sink); #src1<-as.matrix(sources[[1]]); src2<-as.matrix(sources[[2]])
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
  
  # alphas=initalphs; sources=samps; sink=sink; iterations=em_itr; observed = samps
  
  curalphas<-alphas
  newalphas<-alphas
  m_guesses<-c(alphas[1])
  #print(paste("first guess: ", m_guesses[1]))
  #prevll<-LL(alphas=newalphas,sources=sources, sink=sink)
  for(itr in 1:iterations){
    
    curalphas<-E(newalphas, sources)
    tmp <- M(alphas = curalphas, sources = sources, sink = sink, observed = observed)
    newalphas <- tmp$new_alpha
    sources <- tmp$new_sources  
    #print(paste("length newalphas: ", length(newalphas)))
    #print(paste("new alphas: ", newalphas))
    m_guesses<-c(m_guesses, newalphas[1])
    #print(paste("length mguess: ", length(m_guesses)))
    #print(paste("itr: ", itr))
    #print(paste("New alphas: ", newalphas[1]))
    #newll<-LL(alphas=newalphas,sources=sources, sink=sink)
    #print(paste("newlikelihood: ", LL(alphas = newalphas, sources = sources, sink = sink)))
    #print(paste("guessprev: ", m_guesses[length(m_guesses)], "guesscur: ", m_guesses[length(m_guesses)-1]))
    if(abs(m_guesses[length(m_guesses)]-m_guesses[length(m_guesses)-1])<=10^-6) break
    #prevll<-newll
  }         #                             this will return epsilon
  toret<-c(newalphas)
  results <- list(toret = toret, sources = sources)
  
  return(results)
}

"M_basic"<-function(alphas, sources, sink){
  newalphs<-c()
  #sink<-as.matrix(sink); #src1<-as.matrix(sources[[1]]); src2<-as.matrix(sources[[2]])
  #sources<-lapply(sources, as.matrix)
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
  #print(paste("first guess: ", m_guesses[1]))
  #prevll<-LL(alphas=newalphas,sources=sources, sink=sink)
  for(itr in 1:iterations){
    curalphas<-E(newalphas, sources)
    newalphas<-M_basic(curalphas, sources, sink)
    #print(paste("length newalphas: ", length(newalphas)))
    #print(paste("new alphas: ", newalphas))
    m_guesses<-c(m_guesses, newalphas[1])
    #print(paste("length mguess: ", length(m_guesses)))
    #print(paste("itr: ", itr))
    #print(paste("New alphas: ", newalphas[1]))
    #newll<-LL(alphas=newalphas,sources=sources, sink=sink)
    #print(paste("newlikelihood: ", LL(alphas = newalphas, sources = sources, sink = sink)))
    #print(paste("guessprev: ", m_guesses[length(m_guesses)], "guesscur: ", m_guesses[length(m_guesses)-1]))
    if(abs(m_guesses[length(m_guesses)]-m_guesses[length(m_guesses)-1])<=10^-6) break
    #prevll<-newll
  }         #                             this will return epsilon
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

"read_pseudo_data"<-function(dataset, path_to_data){
  
  # if(server_flag == 0)
  #   path_to_data<-"~/Dropbox/Source Tracking/Published_Data/"
  # else
  #  path_to_data <-"/home/orif/CowsTimeSeriesClusterer/run1_to_21_subsampled_10k_with_taxonomy/Paper_analysis/SourceTracker/Source_Tracker_sim/Published_Data/"
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
  
  # num_sources = 50
  
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
    other_ms = resid/(num_sources-3)
    m = rep(NA, num_sources)
    m[index] = c(m_1, m_2, m_3)
    m[is.na(m)] = other_ms
    
  }
  subsum = 0
  idx = 1:length(m)
  while (subsum < EPSILON){
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

unknown_initialize_1 <- function(sources, sink, n_sources){

  # sources = totalsource[c(1:n_sources),]
  # sink = as.numeric(sink)

  unknown_source = rep(0, length(sink))

  #zero all the OTUs with at least 1 known source
  sources_sum = apply(sources, 2 ,sum)
  ind_known_source_abun = which(sources_sum > 0)
  unknown_source[ind_known_source_abun] = 0


  #Select the cor OTUs
  ind_cor = list()
  ind_known_source_abun = c()
  ind_cor_all = which(sources[1,] > 0)
  for(j in 1:n_sources){

    ind_cor[[j]] = which(sources[j,] > 0)


    if(j > 1){

      tmp = ind_cor_all
      ind_cor_all = intersect(ind_cor_all ,ind_cor[[j]])

    }
  }

  cor_abundance = apply(sources[,ind_cor_all], 2, median) #take the median abundnace of the 'cor'
  unknown_source[ind_cor_all] = cor_abundance


  #keep the sink abundance where there is no known source
  ind_no_known_source_abun = which(sources_sum == 0)

  for(j in 1:length(ind_no_known_source_abun)){

    # unknown_source[ind_no_known_source_abun[j]] = max(runif(n = 1, min = 1, max = 100), sink[ind_no_known_source_abun[j]])
    unknown_source[ind_no_known_source_abun[j]] = sink[ind_no_known_source_abun[j]]

  }



  return(unknown_source)

}


unknown_initialize <- function(sources, sink, n_sources){
  
   sources = totalnew_source_old[c(1:n_sources),];
   sink = as.numeric(sinks[i,]);
   n_sources = n_sources 
  
  unknown_source = rep(0, length(sink))
  sum_sources = apply(sources, 2, sum)
  
  unknown_source = c()
  
  for(j in 1:length(sum_sources)){
    
    unknown_source[j] = max(sink[j]-sum_sources[j], 0)
    
  }
  
  
  # is.cor<- function(vec, n = (n_sources-1)){
  # 
  #   if(length(which(vec > 0)) >= n)
  #     return(1)
  #   else
  #     return(0)
  # }
  # 
  # cor_ind = apply(sources, 2, is.cor)
  # 
  # cor_index = which(cor_ind > 0)
  # 
  # cor_abundance = round(apply(sources[,cor_index], 2, min)) #take the median abundnace of the 'cor'
  # unknown_source[cor_index] = cor_abundance
  
  #Select the cor OTUs
  # ind_cor = list()
  # ind_known_source_abun = c()
  # ind_cor_all = which(sources[1,] > 0)
  # for(j in 1:n_sources){
  # 
  #   ind_cor[[j]] = which(sources[j,] > 0)
  # 
  # 
  #   if(j > 1){
  # 
  #     tmp = ind_cor_all
  #     ind_cor_all = intersect(ind_cor_all ,ind_cor[[j]])
  # 
  #   }
  # }
  # 
  # cor_abundance = apply(sources[,ind_cor_all], 2, median) #take the median abundnace of the 'cor'
  # unknown_source[ind_cor_all] = cor_abundance
  
  # unknown_source = t(unknown_source)
  
  # Data = data.frame(t(sources), unknown_source ,sink)
  # str(Data)
  # View(Data)
  
  return(unknown_source)
  
}

create_CI_EM <- function(totalsource_orig, sinks, pred_emnoise, bootstrapping_iterations){
  
  n_cores<-8
  cl<- makeCluster(n_cores)
  registerDoSNOW(cl)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  em_bootstrapping <- c()
  predict_cls <- c(); 
  
  pb <- txtProgressBar(max = bootstrapping_iterations, style = 3)
  em_bootstrapping<-foreach(k=1:bootstrapping_iterations, .combine='c', .export=ls(.GlobalEnv), .options.snow = opts,
                            .packages=c("Rcpp","RcppArmadillo")) %dopar% {
                              
                              Rcpp::cppFunction("arma::mat schur(arma::mat& a, arma::mat& b)
                                                {return(a % b); }", depends="RcppArmadillo")
                              
                              subsample = sample(c(1:dim(totalsource)[2]), replace = T)
                              
                              
                              totalsource= totalsource_orig[,subsample]
                              source<- split(totalsource, seq(nrow(totalsource)))
                              source<-lapply(source, as.matrix)
                              samps <- source
                              samps<-lapply(samps, t)
                              
                              observed_samps <- samps
                              observed_samps[[(num_sources + 1)]] = t(rep(0, dim(samps[[1]])[1]))
                              
                              initalphs<-runif(num_sources+1, 0.0, 1.0)
                              initalphs=initalphs/Reduce("+", initalphs)
                              sink=t(as.matrix(sinks[,subsample]))
                              
                              emnoise_start_time<-as.numeric(proc.time()[3])
                              tmp<-do_EM(alphas=initalphs, sources=samps, sink=sink, iterations=em_itr, observed=observed_samps)
                              pred_emnoise = tmp$toret
                              emnoise_rts = as.numeric(proc.time()[3])-emnoise_start_time
                              c(list(pred_emnoise), list(emnoise_rts)) 
                              
                              
                              
                            }
  
  
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

create_CI_CLS <- function(totalnew_source, sinks, pred_emnoise, bootstrapping_iterations){
  
  # totalnew_source = totalsource;
  # sinks = sinks; 
  # pred_emnoise = pred_emnoise;
  # bootstrapping_iterations = 100
  
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
  

  # source_prop = c(results$proportions, as.numeric(pred_emnoise))
  # source_label = c(rep(c(as.character(unique(envs)), "Unknown"), 2))
  # plot_name = paste(metadata[k[it] ,1], metadata[k[it] ,2] ,length(which(sinks > 0)), "OTUs"
  #                   , "depth", sum(sinks), sep = "_")
  
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

