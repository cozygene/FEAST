
library(glmnet)

#E - step
E <- function(alphas, sources, sink=NA, observed=NA){
  nums<-(sapply(1:length(alphas), function(n) Reduce("+", crossprod(as.numeric(alphas[n]),as.numeric(sources[[n]])))))
  denom<-(Reduce("+", nums))
  return(nums/denom)
}

#M-step
M <- function(alphas, sources, sink, observed){

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
    num_list[[i]]<- num[[1]][1,] + observed[[i]][1,]

    denom <- Reduce("+",unlist(num_list[[i]]))
    source_new[[i]] <- num_list[[i]]/denom
    source_new[[i]][is.na(source_new[[i]])] <- 0
  }

  sources <- source_new

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


M_basic <- function(alphas, sources, sink){

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

alpha.solve.l1 <- function(sink, pij, lambda=1) {
  # sparse constraint on 1:M sources and sum to 1
  # print(paste('Lambda:', lambda))

  # saveRDS(pij, 'pij.rds')
  # saveRDS(sink, 'sink.rds')

  posinds <- sink > 0 & colSums(pij) > 0
  nsources <- dim(pij)[1]
  ntaxa <- sum(posinds)

  cvx.avector <- Variable(nsources)

  terms <- list()
  for (ii in 1:(nsources)) {
    wi <- c(sink)[posinds] * pij[ii,posinds]
    term <- wi * log(cvx.avector[ii])
    terms[[ii]] <- sum(term)
  }
  term <- Reduce('+', terms)

  constraints <- list(
    sum(cvx.avector) <= 1,
    cvx.avector >= 0
  )

  penalty <- lambda * p_norm(cvx.avector[1:(nsources-1)], 1)
  prob <- Problem(Maximize(term - penalty), constraints)
  use_solver <- 'ECOS'
  if (dim(pij)[1] > 100) {
    # ECOS default, SCS scales better for larger inputs
    use_solver <- 'SCS'
  }
  result <- solve(prob, solver = use_solver)
  out <- result$getValue(cvx.avector)

  # print(paste('Max:', result$value))

  return (out)
}

E_stensl <- function(alphas, sources, sink=NA, observed=NA){

  ll <- NA
  if (exists('pij')) {
    # full EM procedure as defined in methods (can be slow)
    nsources <- length(alphas)
    ntaxa <- length(sources[[1]])
    pij <<- matrix(0, nrow=nsources, ncol=ntaxa)
    zero <- 0.000001
    pijNonZero <- matrix(0, nrow=nsources, ncol=ntaxa)


    old.gamma <- do.call(rbind, sources)
    denom <- crossprod(alphas, old.gamma)
    for (ii in 1:nsources) {
      val <- alphas[ii] * old.gamma[ii,] / denom
      pij[ii,] <<- val
      val[is.na(val) | val < zero] <- zero
      pijNonZero[ii,] <- val
    }

    if (!all(is.na(sink)) & !all(is.na(observed))) {
      # we can compute the E[logL] for verification
      term1 <- 0
      for (ii in 1:nsources) {
        logterm <- alphas[ii] * old.gamma[ii,]
        logterm[logterm < zero] <- zero
        val <- sum(sink * pijNonZero[ii,] * log(logterm))
        term1 <- term1 + val
      }
      term2 <- 0
      for (ii in 1:(nsources-1)) {
        # unknown counts are not part of the LL
        logterm <- old.gamma[ii,]
        logterm[logterm < zero] <- zero
        val <- sum(observed[[ii]]*log(logterm))
        term2 <- term2 + val
      }
      ll <- term1 + term2
      ll <- ll - sparse.lambda * sum(alphas[1:(nsources-1)]) # penalty term (from 1:M)
      if (exists('ll.list')) ll.list <<- c(ll.list, ll)
      # print(paste('E[logL]:', ll))
    }
  }

  nums<-(sapply(1:length(alphas), function(n) Reduce("+", crossprod(as.numeric(alphas[n]),as.numeric(sources[[n]])))))
  denom<-(Reduce("+", nums))
  return(nums/denom)
}

M_stensl <- function(alphas, sources, sink, observed){

  alphas.sparse <- NA
  # if (useSparseMax & exists('pij')) {
  if (T) {
    pij[is.na(pij)] <- 0

    # FEAST alphas using max step:
    alphas.closedform <- rep(0, length(alphas))
    for (ii in 1:length(alphas)) {
      pijrow <- pij[ii,]
      pijrow[is.na(pijrow)] <- 0
      val <- sum(sink * pijrow)
      alphas.closedform[ii] <- val
    }
    alphas.closedform <- alphas.closedform / sum(alphas.closedform)

    # Sparse alphas using solver:
    alphas.sparse <- tryCatch({
      alpha.solve.l1(sink, pij, lambda=sparse.lambda)
    },
    error=function(e) {
      print('Convex solver exited.')
      print(e)
      alphas
    })
  }

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
    num_list[[i]]<- num[[1]][1,] + observed[[i]][1,]

    denom <- Reduce("+",unlist(num_list[[i]]))
    source_new[[i]] <- num_list[[i]]/denom
    source_new[[i]][is.na(source_new[[i]])] <- 0
  }

  sources <- source_new

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

  # if (useSparseMax) {
  Results$new_alpha <- alphas.sparse / sum(alphas.sparse)
  feast.alpha <- Results$new_alpha

  return(Results)
}

do_EM <-function(alphas, sources, observed, sink, iterations, E_fn=E, M_fn=M){

  curalphas<-alphas
  newalphas<-alphas
  # m_guesses<-c(alphas[1])
  for(itr in 1:iterations){
    prevalphas = newalphas
    curalphas<-E_fn(newalphas, sources, sink, observed)
    tmp <- M_fn(alphas = curalphas, sources = sources, sink = sink, observed = observed)
    newalphas <- tmp$new_alpha
    sources <- tmp$new_sources

    # m_guesses<-c(m_guesses, newalphas[1])
    # if(abs(m_guesses[length(m_guesses)]-m_guesses[length(m_guesses)-1])<=10^-6)
    #   break
    # More comprehensive stopping conds
    if (itr > 1) {
      if (max(abs(prevalphas-newalphas)) < 10^-6) {
        print('Convergeance criteria.')
        break
      }

      lval = ifelse(exists('ll.list'), ll.list[length(ll.list)], 0)
      lambda.info = ifelse(exists('sparse.lambda'), sprintf('Lambda:%.4f', sparse.lambda), '')
      print(sprintf('[Iter %d/%d] ElogL:%.2f Unk:%.3f %s',
        itr-1, iterations,
        lval,
        newalphas[length(newalphas)],
        lambda.info))
    }
  }
  ll.list <<- ll.list[2:length(ll.list)] # disregard first ll
  toret<-c(newalphas)
  results <- list(toret = toret, sources = sources)

  return(results)
}



do_EM_basic <- function(alphas, sources, sink, iterations){
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

# unknown source initialization
unknown_initialize_1 <- function(sources, sink, n_sources){


  unknown_source <- rep(0, length(sink))
  sources_sum <- apply(sources, 2 ,sum)


  unknown_source <- c()

  if(n_sources > 1){

    for(j in 1:length(sources_sum)){

      unknown_source[j] <- max(sink[j]-sources_sum[j], 0)

    }

    #Select the cor taxa
    ind_cor <- list()
    ind_known_source_abun <- c()
    ind_cor_all <- which(sources[1,] > 0)

    counter <- matrix(0, ncol = dim(sources)[2], nrow =  dim(sources)[1])


    for(j in 1:n_sources){

      ind_cor[[j]] <- which(sources[j,] > 0)

      for(k in 1:length(sources[j,])){

        if(sources[j,k] > 0){

          counter[j,k] <- counter[j,k]+1
        }


      }

    }

    OTU_present_absent <- apply(counter, 2, sum)
    ind_cor_all <- which(OTU_present_absent >= round(n_sources*0.5))

    if(length(ind_cor_all) > 1){

      cor_abundance <- round(apply(sources[,ind_cor_all], 2, min)/2) #take the min abundnace of the 'cor'
      unknown_source[ind_cor_all] <- cor_abundance


    }

  }


  #keep the sink abundance where there is no known source
  ind_no_known_source_abun <- which(sources_sum == 0)

  for(j in 1:length(ind_no_known_source_abun)){

    unknown_source[ind_no_known_source_abun[j]] <- max((sink[ind_no_known_source_abun[j]] - rpois(n = 1, lambda = 0.5)), 0)

  }

  unknown_source[is.na(unknown_source)] <- 0





  return(unknown_source)

}

# unknown source initialization

unknown_initialize <- function(sources, sink, n_sources){

  unknown_source <- rep(0, length(sink))
  sum_sources <- apply(sources, 2, sum)

  unknown_source <- c()

  for(j in 1:length(sum_sources)){

    unknown_source[j] <- max(sink[j]-sum_sources[j], 0)

  }

  if(n_sources == 1)
    unknown_source <- rpois(n = length(sinks), lambda = 0.05)

  return(unknown_source)

}

lsq_glmnet_l1l2 <- function(
  sink, sources,
  l1l2=1, lambda=1e-6,
  normalize=T) {

  if (normalize) {
    sources <- sources / rowSums(sources)
    sink <- sink / sum(sink)
  }

  Amat <- model.matrix(~Sources, list(Sources=t(sources)))
  result <- glmnet(Amat, sink, alpha=l1l2, lambda=lambda, lower.limits=0, intercept=F)

  contr <- result$beta[2:length(result$beta)]

  return (contr)
}

lsq_procedure <- function(sources, sink, default_known=0.99) {
  weights_l1 <- lsq_glmnet_l1l2(
    sink,
    sources,
    normalize=F,
    l1l2=1)

  # 1. The alpha init
  lsq_init <- weights_l1 / sum(weights_l1)
  alpha_init <- c(lsq_init*default_known, (1-default_known)) # appending default unknown p

  # 2. The unknown init
  max_source <- sources[which.max(lsq_init),]

  # We scale the max source according to its given weight
  mixed_max <- max_source * max(lsq_init)
  unknown_init <- sink - mixed_max

  unknown_init[unknown_init < 0] <- 0   # clip at zero

  return (list(
    lsq=weights_l1,
    alpha=alpha_init,
    unknown=unknown_init
  ))
}

LsqResid.Procedure <- function(rare_sink, rare_sources) {
  print('Performing Lsq+Residual Procedure...')

  blob <- lsq_procedure(
    rare_sources,
    rare_sink
  )

  # save the max-estimated unknown for later analysis
  blob$unknown.maxest <- blob$unknown

  # compute Lsq error
  nonzero <- rare_sink != 0
  slen <- length(rare_sink)
  reconst.sink <- matrix(blob$lsq %*% rare_sources, ncol=slen)
  l1_error <- t(matrix(rare_sink, ncol=slen) - reconst.sink)

  # Get the top N% most error taxa
  keep.amount = 0.75
  noneg.error <- sapply(l1_error, function(v) ifelse(v < 0, 0, v))
  error.ord <- rev(order(noneg.error)) # largest to smallest
  cutoff.index <- 1
  if (sum(noneg.error) > 0) {
    for (cutoff.index in 1:length(error.ord)) {
      if (sum(l1_error[error.ord[1:cutoff.index]]) / sum(noneg.error) >= keep.amount) {
        break
      }
    }
  }
  topn_resid_inds <- ifelse(1:length(l1_error) %in% error.ord[1:cutoff.index], T, F)

  keep_inds <- topn_resid_inds

  # save some intermediary values for later inspection
  blob$reconst.sink <- reconst.sink
  blob$l1.error <- l1_error

  return (blob)
}

Infer.SourceContribution <- function(source = sources_data, sinks = sinks, em_itr = 1000, env = rownames(sources_data), include_epsilon = T,
                  COVERAGE, method='feast', unknown_initialize_flag = 1){

  tmp <- source
  test_zeros <- apply(tmp, 1, sum)
  ind_to_use <- as.numeric(which(test_zeros > 0))
  ind_zero <- as.numeric(which(test_zeros == 0))

  source <- tmp[ind_to_use,]
  sinks <- sinks

  alpha_init = NA
  unknown_init = NA
  if (method == 'stensl') {
    print('Finding STENSL init...')

    ll.list <<- c() # reset ll history

    blob <- LsqResid.Procedure(sinks, source)

    # rarefy unknown to the same level as the experiment
    if (sum(blob$unknown) > 0) {
      blob$unknown <- blob$unknown * (COVERAGE/sum(blob$unknown))
    }

    alpha_init <- blob$alpha
    unknown_init <- blob$unknown
  }


  #####adding support for multiple sources#####
  if(length(dim(source)[1]) > 0)
    totalsource<-source
  if(length(dim(source)[1]) == 0)
    totalsource<-t(as.matrix(source))

  totalsource<-as.matrix(totalsource)
  sources <- split(totalsource, seq(nrow(totalsource)))
  sources<-lapply(sources, as.matrix)
  dists<-lapply(sources, function(x) x/(sum(colSums(x))))
  totaldist<-t(Reduce("cbind", dists))
  sinks<-matrix(sinks, nrow = 1, ncol = dim(totalsource)[2])

  if(length(dim(source)[1]) == 0)
    num_sources <- 1
  if(length(dim(source)[1]) > 0)
    num_sources <- dim(source)[1]
  envs_simulation <- c(1:(num_sources))

  source_old <- source
  totalsource_old <- totalsource

  source_old<-lapply(source_old,t)
  source_old<- split(totalsource_old, seq(nrow(totalsource_old)))
  source_old<-lapply(source_old, as.matrix)

  #Creating the unknown source per mixing iteration
  if(include_epsilon == TRUE){

    ##Adding the initial value of the unknown source for CLS and EM
    source_2 <- list()
    totalsource_2 <- matrix(NA, ncol = dim(totalsource_old)[2], nrow = ( dim(totalsource_old)[1] + 1))

    for(j in 1:num_sources){

      source_2[[j]] <- source_old[[j]]
      totalsource_2[j,] <- totalsource_old[j,]
    }

    #create unknown for each sink i

    sinks_rarefy <- FEAST_rarefy(matrix(sinks, nrow = 1), maxdepth = apply(totalsource_old, 1, sum)[1]) #make

    if(num_sources > 1){

      if(unknown_initialize_flag == 1)
        unknown_source <- unknown_initialize_1(sources = totalsource[c(1:num_sources),], sink = as.numeric(sinks),
                                              n_sources = num_sources)


      if(unknown_initialize_flag == 0)
        unknown_source <- unknown_initialize(sources = totalsource[c(1:num_sources),], sink = as.numeric(sinks),
                                            n_sources = num_sources)
    }

    if(num_sources == 1){

      if(unknown_initialize_flag == 1)
        unknown_source <- unknown_initialize_1(sources = t(as.matrix(totalsource[c(1:num_sources),])), sink = as.numeric(sinks),
                                              n_sources = num_sources)


      if(unknown_initialize_flag == 0)
        unknown_source <- unknown_initialize(sources = t(as.matrix(totalsource[c(1:num_sources),])), sink = as.numeric(sinks),
                                            n_sources = num_sources)
    }

    if(unknown_initialize_flag == 2 & !is.na(sum(unknown_init))) {
      print('Using STENSL unknown')
      unknown_source <- unknown_init
    }

    unknown_source_rarefy <- FEAST_rarefy(matrix(unknown_source, nrow = 1), maxdepth = COVERAGE)
    source_2[[j+1]] <- t(unknown_source_rarefy)
    totalsource_2[(j+1),] <- t(unknown_source_rarefy)
    totalsource <- totalsource_2

    source=lapply(source_2,t)
    source<- split(totalsource, seq(nrow(totalsource_2)))
    source<-lapply(source_2, as.matrix)

    envs_simulation <- c(1:(num_sources+1))

  }


  samps <- source
  samps<-lapply(samps, t)

  observed_samps <- samps
  observed_samps[[(num_sources + 1)]] <- t(rep(0, dim(samps[[1]])[2]))


  # initalphs<-runif(num_sources+1, 0.0, 1.0)
  # for reproducibility, recommending fixed uniform init
  initalphs<-rep(1/(num_sources+1), num_sources+1)
  initalphs=initalphs/Reduce("+", initalphs)
  sink_em <- as.matrix(sinks)

  if (!is.na(sum(alpha_init))) {
    print('Using STENSL alpha')
    alpha_init <- blob$alpha
  }

  # pred_em<-do_EM_basic(alphas=initalphs, sources=samps, sink=sink_em, iterations=em_itr)

  tmp<-do_EM(
    alphas=initalphs,
    sources=samps,
    sink=sink_em,
    iterations=em_itr,
    observed=observed_samps,
    E_fn=ifelse(method == 'stensl', E_stensl, E),
    M_fn=ifelse(method == 'stensl', M_stensl, M)
  )
  pred_emnoise <- tmp$toret

  k <- 1
  pred_emnoise_all <- c()
  # pred_em_all <- c()

  for(j in 1:length(env)){

    if(j %in% ind_to_use){

      pred_emnoise_all[j] <- pred_emnoise[k]
      # pred_em_all[j] <- pred_em[k]
      k <- k+1

    }

    else{

      pred_emnoise_all[j] <- 0
      # pred_em_all[j] <- 0
    }

  }

  pred_emnoise_all[j+1] <- pred_emnoise[k]
  # pred_em_all[j+1] <- pred_em[k]



  names(pred_emnoise_all) <- c(env,"unknown")
  # names(pred_em_all) <- c(env,"unknown")


  Results <- list(unknown_source = unknown_source, unknown_source_rarefy = unknown_source_rarefy,
                 data_prop = data.frame(pred_emnoise_all))
  return(Results)

}


SharedTaxa_function <- function(count_data = otus,
                                Taxonomy = taxonomy,
                                all_sources = train.ix,
                                SpecificSources_idx = SourcesInd,
                                sink = sinks,
                                all_sources_flag = Shared_all_sources,
                                specific_sources_flag = Shared_specific_sources){
  if(all_sources_flag == 1)
    all_sources_tmp <- count_data[train.ix, ]

  if(all_sources_flag == 0 & specific_sources_flag == 1)
    all_sources_tmp <- count_data[SpecificSources_idx, ]

  if(all_sources_flag == 1 & specific_sources_flag == 1){
    print("Error, please set only one flag to 1 (specific_sources_flag or all_sources_flag)")
    break;
  }
  if(all_sources_flag == 1 & length(all_sources) > 1)
    agg_source_vec <- as.numeric(apply(all_sources_tmp, 2, sum))
  if(all_sources_flag == 1 & length(all_sources) == 1){
    agg_source_vec <- as.matrix(count_data[all_sources, ]) #example train.ix[1]
  }
  if(specific_sources_flag == 1 & length(SpecificSources_idx) == 1)
    agg_source_vec <- as.matrix(count_data[SpecificSources_idx, ]) #example train.ix[1]
  if(specific_sources_flag == 1 & length(SpecificSources_idx) > 1){
    agg_source_vec <- as.numeric(apply(all_sources_tmp, 2, sum))
  }


  Data_all_sources_taxonomy <- data.frame(Taxonomy, t(sink)/sum(sink), agg_source_vec/sum(agg_source_vec))
  names(Data_all_sources_taxonomy) <- c("Kingdom" ,"Phylum", "Class", "Order", "Family", "Genus", "Species", "Sink_relative_abundance" , "Source/s_relative_abundance")
  if(dim(Data_all_sources_taxonomy)[2] == 9)
    return(Data_all_sources_taxonomy)

}


