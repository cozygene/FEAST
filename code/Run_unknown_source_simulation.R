
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



source("SourceTracker.R")
source("src.R")


SAMPLINGITERATIONS=4
NUMSOURCES=4
BODYPART=T
HOUSEONLY=T
POOLBODYPARTS=T
COVERAGE=10000
EM_ITERATIONS<-100
MIXINGITERATIONS<-10
USESTEP=F
include_epsilon <- T
STEPS=4


metadata<-read.table(paste0("../data/house1_t", 1,"_metadata.tsv"),
                     sep= '\t', header=T, comment='', stringsAsFactors = F)
idx<-apply(metadata, 1, function(var) all(var != ''))
metadata<-metadata[idx,]
rownames(metadata)<-metadata[,1]
metadata<-metadata[,-1]                 
otus<-read.table(paste0("../data/house1_t", 1,"_otus.txt"),
                 comment='', skip=1, header=T, sep='\t', check=F )

houseidx<-as.numeric(grep(pattern='2192', row.names(metadata)))
metadata<-metadata[houseidx,]
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
}
metadata<-metadata[which(metadata$SourceSink == 'source'),]


rownames(otus)<-otus[,1]
otus<-otus[,-1]
otus<-t(otus)
common_names<-intersect(rownames(metadata), rownames(otus))
metadata<-metadata[common_names,]
otus<-otus[common_names,]
sources<-rownames(metadata)[which(metadata$SourceSink == "source")]
envs<-as.character(metadata[sources, 2])
sources<-otus[sources,]
holdsources<-sources

eps_list = seq(0,1,0.05)

for(eps in eps_list){

print(eps)
  srcidx<-sample(dim(holdsources)[1], NUMSOURCES)
  sources<-holdsources[srcidx,]
  sources<-rarefy(sources, maxdepth = COVERAGE)
  envs<-metadata$Env[srcidx]
  unk_source<-rarefy(t(holdsources[sample(seq(1,dim(holdsources)[1])[-srcidx], 1),]), maxdepth=COVERAGE)
  ms<-list()
  sinks<-list()

  for(m in 1:MIXINGITERATIONS){
    ms[[m]]<-runif(NUMSOURCES, 0, 1)
    ms[[m]]<-c(ms[[m]], eps*sum(ms[[m]])/(1-eps)) / sum(c(ms[[m]], eps*sum(ms[[m]])/(1-eps)))
    sinks[[m]]<-round(t(ms[[m]]) %*% rbind(sources,unk_source))
  }

  st<-sourcetracker(sources, envs, COVERAGE)
  sources<-st$sources[-which(rownames(st$sources) == 'Unknown'),]
  sources<-rarefy(sources, COVERAGE)
  js_values<-c()
  x <- sources[c(1:dim(sources)[1]),]
  JSDMatrix <- jsdmatrix(x)
  JSDMatrix <- JSDMatrix/apply(sources, 1, sum)[1]
  JS = mean(JSDMatrix[-which(JSDMatrix == 0)])
  js_values = append(js_values, JS)
  print(js_values)
  st_df<-data.frame()
  feast_df<-data.frame()
  for(i in 1:length(sinks)){
    sink<-rarefy(sinks[[i]], apply(sources, 1, sum)[1])
    
    ####  Run SourceTracker ####
    if(ST_unknown_flag == 1){
        
        results <- predict(st,sink)
        st_df<-rbind(st_df, results$proportions)
        print("st")
        print(results$proportions)
    
    }
    
    ####  Initialize the unknown source using raw sources and sink ####
    unk_source<-unknown__initialize_1(sources, sink, n_sources=dim(sources)[1])
    unk_source<-rarefy(t(unk_source), apply(sources, 1, sum)[1])
    em_source<-rbind(sources, unk_source)
    labels<-c(levels(st$envs), "Unknown")
    obs_source<-st$sources
    
    
    initalphs<-runif(dim(em_source)[1], min = 0.0, max = 1.0)
    initalphs <- initalphs / Reduce('+', initalphs)
  
    
    ####  Run EM  ####
    em_source<-split(em_source, seq(nrow(em_source)))
    obs_source<-split(obs_source, seq(nrow(obs_source)))
    em_source<-lapply(em_source,t)
    obs_source<-lapply(obs_source,t)
    tmp<-do_EM(alphas=initalphs, sources=em_source, 
               observed=obs_source, sink = sink, iterations = EM_ITERATIONS)
    names(tmp$toret)<-labels
    feast_df<-rbind(feast_df, tmp$toret)
    print("feast")
    print(tmp$toret)
  
  }
  
  print(getwd())



}
