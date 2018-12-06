library("reshape2")
library("ggplot2")
library("Rcpp")
library("vegan")
library("dplyr")
library("doParallel")
library("foreach")
library("mgcv")
library("cowplot")

 #Plot----
  setwd("../data/unknown_simulations_scatter_mse/")
  ls<-list.files()
  feasts<-ls[grep("unkinit_1_feast.*", ls)]
  sts<-ls[grep("unkinit_1_st.*", ls)]
  ms_<-ls[grep("unkinit_1_ms.*", ls)]
  # JSD: .5889564
  feastunks<-c()
  stunks<-c()
  msunks<-c()
  for(ffile in feasts){
    x<-readRDS(ffile)
    unk<-mean(x[,dim(x)[2]])
    feastunks<-c(feastunks, unk)
  }
  for(sfile in sts){
    x<-readRDS(sfile)
    unk<-mean(x[,dim(x)[2]])
    stunks<-c(stunks, unk)
  }
  for(mfile in ms_){
    x<-readRDS(mfile)
    unk<-mean(unlist(lapply(x, function(z) z[[length(z)]])))
    msunks<-c(msunks, unk)
  }

  toplot<-data.frame(cbind(msunks, feastunks, stunks))
  scatplt<-ggplot(toplot, aes(msunks)) + geom_point(aes(y=feastunks, colour=variable), colour="#E69F00", size=2.5,
                                                    show.legend = T)
  scatplt = scatplt + geom_point(aes(y= stunks, colour=variable), colour="#0072B2", size=2.5, 
                                 show.legend = T)
  scatplt = scatplt  + geom_abline(intercept=0, slope=1, alpha = .3, size = 1)
  scatplt = scatplt + scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits=c(0,1))
  scatplt = scatplt + theme(legend.position="bottom")
  scatplt = scatplt + xlab("True unknown source proportion") + ylab("Predicted unknown source proportion")
  scatplt = scatplt + theme(axis.title=element_text(size=14)) + theme(axis.text=element_text(size=12))
  scatplt = scatplt + ggtitle("Unknown source predicted mixing proportions")
  scatplt = scatplt + theme_minimal()
  scatplt = scatplt + theme(plot.title = element_text(hjust = 0.5,size=16))
  scatplt = scatplt + scale_colour_manual(name="Method", values=c("#E69F00", "#0072B2"), 
                                          breaks=c("FEAST","SourceTracker"), labels=c("FEAST","SourceTracker"))
  names(toplot)<-c("Ms", "FEAST", "SourceTracker")
  toplot.melt <- melt(toplot, id.vars = "Ms", variable.name= "Method",
                      value.name = "Predicted")
  newplot <- ggplot(toplot.melt, aes(x = Ms, y = Predicted, color = Method)) + geom_point(size=2.5)
  newplot = newplot + theme_bw()
  newplot = newplot  + geom_abline(intercept=0, slope=1, alpha = .3, size = 1)
  newplot = newplot + scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits=c(0,1))
  newplot = newplot + xlab("True unknown source proportion") + ylab("Predicted unknown source proportion")
  newplot = newplot + ggtitle("Unknown source predicted mixing proportions")
  newplot = newplot + theme(plot.title = element_text(hjust = 0.5,size=20))
  newplot = newplot + theme(axis.title=element_text(size=16)) + theme(axis.text=element_text(size=16))
  newplot = newplot + theme(legend.title=element_text(size=16), legend.text=element_text(size=15)) 

  ggsave(filename="/results/Unknown_source_predicted_mixing_proportions.png", plot = newplot , dpi = 600, width = 8.75, height = 6.1, units = "in")

  
