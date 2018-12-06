library("reshape2")
library("ggplot2")
library("Rcpp")
library("vegan")
library("dplyr")
library("doParallel")
library("foreach")
library("mgcv")
library("cowplot")
library("pROC")
library("plotROC")


 #Plot----

    setwd("/data")
    load("FEAST_healthy.RData")
    load("FEAST_icu.RData")
    
    boxplot(em_noise_healthy$em_noise_results_healthy[,2],
            em_noise_results_icu$em_noise_results_icu[,2])


survived <- c(rep(0, length(na.omit(em_noise_healthy$em_noise_results_healthy[,2]))),
              rep(1, length(na.omit(em_noise_results_icu$em_noise_results_icu[,2]))))
predicted <- c(em_noise_healthy$em_noise_results_healthy[,2], 
               em_noise_results_icu$em_noise_results_icu[,2])
predicted <- na.omit(predicted)


print(auc(survived, predicted))

roc_curves_1 = data.frame(survived, predicted)

ROC_plot <- ggplot() + 
  geom_roc(aes(d = survived, m = predicted, color="FEAST"), roc_curves_1, n.cuts = 0) + 
  geom_line(aes(predicted,predicted), size = 0.2)

ggsave(filename="/results/ICU_ROC.png", plot = ROC_plot , dpi = 600, width = 8.75, height = 6.1, units = "in")


t.test(na.omit(em_noise_healthy$em_noise_results_healthy[,2]), 
       na.omit(em_noise_results_icu$em_noise_results_icu[,2])) #significant

wilcox.test(na.omit(em_noise_healthy$em_noise_results_healthy[,2]), 
            na.omit(em_noise_results_icu$em_noise_results_icu[,2])) #significant