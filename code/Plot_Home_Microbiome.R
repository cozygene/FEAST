library("Rcpp")
library("vegan")
library("dplyr")
library("doParallel")
library("foreach")
library("mgcv")
library("reshape2")
library("ggplot2")
library("cowplot")
library("RColorBrewer")
library("RcppArmadillo")
cppFunction("arma::mat schur(arma::mat& a, arma::mat& b) 
            {return(a % b); }", depends="RcppArmadillo")


set.seed(42)
dataset = "FEAST_Results_house1_10000_init1_nofilter"

if(grepl("FEAST",dataset)){
  tool = "FEAST"
}else{
  tool = "ST"
}
if(grepl("house1",dataset)){
  house_num = 1
}
if(grepl("house7b",dataset)){
  house_num = "7b"
}
if(grepl("house7_b",dataset)){
  house_num = "7_b"
}

num_times = 14
BODYPARTS = T
color_chunks = T

ST_sinks <- c("kitchen_counter","bathroom_door_knob","front_door_knob","kitchen_light_switch","bedroom_floor","kitchen_floor")
FEAST_sinks <- c("KitchenCounter","BathroomDoorKnob","FrontDoorKnob","KitchenLightSwitch","BedroomFloor","KitchenFloor")
my_sinks <- c("Kitchen Counter")


if(BODYPARTS){
  if(grepl("house1_",dataset)){
    person_label <- c("person1_1_Foot","person1_1_Hand","person1_1_Nose","person2_1_Foot","person2_1_Hand", "person2_1_Nose","person3_1_Foot",
                      "person3_1_Hand","person3_1_Nose","person4_1_Foot","person4_1_Hand","person4_1_Nose","Unknown")
    sources_names <- c("person1_1_Foot","person1_1_Hand","person1_1_Nose","person2_1_Foot","person2_1_Hand", "person2_1_Nose","person3_1_Foot",
                       "person3_1_Hand","person3_1_Nose","person4_1_Foot","person4_1_Hand","person4_1_Nose","Unknown")
  }
  if(grepl("house7b_",dataset)){
    person_label <- c("cat_7_Fur","cat_7_Paw","person1_7_Foot","person1_7_Hand","person1_7_Nose", 
                      "person2_7_Foot","person2_7_Hand","person2_7_Nose","Unknown" )
    sources_names <- c("cat_7_Fur","cat_7_Paw","person1_7_Foot","person1_7_Hand","person1_7_Nose", 
                       "person2_7_Foot","person2_7_Hand","person2_7_Nose","Unknown" )
  }
}else{
  person_label <- c("person1","person2","person3","person4","unknown")
  sources_names <- c("person 1","person 2","person 3","person 4","unknown")
  
}
num_sources = length(sources_names)


getpal<-colorRampPalette(brewer.pal(8,"Set3"))
cols<-getpal(length(sources_names))
cols<-c(cols[-length(cols)],"#556D7C") #"#40525E"

if(color_chunks == T){
  if(grepl("house1_",dataset)){
    cols<-c(brewer.pal(4,"Greens"),brewer.pal(4,"Purples"),brewer.pal(4,"Oranges"))
  }
  if(grepl("house7b_",dataset)){
    #cols<-c(brewer.pal(2,"Greens"),brewer.pal(2,"Purples"),brewer.pal(2,"Oranges"),brewer.pal(1,"Reds"))
    cols<-brewer.pal(8,"YlOrBr")
    #cols<-c(cols,"#2c7fb8")
  }
  
  cols<-c(cols,"#556D7C") #"#40525E"
}

setwd("/data/Home_Microbiome_house1/")

for(s in 1:length(my_sinks)){

  my_sink_frame <- matrix(0,num_times,num_sources)
  colnames(my_sink_frame) = sources_names
  for(t in 1:num_times){
    print(paste0("time",t))
    

    if(tool == "ST"){
      
      ls = list.files()
      my_file = ls[grep(paste0("House",house_num,"_t_",t,"_",ST_sinks[s],".txt"),ls)]
    }
    
    if(tool == "FEAST"){
    
      
      ls = list.files()
      my_file = ls[grep(paste0("^em_baseline_unkinit.*",FEAST_sinks[s],".*t",t,"_.txt"), ls)]
    }
    
    if(length(my_file)>0){
      input <- unlist(read.table(my_file,sep="\t",header=T,row.names=1,check=FALSE))
      print(input)
    }else{
      break
    }
    for(i in 1:(num_sources-1)){
      if( length(grep(person_label[i],names(input))) > 0){
        my_sink_frame[t,i] <- input[grep(person_label[i],names(input))]
      }
      
    }
    my_sink_frame[t,ncol(my_sink_frame)] <- input["Unknown"]
    
  }
  if(color_chunks == T){
    # desired_order <- c("person1_1_Foot","person1_1_Hand","person1_1_Nose","person2_1_Foot","person2_1_Hand", "person2_1_Nose","person3_1_Foot",
    #                   "person3_1_Hand","person3_1_Nose","person4_1_Foot","person4_1_Hand","person4_1_Nose","Unknown")
    if(grepl("house1_",dataset)){
      new_order = c(grep("Nose",person_label),grep("Hand",person_label),grep("Foot",person_label),ncol(my_sink_frame))
      
    }
    if(grepl("house7b_",dataset)){
      new_order = c(grep("Nose",person_label),grep("Hand",person_label),grep("Foot",person_label),grep("cat",person_label),ncol(my_sink_frame))
      
    }
    
    my_sink_frame <- my_sink_frame[,new_order]
  }
  
  
  my_sink_frame <- my_sink_frame[,new_order]
  my_melt <- melt(my_sink_frame)
  my_melt$Var2 <- as.character(my_melt$Var2)
  
  Humans = rep("Humans", length(which(my_melt$Var2 != "Unknown")))
  
  my_melt$Var2[my_melt$Var2 != "Unknown"] = Humans
  my_melt$Var2 = as.factor(my_melt$Var2)
  
  my_plot <- ggplot(my_melt, aes(x=Var1,y=value,group=Var2,fill=Var2)) + 
    geom_bar(stat="identity") +
    guides(fill=guide_legend(title="Source")) +
    xlab("Time point") +
    scale_x_continuous(breaks=1:num_times) +
    ylab("Proportion") +
    theme(plot.title = element_text(hjust=1)) +
    scale_fill_manual(values = cols, guide = guide_legend(nrow=1)) +
    ggtitle(my_sinks[s]) + theme_bw()

    ggsave(filename="/results/Home_microbiome_analysis.png", plot = my_plot , 
    dpi = 600, width = 8.75, height = 6.1, units = "in")

}



