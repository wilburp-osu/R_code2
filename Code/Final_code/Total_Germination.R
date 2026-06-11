rm(list=ls()) # clear workspace
# read in libraries ####
library(raster)
library(terra)
library(lme4)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(tidyverse)
library(plotly)
library(Cairo)
library(dplyr)

######################################
#create open list
data.list<-list()
#create species names
species<-c("POSE","ARTR","ACMI","ELEL")

##set wd
setwd("C:/Users/18034/Dropbox/PC/Desktop/R_code/R_code2")

# Read in data ####
for (s in species ) {
  raw.read<-read.csv(paste0("Data/",s,"_Data/",s,"_Data_FINAL.csv"))
  if (any(is.na(raw.read$Group))){
    raw.read2<-raw.read[-which(is.na(raw.read$Group)),]
  } else {raw.read2<-raw.read}
  data.list[[s]]<-raw.read2
}

# separate data into just day 14

final.day.list <- list()

for (s in species ) {
  full.day <- data.list[[s]]
  
  final.day <- full.day[which(full.day$Days_Since_Install==14),]
  
  final.day.list[[s]] <-final.day
}

# create a binary for yes or no germination/ emergence

Total.germ.list <-list()

for (s in species){ 
binary.start <- final.day.list[[s]]
  
binary.start$germ_yn<- ifelse(binary.start$Root_Pheno_Code < 0.5, 0, 1)
binary.start$emerg_yn<- ifelse(binary.start$Shoot_Pheno_Code < 0.5, 0, 1)

Total.germ.list[[s]] <- binary.start

}

# graph

plots.list <- list()

for (s in species){
  germ_agg <- Total.germ.list[[s]] 
  
  total.germ.boxplots <- ggplot(germ_agg, aes(Treatment, germ_yn,fill=Treatment)) +
    geom_bar(stat = "identity") +
    #geom_errorbar(aes(ymin=germ_yn, ymax=germ_yn+sd), width=.2,
          #                   position=position_dodge(.9))+
    labs( title = paste0("Total Germination x Treatment (", s, ")"),
          x = "Treatment",
          y = "Total Germination")+
    theme_minimal(base_size = 15)+
    scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))
  
  print(total.germ.boxplots)
  
  plots.list[[s]] <- total.germ.boxplots
  
}

(plots.list[["POSE"]] + plots.list[["ACMI"]]) +
  (plots.list[["ARTR"]] + plots.list[["ELEL"]])

Total_germ.plotsall <- plots.list[["POSE"]] + plots.list[["ACMI"]] + 
  plots.list[["ARTR"]] + plots.list[["ELEL"]]+
  plot_layout(ncol= 2, axis_titles = 'collect', guides = 'collect')&
  theme(legend.position="bottom")
# run stats
