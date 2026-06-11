########################################
rm(list=ls()) # clear workspace
# read in libraries ####
library(terra)
library(tidyverse)
library(patchwork)

######################################
#create open list
data.list<-list()
#create species names
species<-c("POSE","ARTR","ACMI","ELEL")

# Read in data ####
for (s in species ) {
  raw.read<-read.csv(paste0("Data/",s,"_Data/",s,"_Data_FINAL.csv"))
  if (any(is.na(raw.read$Group))){
    raw.read2<-raw.read[-which(is.na(raw.read$Group)),]
  } else {raw.read2<-raw.read}
  data.list[[s]]<-raw.read2
}


source("code/Final_Code/Days_til_Germ_PMW.R")
source("code/Final_Code/Days_til_Emerg_PMW.R")
source("code/Final_Code/Shoot_Growth_Rate_PMW.R")
source("code/Final_Code/Root_Growth_Rate_PMW.R")
source("code/Final_Code/BoxPlots_PMW.R")
source("code/Final_Code/Interaction_Plots_PMW.R")
source("code/Final_Code/Fiftee_Twentyfive_box.R")

Root_Rasters <- Germination_Rasters + Root_Growth_Raster + 
  plot_layout(ncol=2, axis_titles = 'collect_y', axes = 'collect_y')+
  plot_annotation(tag_levels = 'a')

ggsave("Figures/Side_by_Root.png", Root_Rasters)

Shoot_Rasters <- Emergence_Rasters + Shoot_Growth_Raster + 
  plot_layout(ncol=2, axis_titles = 'collect_y', axes = 'collect_y')+
  plot_annotation(tag_levels = 'a')

ggsave("Figures/Side_by_Shoot.png", Shoot_Rasters)

DayGerm_DayEmerg <- Germination_Rasters + Emergence_Rasters +
  plot_layout(ncol= 2, axis_titles = 'collect_y', axes = 'collect_y')+ 
  plot_annotation(tag_levels = 'a')

GRate_ERate <- Root_Growth_Raster + Shoot_Growth_Raster +
  plot_layout(ncol= 2, axis_titles = 'collect_y', axes = 'collect_y')+ 
  plot_annotation(tag_levels = 'a')

ggsave("Figures/Root_by_shoot_days.png", DayGerm_DayEmerg)

ggsave("Figures/Root_by_shoot_rate.png", GRate_ERate)
