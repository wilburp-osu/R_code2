######################################
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

source("code/forloops/for_loops_growth_rate_TEM.R")
source("code/forloops/testing_models_zero_nbimo_TEM.R")