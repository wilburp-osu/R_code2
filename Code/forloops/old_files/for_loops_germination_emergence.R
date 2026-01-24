rm(list=ls())
#### read in libraries ####
library(raster)
library(terra)
library(lme4)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(tidyverse)
library(plotly)

#create open list
data.list<-list()
#create species names
species<-c("POSE","ARTR","ACMI","ELEL")

###Read in data
for (s in species ) {
  data.list[[s]]<-read.csv(paste0("Data/",s,"_Data/",s,"_Data_FINAL.csv"))
}

#days<-unique(POSE$Day_Collected)

str(data.list)
clean.data.list<-list()

for (s in species){

### clean data ####
POSE$Root_Pheno_Code<-as.numeric(paste0(POSE$Root_Pheno_Code))
POSE$Shoot_Pheno_Code<-as.numeric(paste0(POSE$Shoot_Pheno_Code))

temp$germ_rate<-NA
temp$emerg_rate<-NA
####1
### calculate days til emergence and germination

for (i in 1:max(temp$Group)){
  for (t in unique(temp$Treatment)) {
    for (c in 1:4) {
      temp_seed<-temp[which(temp$Group==i & temp$Treatment==t & temp$Cell==c),]
      if (any(temp_seed$Root_Pheno_Code==1)){
        Days<-temp_seed$Days_Since_Install[which(temp_seed$Root_Pheno_Code==1)]
        temp$germ_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-1/Days
      } else {
        temp$germ_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-0
      }
      if (any(temp_seed$Shoot_Pheno_Code==6)){
        Days2<-temp_seed$Days_Since_Install[which(temp_seed$Shoot_Pheno_Code==6)]
        temp$emerg_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-1/Days2  
      } else {
        temp$emerg_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-0
      }
    }
  }
}

}
