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

data.list<-list()
species<-c("POSE","ARTR","ACMI","ELEL")

for (s in species ) {
  data.list[[s]]<-read.csv(paste0("Data/",s,"_Data/SA_",s,"_Data_FINAL.csv"))
}

str(data.list)
clean.data.list<-list()

for (s in species) {
  temp<-data.list[[s]]
  ### clean data ####
  temp$Root_Pheno_Code<-as.numeric(paste0(temp$Root_Pheno_Code))
  temp$Shoot_Pheno_Code<-as.numeric(paste0(temp$Shoot_Pheno_Code))
  
  temp$Root_Pheno_Code[which(temp$Root_Pheno_Code==999)]<-NA
  temp$Shoot_Pheno_Code[which(temp$Shoot_Pheno_Code==999)]<-NA
  
  temp$Avg_Root_rate<-NA
  temp$Avg_Shoot_rate<-NA
  
  ## Calculate growth rates
  Root_Code<-c(1,2,3,4)
  Root_Length_mid<-c(0.5,5,10,20)
  Root_match<-data.frame(Code=Root_Code,Length=Root_Length_mid)
  temp$Root_Length<-Root_match$Length[match(temp$Root_Pheno_Code,Root_match$Code)]
  
  Shoot_Code<-c(6,7,8)
  Shoot_Length_mid<-c(0.5,10,20)
  Shoot_match<-data.frame(Code=Shoot_Code,Length=Shoot_Length_mid)
  temp$Shoot_Length<-Shoot_match$Length[match(temp$Shoot_Pheno_Code,Shoot_match$Code)]
  
  for (i in 1:max(temp$Group)){
    for (t in unique(temp$Treatment)) {
      for (c in 1:4) {
        temp_seed<-temp[which(temp$Group==i & temp$Treatment==t & temp$Cell==c),]
        if (any(temp_seed$Root_Pheno_Code==1)){
          root_days<-temp_seed$Root_Length[which(temp_seed$Root_Pheno_Code>0&temp_seed$Root_Pheno_Code<5)]
          #root_days_dup<-unique(root_days)
          #root_rate<-sum(root_days_dup)/length(root_days)
          root_rate<-max(root_days)/length(root_days)
        } else {root_rate<-0}
        if (any(temp_seed$Shoot_Pheno_Code==6)){
          shoot_days<-(temp_seed$Shoot_Length[which(temp_seed$Shoot_Pheno_Code>5&temp_seed$Shoot_Pheno_Code<9)])
          #shoot_days_dup<-unique(shoot_days)
          #shoot_rate<-sum(shoot_days_dup)/length(shoot_days)
          shoot_rate<-max(shoot_days)/length(shoot_days)
        } else {shoot_rate<-0}
        temp$Avg_Root_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-root_rate
        temp$Avg_Shoot_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-shoot_rate
      }
    }
  }
  
  clean.data.list[[s]]<-temp
}

#### read in data ####
#POSE<-read.csv("Data/POSE_Data/SA_POSE_Data_01_08_25.csv")
#days<-unique(POSE$Day_Collected)
