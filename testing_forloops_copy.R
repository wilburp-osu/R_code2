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

###Peter Uncomment to set WD
#setwd("C:/Users/18034/Dropbox/PC/Desktop/Martyn_Lab/Peter_Thermodatatable_entry/4_SA")

#create open list
data.list<-list()
#create species names
species<-c("POSE","ARTR","ACMI","ELEL")

###Read in data
for (s in species ) {
  data.list[[s]]<-read.csv(paste0("Data/",s,"_Data/",s,"_Data_FINAL.csv"))
}

####GROWTH RATES####

###Join Strings
str(data.list)
clean.data.list<-list()

for (s in species) {
  temp<-data.list[[s]]
  ### clean data ####
  ##convert codes to numeric type and save in temp$R_P_C
  temp$Root_Pheno_Code<-as.numeric(paste0(temp$Root_Pheno_Code))
  temp$Shoot_Pheno_Code<-as.numeric(paste0(temp$Shoot_Pheno_Code))
  
  ##make NAs
  temp$Root_Pheno_Code[which(temp$Root_Pheno_Code==999)]<-NA
  temp$Shoot_Pheno_Code[which(temp$Shoot_Pheno_Code==999)]<-NA
  
  ##Make avg root rates NAs
  temp$Avg_Root_rate<-NA
  temp$Avg_Shoot_rate<-NA
  
  ## Calculate growth rates
  ##create root code vector
  Root_Code<-c(1,2,3,4) ##Why no 5?##
  ## create root length mid vector 
  Root_Length_mid<-c(0.5,5,10,20)
  #combine vectors
  Root_match<-data.frame(Code=Root_Code,Length=Root_Length_mid)
  ##match root pheno codes to newly created codes
  temp$Root_Length<-Root_match$Length[match(temp$Root_Pheno_Code,Root_match$Code)]
  
  ##create shoot code vector
  Shoot_Code<-c(6,7,8) ##9?##
  ## create shoot length mid vector
  Shoot_Length_mid<-c(0.5,10,20)
  ##combine vectors
  Shoot_match<-data.frame(Code=Shoot_Code,Length=Shoot_Length_mid)
  ##match shoot pheno codes to new codes
  temp$Shoot_Length<-Shoot_match$Length[match(temp$Shoot_Pheno_Code,Shoot_match$Code)]
  
  ##open group for loop
  for (i in 1:max(temp$Group)){
    ##open treatment for loop
    for (t in unique(temp$Treatment)) {
      ##Open seed for loop
      for (c in 1:4) {
        #create data set for each seed
        temp_seed<-temp[which(temp$Group==i & temp$Treatment==t & temp$Cell==c),]
        if (any(temp_seed$Root_Pheno_Code==1, na.rm = TRUE)){ ##only way I could get it to run##
          root_days<-temp_seed$Root_Length[which(temp_seed$Root_Pheno_Code>0&temp_seed$Root_Pheno_Code<5)]
          #root_days_dup<-unique(root_days)
          #root_rate<-sum(root_days_dup)/length(root_days)
          root_rate<-max(root_days)/length(root_days)
        } else {root_rate<-0}
        if (any(temp_seed$Shoot_Pheno_Code==6, na.rm = TRUE)){ ##Ditto##
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

#Still getting these errors/ warnings
#Error in 1:max(temp$Group) : NA/NaN argument
#In addition: Warning messages:
#1: In as.numeric(paste0(temp$Root_Pheno_Code)) :
  #NAs introduced by coercion
#2: In as.numeric(paste0(temp$Shoot_Pheno_Code)) :
  #NAs introduced by coercion

#### read in data ####
#POSE<-read.csv("Data/POSE_Data/POSE_Data_FINAL.csv")
#days<-unique(POSE$Day_Collected)
