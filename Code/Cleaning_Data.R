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

#### read in POSE data ####

POSE<-read.csv("Data/POSE_Data/SA_POSE_Data_08_18_25.csv")

str(POSE)

POSE$Uni_ID<-paste0(POSE$Group,"_",POSE$Treatment,"_",POSE$Cell)
str(POSE)

for (i in 1:49){
seed_plot<-ggplot(POSE[which(POSE$Group==i),],
                  aes(x=Days_Since_Install,y=Shoot_Pheno_Code,col=Uni_ID))+
  geom_line()+
  geom_point(size=3)+
  geom_hline(yintercept = 6)+
  facet_grid(.~Uni_ID)+
  theme(legend.position="none")+
  ylim(c(0,10))
print(seed_plot)
}


#### read in ACMI data ####
ACMI<-read.csv("Data/ACMI_Data/SA_ACMI_Data_08_18_25.csv")

str(ACMI)

ACMI$Uni_ID<-paste0(ACMI$Group,"_",ACMI$Treatment,"_",ACMI$Cell)
str(ACMI)

for (i in 1:49){
  seed_plot<-ggplot(ACMI[which(ACMI$Group==i),],
                    aes(x=Days_Since_Install,y=Shoot_Pheno_Code,col=Uni_ID))+
    geom_line()+
    geom_point(size=3)+
    geom_hline(yintercept = 6)+
    facet_grid(.~Uni_ID)+
    theme(legend.position="none")+
    ylim(c(0,10))
  print(seed_plot)
}
