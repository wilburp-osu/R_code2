### Bio mass
## Clear Workspace
rm(list=ls())
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
library(tidyr)

#could not figure out how to do the subtraction in r, so just did it in excel
#read in data
 
raw.readbiomass<-read.csv("Data/POSE_Data/SA_POSE_Biomass_mathinsheet.csv")

### Shoot Mass

shootmass_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","Corrected_shootweight") #select the columns I want for the germ model
shootmass<-raw.readbiomass[,shootmass_cols] # make a new dataframe with just those values

null.model<-glm(Corrected_shootweight~. ,
                family=gaussian(link = "identity"), ##is identity right?
                data=shootmass)

#fit the global models
global.model <- glm(Corrected_shootweight~Treatment+Night_Scale_7+Day_Scale_7+ 
                   Treatment:Night_Scale_7+Treatment:Day_Scale_7+Night_Scale_7:Day_Scale_7,
                   gaussian(link = "identity"),
    data=shootmass)

#run the step function in both directions bounded by the null and global model
models.step<- step(null.model,direction="both",scope = list(lower=null.model, upper=global.model))

summary(models.step)

shootWeightfull<-glm(Corrected_shootweight~Treatment+Day_Scale_7+Night_Scale_7+
                       Day_Scale_7:Night_Scale_7,
                  family=gaussian(link = "identity"),data=shootmass)

summary(shootWeightfull)

## Root Mass

rootmass_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","Corrected_rootweight") #select the columns I want for the germ model
rootmass<-raw.readbiomass[,rootmass_cols] # make a new dataframe with just those values

#fit the null model
null.model<-glm(Corrected_rootweight~. ,
                family=gaussian(link = "identity"), ##is identity right?
                data=rootmass)

#fit the global models
global.model <- glm(Corrected_rootweight~Treatment+Night_Scale_7+Day_Scale_7+ 
                      Treatment:Night_Scale_7+Treatment:Day_Scale_7+Night_Scale_7:Day_Scale_7,
                    gaussian(link = "identity"),
                    data=rootmass)

#run the step function in both directions bounded by the null and global model
models.step<- step(null.model,direction="both",scope = list(lower=null.model, upper=global.model))

summary(models.step)

ggplot(raw.readbiomass, aes(Treatment, AVG_RT_BMSS_g)) +
  #geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,size=2)+ ##put the line int he background
  #geom_smooth(se=F,aes(group=Treatment,col="black"),span = 0.9,size=1)+ ##put the line int he background
  geom_boxplot() +
  stat_boxplot(geom ='errorbar',col="black") +
  stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
               aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
  labs( title = "Root Biomass (g) x Treatment",
        x = "Treatment",
        y = "Root Biomass")+
  theme_minimal(base_size = 15)#+
# scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
#scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))

raw.readbiomass$interaction_axis<- log(raw.readbiomass$Day_Scale_7/raw.readbiomass$Night_Scale_7)

raw.readbiomass$interaction_groups[raw.readbiomass$interaction_axis >= 0.9737892] <- "2.648~3.909 \n Hot Day, Cold Night"
raw.readbiomass$interaction_groups[raw.readbiomass$interaction_axis >= 0.5842735 & raw.readbiomass$interaction_axis < 0.9737892] <- "1.794~2.648 \n Warm Day, Cool Night"
raw.readbiomass$interaction_groups[raw.readbiomass$interaction_axis >= 0.1947578 & raw.readbiomass$interaction_axis < 0.5842735] <- "1.215~1.794 \n Tepid Day, Chilly Night"
raw.readbiomass$interaction_groups[raw.readbiomass$interaction_axis >= -0.1947578 & raw.readbiomass$interaction_axis < 0.1947578] <- "0.823~1.215 \n Similar Day & Night"
raw.readbiomass$interaction_groups[raw.readbiomass$interaction_axis >= -0.5842735 & raw.readbiomass$interaction_axis < -0.1947578] <- "0.558~0.823 \n Chilly Day, Tepid Night"
raw.readbiomass$interaction_groups[raw.readbiomass$interaction_axis >= -0.9737892 & raw.readbiomass$interaction_axis < -0.5842735] <- "0.378~0.558 \n Cool Day, Warm Night"
raw.readbiomass$interaction_groups[raw.readbiomass$interaction_axis < -0.9737892] <- "0.256~0.378 \n Cold Day, Hot Night"

Temp.int.rmass.box <- ggplot(raw.readbiomass, aes(x = interaction_groups, y = AVG_RT_BMSS_g, fill = Treatment))+
  geom_smooth(#method="lm",
    se=F,aes(group=Treatment,col=Treatment),linewidth=1.5,span=1.5)+ ##put the line int he background
  geom_boxplot()+
  stat_boxplot(geom ='errorbar',col="black") +
  stat_summary(fun = "mean", geom = "point", size = 4, shape = 17,col="black",
               aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"),
               show.legend = F)+
  labs( title = "Average Root Mass (g) x Day Scale:Night Scale",
        x = "Day Temperature:Night Temperature Ratio",
        y = "Average Root Mass (g)")+
  theme_minimal(base_size = 15)+
  scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
  scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))

Temp.int.rmass.box

ggplot(raw.readbiomass, aes(Treatment, AVG_SH_BMSS_g)) +
  #geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,size=2)+ ##put the line int he background
  #geom_smooth(se=F,aes(group=Treatment,col="black"),span = 0.9,size=1)+ ##put the line int he background
  geom_boxplot() +
  stat_boxplot(geom ='errorbar',col="black") +
  stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
               aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
  labs( title = "Shoot Biomass (g) x Treatment",
        x = "Treatment",
        y = "Shoot Biomass")+
  theme_minimal(base_size = 15)

Temp.int.shmass.box <- ggplot(raw.readbiomass, aes(x = interaction_groups, y = AVG_SH_BMSS_g, fill = Treatment))+
  geom_smooth(#method="lm",
    se=F,aes(group=Treatment,col=Treatment),linewidth=1.5,span=1.5)+ ##put the line int he background
  geom_boxplot()+
  stat_boxplot(geom ='errorbar',col="black") +
  stat_summary(fun = "mean", geom = "point", size = 4, shape = 17,col="black",
               aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"),
               show.legend = F)+
  labs( title = "Average Shoot Mass (g) x Day Scale:Night Scale",
        x = "Day Temperature:Night Temperature Ratio",
        y = "Average Shoot Mass (g)")+
  theme_minimal(base_size = 15)+
  scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
  scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))

Temp.int.shmass.box

# #read in data
# 
# raw.readbiomass<-read.csv("Data/POSE_Data/SA_POSE_Biomass_08_18_25.csv")
# 
# glimpse(raw.readbiomass)
# 
# raw.readbiomass$unique_ID<-paste0(raw.readbiomass$Group,"_",raw.readbiomass$Treatment)
# 
# raw.readbiomass
# 
# raw.readbiomass %>%
#   mutate(across(where(is.character), as.factor))
# 
# raw.readbiomass$Shoot_difference <- raw.readbiomass %>%
#   group_by(Group) %>%
#   mutate(across(raw.readbiomass$Treatment, ~ raw.readbiomass$AVG_SH_BMSS_g-
#                   raw.readbiomass$AVG_SH_BMSS_g[raw.readbiomass$Treatment == "C"]))
#   
#   
# 
# if (raw.readbiomass$Group == 1) {
#   raw.readbiomass$corrected_shootbio <- raw.readbiomass$AVG_SH_BMSS_g-
#    raw.readbiomass$AVG_SH_BMSS_g[raw.readbiomass$Treatment == C]
#   }
#     


# remove all columns not needed
# what columns we need are... c("Group","days_til_emerg",
#         "days_til_germ","Treatment","Day_Scale_7","Night_Scale_7","unique_ID)


##Take out each cuvette value
# 
# raw.readbiomass$Treatment
# 
# ##Subtract cuvette from control
# 
# ##Put into A new column
# 
# for (i in 1:max(temp$Group)){
#   for (t in unique(temp$Treatment)) {
#     for (c in 1:4) {
#       temp_seed<-temp[which(temp$Group==i & temp$Treatment==t & temp$Cell==c),]
#       if (any(temp_seed$Root_Pheno_Code==1)){
#         Days<-temp_seed$Days_Since_Install[which(temp_seed$Root_Pheno_Code==1)]
#         if (length(Days)>1){
#           print(paste0(c("1s",s,Days,i,t,c),collapse="_"))} #print errors for doubling up 1's
#         temp$germ_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-1/Days
#         temp$days_til_germ[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-Days
#       } else {
#         temp$germ_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-0
#         temp$days_til_germ[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-0
#       }
#       if (any(temp_seed$Shoot_Pheno_Code==6)){
#         Days2<-temp_seed$Days_Since_Install[which(temp_seed$Shoot_Pheno_Code==6)]
#         if (length(Days2)>1){
#           print(paste0(c("6s",s,Days2,i,t,c),collapse="_"))} #print errors for doubling up 6's
#         temp$emerg_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-1/Days2  
#         temp$days_til_emerg[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-Days2  
#       } else {
#         temp$emerg_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-0
#         temp$days_til_emerg[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-0
#       }
#     }
#   }
# }
# 
# Above.Biomass.cols<-c("Treatment", "AVG_SH_BMSS_g") #select the columns I want for the germ model
# temp_agg2_germ<-temp_agg2[,germ_cols] # make a new dataframe with just those values
# 
# 
# 
# #fit the null model
# null.model<-glm(days_til_germ~. ,
#                 family=poisson(link = "log"),data=temp_agg2_germ)

# 
# biomass.temp <- 
#   
# temp_agg2 <- short.data.list2[[s]] # pull out the data
# germ_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","days_til_germ") #select the columns I want for the germ model
# temp_agg2_germ<-temp_agg2[,germ_cols] # make a new dataframe with just those values
# 
# 
# 
# #fit the null model
# null.model<-glm(AVG_RT_BMSS_g~. ,
#                 family=gaussian(link = "log"),data=raw.readbiomass)
# 
# ###Why would this not work?
# #fit the global models
# glm(AVG_RT_BMSS_g~Treatment + Night_Scale_7 + Day_Scale_7,
#         family=gaussian(link = "identity"),data=raw.readbiomass)
# 
# #run the step function in both directions bounded by the null and global model
# models.step<- step(null.model,direction="both",scope = list(lower=null.model, upper=global.model))
# 
# print(summary(models.step))