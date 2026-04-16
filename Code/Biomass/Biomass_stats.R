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
 
raw.readbiomass<-read.csv("Data/POSE_Data/SA_POSE_Biomass_08_18_25.csv")

raw.readbiomass2<- raw.readbiomass[-c(69:84),]

### Shoot Mass

for (i in 1:max(raw.readbiomass2$Group, na.rm = TRUE)){
  for (t in unique(raw.readbiomass2$Treatment)) {
    
      raw.readbiomass2$mass_RT_diff[raw.readbiomass2$Group==i & raw.readbiomass2$Treatment==t]<- 
        raw.readbiomass2$AVG_RT_BMSS_g[raw.readbiomass2$Group==i & raw.readbiomass2$Treatment==t] -
        raw.readbiomass2$AVG_RT_BMSS_g[raw.readbiomass2$Group==i & raw.readbiomass2$Treatment=="C"]
      
      raw.readbiomass2$mass_SH_diff[raw.readbiomass2$Group==i & raw.readbiomass2$Treatment==t]<- 
        raw.readbiomass2$AVG_SH_BMSS_g[raw.readbiomass2$Group==i & raw.readbiomass2$Treatment==t] - 
        raw.readbiomass2$AVG_SH_BMSS_g[raw.readbiomass2$Group==i & raw.readbiomass2$Treatment=="C"]
       
  }
}

shootmass_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","mass_SH_diff") #select the columns I want for the germ model
shootmass<-raw.readbiomass2[,shootmass_cols] # make a new dataframe with just those values

null.model<-glm(mass_SH_diff~. ,
                family=gaussian(link = "identity"), ##is identity right?
                data=shootmass)

#fit the global models
global.model <- glm(mass_SH_diff~Treatment+Night_Scale_7+Day_Scale_7+ 
                   Treatment:Night_Scale_7+Treatment:Day_Scale_7+Night_Scale_7:Day_Scale_7,
                   gaussian(link = "identity"),
    data=shootmass)

#run the step function in both directions bounded by the null and global model
models.step<- step(null.model,direction="both",scope = list(lower=null.model, upper=global.model))

summary(models.step)

shootWeightfull<-glm(mass_SH_diff~Treatment+Day_Scale_7+Night_Scale_7+
                       Day_Scale_7:Night_Scale_7,
                  family=gaussian(link = "identity"),data=shootmass)

summary(shootWeightfull)

shootweightrast <- data.frame(summary(shootWeightfull)$coefficients)
shootweightrast$Species <- "POSE" #change to s inside loop
shootweightrast$coefficients <- rownames(shootweightrast)

#Save for later when we have all
#output.all.shootweight <- do.call(rbind, shootweight.output)

shootweightrast$coefficients<-factor(shootweightrast$coefficients,
                                         levels = c("Day_Scale_7:Night_Scale_7",
                                                    "(Intercept)", "TreatmentW",
                                                    "TreatmentL","TreatmentH",
                                                    "Day_Scale_7", "Night_Scale_7"),
                                         labels = c("Day Temperature: Night Temperature",
                                                    "(Intercept)", "Water Treatment", 
                                                    "Low SA Treatment", "High SA Treatment", 
                                                    "Day Temperature", "Night Temperature"))

###is this right??
shootweightrast$sigvar <-
  ifelse(shootweightrast$Pr...t.. < 0.05, 1, 0)

shootweight_sig<-shootweightrast[which(shootweightrast$sigvar==1),]         

ggplot(shootweightrast[-which(shootweightrast$coefficients=="(Intercept)"),],
       aes(x=Species,y=coefficients))+
  geom_raster(aes(fill=Estimate))+
  geom_tile(data =shootweight_sig,
            #[-which(shootweight_sig$coefficients=="(Intercept)"),], 
            aes(x=Species,y=coefficients, fill = sigvar),fill="transparent",
            col="black",size=2)+
  scale_fill_gradient2(low = "red", high = "blue", mid = "white",
                       midpoint = 0, space = "Lab",
                       name="Effect")+
  geom_text(aes(label = paste0(round(Estimate,10))))+
  ylab("Variables & Interactions")+
  xlab("Species")+
  ggtitle("Temperature & Treatments Interactions with Shoot Biomass")+
  theme_grey()

####################
## Root Mass

## Model

rootmass_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","mass_RT_diff") #select the columns I want for the germ model
rootmass<-raw.readbiomass2[,rootmass_cols] # make a new dataframe with just those values

#fit the null model
null.model<-glm(mass_RT_diff~. ,
                family=gaussian(link = "identity"), ##is identity right?
                data=rootmass)

#fit the global models
global.model <- glm(mass_RT_diff~Treatment+Night_Scale_7+Day_Scale_7+ 
                      Treatment:Night_Scale_7+Treatment:Day_Scale_7+Night_Scale_7:Day_Scale_7,
                    gaussian(link = "identity"),
                    data=rootmass)

#run the step function in both directions bounded by the null and global model
models.step<- step(null.model,direction="both",scope = list(lower=null.model, upper=global.model))

summary(models.step)

rootWeightfull<-glm(mass_RT_diff~Treatment+Day_Scale_7+Night_Scale_7+
                       Treatment:Day_Scale_7,
                     family=gaussian(link = "identity"),data=rootmass)

summary(rootWeightfull)

####################
## Raster

rootweightrast <- data.frame(summary(rootWeightfull)$coefficients)
rootweightrast$Species <- "POSE" #change to s inside loop
rootweightrast$coefficients <- rownames(rootweightrast)

#Save for later when we have all
#output.all.shootweight <- do.call(rbind, shootweight.output)

rootweightrast$coefficients<-factor(rootweightrast$coefficients,
                                     levels = c("TreatmentH:Day_Scale_7", "TreatmentL:Day_Scale_7",
                                                "TreatmentW:Day_Scale_7",
                                                "Day_Scale_7:Night_Scale_7",
                                                "(Intercept)", "TreatmentW",
                                                "TreatmentL","TreatmentH",
                                                "Day_Scale_7", "Night_Scale_7"),
                                     labels = c("High SA Treatment: Day Temperature",
                                                "Low SA Treatment: Day Temperature",
                                                "Water Treatment: Day Temperature",
                                                "Day Temperature: Night Temperature",
                                                "(Intercept)", "Water Treatment", 
                                                "Low SA Treatment", "High SA Treatment", 
                                                "Day Temperature", "Night Temperature"))

###is this right??
rootweightrast$sigvar <-
  ifelse(rootweightrast$Pr...t.. < 0.05, 1, 0)

rootweight_sig<-rootweightrast[which(rootweightrast$sigvar==1),]         

ggplot(rootweightrast[-which(rootweightrast$coefficients=="(Intercept)"),],
       aes(x=Species,y=coefficients))+
  geom_raster(aes(fill=Estimate))+
  geom_tile(data =rootweight_sig,
            #[-which(shootweight_sig$coefficients=="(Intercept)"),], 
            aes(x=Species,y=coefficients, fill = sigvar),fill="transparent",
            col="black",size=2)+
  scale_fill_gradient2(low = "red", high = "blue", mid = "white",
                       midpoint = 0, space = "Lab",
                       name="Effect")+
  geom_text(aes(label = paste0(round(Estimate,10))))+
  ylab("Variables & Interactions")+
  xlab("Species")+
  ggtitle("Temperature & Treatments Interactions with Root Biomass")+
  theme_grey()

###########
##Plots of just the raw data
###########

# ggplot(raw.readbiomass2, aes(Treatment, AVG_RT_BMSS_g)) +
#   #geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,size=2)+ ##put the line int he background
#   #geom_smooth(se=F,aes(group=Treatment,col="black"),span = 0.9,size=1)+ ##put the line int he background
#   geom_boxplot() +
#   stat_boxplot(geom ='errorbar',col="black") +
#   stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
#                aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
#   labs( title = "Root Biomass (g) x Treatment",
#         x = "Treatment",
#         y = "Root Biomass")+
#   theme_minimal(base_size = 15)#+
# # scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
# #scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))
# 
# unique(raw.readbiomass2$Treatment)
# 
# raw.readbiomass2$interaction_axis<- log(raw.readbiomass2$Day_Scale_7/raw.readbiomass2$Night_Scale_7)
# 
# raw.readbiomass2$interaction_groups[raw.readbiomass2$interaction_axis >= 0.9737892] <- "2.648~3.909 \n Hot Day, Cold Night"
# raw.readbiomass2$interaction_groups[raw.readbiomass2$interaction_axis >= 0.5842735 & raw.readbiomass2$interaction_axis < 0.9737892] <- "1.794~2.648 \n Warm Day, Cool Night"
# raw.readbiomass2$interaction_groups[raw.readbiomass2$interaction_axis >= 0.1947578 & raw.readbiomass2$interaction_axis < 0.5842735] <- "1.215~1.794 \n Tepid Day, Chilly Night"
# raw.readbiomass2$interaction_groups[raw.readbiomass2$interaction_axis >= -0.1947578 & raw.readbiomass2$interaction_axis < 0.1947578] <- "0.823~1.215 \n Similar Day & Night"
# raw.readbiomass2$interaction_groups[raw.readbiomass2$interaction_axis >= -0.5842735 & raw.readbiomass2$interaction_axis < -0.1947578] <- "0.558~0.823 \n Chilly Day, Tepid Night"
# raw.readbiomass2$interaction_groups[raw.readbiomass2$interaction_axis >= -0.9737892 & raw.readbiomass2$interaction_axis < -0.5842735] <- "0.378~0.558 \n Cool Day, Warm Night"
# raw.readbiomass2$interaction_groups[raw.readbiomass2$interaction_axis < -0.9737892] <- "0.256~0.378 \n Cold Day, Hot Night"
# 
# Temp.int.rmass.box <- ggplot(raw.readbiomass2, aes(x = interaction_groups, y = AVG_RT_BMSS_g, fill = Treatment))+
#   geom_smooth(#method="lm",
#     se=F,aes(group=Treatment,col=Treatment),linewidth=1.5,span=1.5)+ ##put the line int he background
#   geom_boxplot()+
#   stat_boxplot(geom ='errorbar',col="black") +
#   stat_summary(fun = "mean", geom = "point", size = 4, shape = 17,col="black",
#                aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"),
#                show.legend = F)+
#   labs( title = "Average Root Mass (g) x Day Scale:Night Scale",
#         x = "Day Temperature:Night Temperature Ratio",
#         y = "Average Root Mass (g)")+
#   theme_minimal(base_size = 15)+
#   scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
#   scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))
# 
# Temp.int.rmass.box
# 
# ggplot(raw.readbiomass2, aes(Treatment, AVG_SH_BMSS_g)) +
#   #geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,size=2)+ ##put the line int he background
#   #geom_smooth(se=F,aes(group=Treatment,col="black"),span = 0.9,size=1)+ ##put the line int he background
#   geom_boxplot() +
#   stat_boxplot(geom ='errorbar',col="black") +
#   stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
#                aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
#   labs( title = "Shoot Biomass (g) x Treatment",
#         x = "Treatment",
#         y = "Shoot Biomass")+
#   theme_minimal(base_size = 15)
# 
# Temp.int.shmass.box <- ggplot(raw.readbiomass2, aes(x = interaction_groups, y = AVG_SH_BMSS_g, fill = Treatment))+
#   geom_smooth(#method="lm",
#     se=F,aes(group=Treatment,col=Treatment),linewidth=1.5,span=1.5)+ ##put the line int he background
#   geom_boxplot()+
#   stat_boxplot(geom ='errorbar',col="black") +
#   stat_summary(fun = "mean", geom = "point", size = 4, shape = 17,col="black",
#                aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"),
#                show.legend = F)+
#   labs( title = "Average Shoot Mass (g) x Day Scale:Night Scale",
#         x = "Day Temperature:Night Temperature Ratio",
#         y = "Average Shoot Mass (g)")+
#   theme_minimal(base_size = 15)+
#   scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
#   scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))
# 
# Temp.int.shmass.box

#################################
#### PLOTS WITH DIFFERENCE VALUES
#################################

##SHOOT

shootmass2<- shootmass[-c(1:4,65:68),]

ggplot(shootmass2, aes(Treatment, mass_SH_diff)) +
  #geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,size=2)+ ##put the line int he background
  #geom_smooth(se=F,aes(group=Treatment,col="black"),span = 0.9,size=1)+ ##put the line int he background
  geom_boxplot() +
  stat_boxplot(geom ='errorbar',col="black") +
  stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
               aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
  geom_hline(yintercept = 0, lwd= 0.75, linetype = "dashed", col = "blue")+
  labs( title = "Shoot Biomass (g) x Treatment",
        x = "Treatment",
        y = "Shoot Biomass")+
  theme_minimal(base_size = 15)

shootmass2$interaction_axis<- log(shootmass2$Day_Scale_7/shootmass2$Night_Scale_7)

shootmass2$interaction_groups[shootmass2$interaction_axis >= 0.9737892] <- "2.648~3.909 \n Hot Day, Cold Night"
shootmass2$interaction_groups[shootmass2$interaction_axis >= 0.5842735 & shootmass2$interaction_axis < 0.9737892] <- "1.794~2.648 \n Warm Day, Cool Night"
shootmass2$interaction_groups[shootmass2$interaction_axis >= 0.1947578 & shootmass2$interaction_axis < 0.5842735] <- "1.215~1.794 \n Tepid Day, Chilly Night"
shootmass2$interaction_groups[shootmass2$interaction_axis >= -0.1947578 & shootmass2$interaction_axis < 0.1947578] <- "0.823~1.215 \n Similar Day & Night"
shootmass2$interaction_groups[shootmass2$interaction_axis >= -0.5842735 & shootmass2$interaction_axis < -0.1947578] <- "0.558~0.823 \n Chilly Day, Tepid Night"
shootmass2$interaction_groups[shootmass2$interaction_axis >= -0.9737892 & shootmass2$interaction_axis < -0.5842735] <- "0.378~0.558 \n Cool Day, Warm Night"
shootmass2$interaction_groups[shootmass2$interaction_axis < -0.9737892] <- "0.256~0.378 \n Cold Day, Hot Night"

Temp.int.shmass.box <- ggplot(shootmass2, aes(x = interaction_groups, y = mass_SH_diff, fill = Treatment))+
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
  geom_hline(yintercept = 0, lwd= 0.75, linetype = "dashed", col = "blue")+
  scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
  scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))

Temp.int.shmass.box

##ROOT
##NO Sigfnicant interaction w DAY:NIGHT 

rootmass2<- rootmass[-c(1:4,65:68),]

ggplot(rootmass2, aes(Treatment, mass_RT_diff)) +
  #geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,size=2)+ ##put the line int he background
  #geom_smooth(se=F,aes(group=Treatment,col="black"),span = 0.9,size=1)+ ##put the line int he background
  geom_boxplot() +
  geom_hline(yintercept = 0, lwd= 0.75, linetype = "dashed", col = "blue")+
  stat_boxplot(geom ='errorbar',col="black") +
  stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
               aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
  labs( title = "Root Biomass Difference from Control (g) x Treatment",
        x = "Treatment",
        y = "Difference in Root Biomass from Control")+
  theme_minimal(base_size = 15)#+

##Interaction-- not significant
rootmass2$interaction_axis<- log(rootmass2$Day_Scale_7/rootmass2$Night_Scale_7)

rootmass2$interaction_groups[rootmass2$interaction_axis >= 0.9737892] <- "2.648~3.909 \n Hot Day, Cold Night"
rootmass2$interaction_groups[rootmass2$interaction_axis >= 0.5842735 & rootmass2$interaction_axis < 0.9737892] <- "1.794~2.648 \n Warm Day, Cool Night"
rootmass2$interaction_groups[rootmass2$interaction_axis >= 0.1947578 & rootmass2$interaction_axis < 0.5842735] <- "1.215~1.794 \n Tepid Day, Chilly Night"
rootmass2$interaction_groups[rootmass2$interaction_axis >= -0.1947578 & rootmass2$interaction_axis < 0.1947578] <- "0.823~1.215 \n Similar Day & Night"
rootmass2$interaction_groups[rootmass2$interaction_axis >= -0.5842735 & rootmass2$interaction_axis < -0.1947578] <- "0.558~0.823 \n Chilly Day, Tepid Night"
rootmass2$interaction_groups[rootmass2$interaction_axis >= -0.9737892 & rootmass2$interaction_axis < -0.5842735] <- "0.378~0.558 \n Cool Day, Warm Night"
rootmass2$interaction_groups[rootmass2$interaction_axis < -0.9737892] <- "0.256~0.378 \n Cold Day, Hot Night"

Temp.int.rmass.box <- ggplot(rootmass2, aes(x = interaction_groups, y = mass_RT_diff, fill = Treatment))+
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
  geom_hline(yintercept = 0, lwd= 0.75, linetype = "dashed", col = "blue")+
  scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
  scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))

Temp.int.rmass.box

