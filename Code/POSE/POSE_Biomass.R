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

#### read in data ####
POSE.bmass<-read.csv("Data/POSE_Data/SA_POSE_Biomass.csv")
POSE<-read.csv("Data/POSE_Data/SA_POSE_Data_TEM.csv")
days<-unique(POSE.bmass$Day_Collected)

ggplot(POSE.bmass, aes(y=AVG_RT_BMSS_g,x=Treatment))+
  geom_boxplot()
ggplot(POSE.bmass, aes(y=AVG_SH_BMSS_g,x=Treatment))+
  geom_boxplot()
ggplot(POSE.bmass, aes(y=Root_Biomass_g,x=Treatment))+
  geom_boxplot()
ggplot(POSE.bmass, aes(y=Shoot_Biomass_g,x=Treatment))+
  geom_boxplot()

## add temperature to dataset
POSE.temp<-POSE[,c(15,17,19,20,21)]
POSE.bmass2<-merge(POSE.bmass,POSE.temp)

# clean data
POSE.bmass2$Temp_Ratio<-POSE.bmass2$Day_Scale_7/POSE.bmass2$Night_Scale_7
POSE.bmass2$Treatment<-as.factor(paste0(POSE.bmass2$Treatment))
POSE.bmass2$Day_Temp_Scale<-scale(POSE.bmass2$Day_Scale_7)
POSE.bmass2$Night_Temp_Scale<-scale(POSE.bmass2$Night_Scale_7)
#POSE.bmass2$Temp_Ratio<-POSE.bmass2$Day_Temp_Scale/POSE.bmass2$Night_Temp_Scale
#prelim figures
ggplot(POSE.bmass2, aes(y=AVG_RT_BMSS_g,x=Day_Temp_Scale,col=Treatment))+
  geom_point()+geom_smooth()
ggplot(POSE.bmass2, aes(y=AVG_SH_BMSS_g,x=Day_Temp_Scale,col=Treatment))+
  geom_point()+geom_smooth()
ggplot(POSE.bmass2, aes(y=Root_Biomass_g,x=Day_Temp_Scale,col=Treatment))+
  geom_point()+geom_smooth()
ggplot(POSE.bmass2, aes(y=Shoot_Biomass_g,x=Day_Temp_Scale,col=Treatment))+
  geom_point()+geom_smooth()

ggplot(POSE.bmass2, aes(y=AVG_RT_BMSS_g,x=Night_Temp_Scale,col=Treatment))+
  geom_point()+geom_smooth()
ggplot(POSE.bmass2, aes(y=AVG_SH_BMSS_g,x=Night_Temp_Scale,col=Treatment))+
  geom_point()+geom_smooth()
ggplot(POSE.bmass2, aes(y=Root_Biomass_g,x=Night_Temp_Scale,col=Treatment))+
  geom_point()+geom_smooth()
ggplot(POSE.bmass2, aes(y=Shoot_Biomass_g,x=Night_Temp_Scale,col=Treatment))+
  geom_point()+geom_smooth()

summary(POSE.bmass2$Root_Biomass_g)

ggplot(POSE.bmass2,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
  geom_tile(aes(fill=Root_Biomass_g))+
  facet_grid(.~Treatment)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.0006, space = "Lab",
                       name="Total Root Biomass") +
  theme_minimal()+
  theme(legend.position = "bottom")

summary(POSE.bmass2$Shoot_Biomass_g)

ggplot(POSE.bmass2,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
  geom_tile(aes(fill=Shoot_Biomass_g))+
  facet_grid(.~Treatment)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.001, space = "Lab",
                       name="Total Shoot Biomass") +
  theme_minimal()+
  theme(legend.position = "bottom")

summary(POSE.bmass2$AVG_RT_BMSS_g)

ggplot(POSE.bmass2,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
  geom_tile(aes(fill=AVG_RT_BMSS_g))+
  facet_grid(.~Treatment)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 1.8e-4, space = "Lab",
                       name="Individual Root Biomass") +
  theme_minimal()+
  theme(legend.position = "bottom")

summary(POSE.bmass2$AVG_SH_BMSS_g)

ggplot(POSE.bmass2,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
  geom_tile(aes(fill=AVG_SH_BMSS_g))+
  facet_grid(.~Treatment)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.0003, space = "Lab",
                       name="Individual Shoot Biomass") +
  theme_minimal()+
  theme(legend.position = "bottom")


#### shoot individual biomass glmmTMB ####
#rootratefit<-glmmTMB(Avg_Root_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment,
#                     family=Gamma(link = "log"),data=POSE_agg)

shoot.bmass.fit.full<-glmmTMB(AVG_SH_BMSS_g~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                            Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                            Day_Temp_Scale:Night_Temp_Scale+
                            Day_Temp_Scale:Night_Temp_Scale:Treatment,
                          family=Gamma(link="log"),data=POSE.bmass2)
summary(shoot.bmass.fit.full)


#### prediction ####
pnew.data<-POSE.bmass2[,c(1,5,19,24,25)]

pnew.data$Prediction<-predict(shoot.bmass.fit.full,pnew.data,allow.new.levels=T)
pnew.data$Prediction_back<-exp(pnew.data$Prediction)
ggplot(pnew.data,aes(x=POSE.bmass2$AVG_SH_BMSS_g,y=Prediction_back,col=Treatment))+
  geom_point()+geom_abline(slope=1)+theme_minimal()+geom_smooth(method="lm")
# plotly_day<-matrix(pnew.data$Day_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
# plotly_night<-matrix(pnew.data$Night_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
# plotly_zH<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="H")],byrow = T,ncol=7)
# plotly_zL<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="L")],byrow = T,ncol=7)
# plotly_zW<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="W")],byrow = T,ncol=7)
# plotly_zC<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="C")],byrow = T,ncol=7)
# fig<-plot_ly(showscale=F,type="surface", x=~plotly_day,y=~plotly_night)
# fig <- fig %>% add_surface(z = ~plotly_zC, opacity = 0.5, colorscale = list(c(0,1),c("rgb(107,255,184)","rgb(0,124,90)"))) #green
# fig <- fig %>% add_surface(z = ~plotly_zW, opacity = 0.5, colorscale = list(c(0,1),c("rgb(129,212,247)","rgb(12,177,247)"))) #blue
# fig <- fig %>% add_surface(z = ~plotly_zL, opacity = 0.5, colorscale = list(c(0,1),c("rgb(247,226,157)","rgb(255,192,7)"))) #yellow
# fig <- fig %>% add_surface(z = ~plotly_zH, opacity = 0.5, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) #pink
# fig

summary(pnew.data$Prediction_back)

ggplot(pnew.data,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
  geom_tile(aes(fill=Prediction_back))+
  facet_grid(.~Treatment)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 3.11e-4, space = "Lab",
                       name="Indiv Root Biomass") +
  theme_minimal()+
  theme(legend.position = "bottom")


#### root individual biomass glmmTMB ####
#rootratefit<-glmmTMB(Avg_Root_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment,
#                     family=Gamma(link = "log"),data=POSE_agg)

root.bmass.fit.full<-glmmTMB(AVG_RT_BMSS_g~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                                Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                                Day_Temp_Scale:Night_Temp_Scale+
                                Day_Temp_Scale:Night_Temp_Scale:Treatment,
                              family=Gamma(link="log"),data=POSE.bmass2)
summary(root.bmass.fit.full)


#### prediction ####
pnew.data<-POSE.bmass2[,c(1,5,18,24,25)]

pnew.data$Prediction<-predict(root.bmass.fit.full,pnew.data,allow.new.levels=T)
pnew.data$Prediction_back<-exp(pnew.data$Prediction)
ggplot(pnew.data,aes(x=POSE.bmass2$AVG_RT_BMSS_g,y=Prediction_back,col=Treatment))+
  geom_point()+geom_abline(slope=1)+theme_minimal()+geom_smooth(method="lm")
# plotly_day<-matrix(pnew.data$Day_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
# plotly_night<-matrix(pnew.data$Night_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
# plotly_zH<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="H")],byrow = T,ncol=7)
# plotly_zL<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="L")],byrow = T,ncol=7)
# plotly_zW<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="W")],byrow = T,ncol=7)
# plotly_zC<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="C")],byrow = T,ncol=7)
# fig<-plot_ly(showscale=F,type="surface", x=~plotly_day,y=~plotly_night)
# fig <- fig %>% add_surface(z = ~plotly_zC, opacity = 0.5, colorscale = list(c(0,1),c("rgb(107,255,184)","rgb(0,124,90)"))) #green
# fig <- fig %>% add_surface(z = ~plotly_zW, opacity = 0.5, colorscale = list(c(0,1),c("rgb(129,212,247)","rgb(12,177,247)"))) #blue
# fig <- fig %>% add_surface(z = ~plotly_zL, opacity = 0.5, colorscale = list(c(0,1),c("rgb(247,226,157)","rgb(255,192,7)"))) #yellow
# fig <- fig %>% add_surface(z = ~plotly_zH, opacity = 0.5, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) #pink
# fig
summary(pnew.data$Prediction_back)

ggplot(pnew.data,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
  geom_tile(aes(fill=Prediction_back))+
  facet_grid(.~Treatment)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 1.7e-4, space = "Lab",
                       name="Indiv Root Biomass") +
  theme_minimal()+
  theme(legend.position = "bottom")




#### shoot total biomass glmmTMB ####
#rootratefit<-glmmTMB(Avg_Root_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment,
#                     family=Gamma(link = "log"),data=POSE_agg)

shoot.bmass.fit.full<-glmmTMB(Shoot_Biomass_g~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                                Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                                Day_Temp_Scale:Night_Temp_Scale+
                                Day_Temp_Scale:Night_Temp_Scale:Treatment,
                              family=Gamma(link="log"),data=POSE.bmass2)
summary(shoot.bmass.fit.full)


#### prediction ####
pnew.data<-POSE.bmass2[,c(1,5,15,24,25)]

pnew.data$Prediction<-predict(shoot.bmass.fit.full,pnew.data,allow.new.levels=T)
pnew.data$Prediction_back<-exp(pnew.data$Prediction)
ggplot(pnew.data,aes(x=POSE.bmass2$Shoot_Biomass_g,y=Prediction_back,col=Treatment))+
  geom_point()+geom_abline(slope=1)+theme_minimal()
plotly_day<-matrix(pnew.data$Day_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
plotly_night<-matrix(pnew.data$Night_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
plotly_zH<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="H")],byrow = T,ncol=7)
plotly_zL<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="L")],byrow = T,ncol=7)
plotly_zW<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="W")],byrow = T,ncol=7)
plotly_zC<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="C")],byrow = T,ncol=7)
fig<-plot_ly(showscale=F,type="surface", x=~plotly_day,y=~plotly_night)
fig <- fig %>% add_surface(z = ~plotly_zC, opacity = 0.5, colorscale = list(c(0,1),c("rgb(107,255,184)","rgb(0,124,90)"))) #green
fig <- fig %>% add_surface(z = ~plotly_zW, opacity = 0.5, colorscale = list(c(0,1),c("rgb(129,212,247)","rgb(12,177,247)"))) #blue
fig <- fig %>% add_surface(z = ~plotly_zL, opacity = 0.5, colorscale = list(c(0,1),c("rgb(247,226,157)","rgb(255,192,7)"))) #yellow
fig <- fig %>% add_surface(z = ~plotly_zH, opacity = 0.5, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) #pink
fig

ggplot(pnew.data,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
  geom_tile(aes(fill=Prediction_back))+
  facet_grid(.~Treatment)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.001, space = "Lab",
                       name="Total Shoot Biomass") +
  theme_minimal()+
  theme(legend.position = "bottom")


#### root total biomass glmmTMB ####
#rootratefit<-glmmTMB(Avg_Root_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment,
#                     family=Gamma(link = "log"),data=POSE_agg)

root.bmass.fit.full<-glmmTMB(Root_Biomass_g~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                               Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                               Day_Temp_Scale:Night_Temp_Scale+
                               Day_Temp_Scale:Night_Temp_Scale:Treatment,
                             family=Gamma(link="log"),data=POSE.bmass2)
summary(root.bmass.fit.full)


#### prediction ####
pnew.data<-POSE.bmass2[,c(1,5,14,24,25)]

pnew.data$Prediction<-predict(root.bmass.fit.full,pnew.data,allow.new.levels=T)
pnew.data$Prediction_back<-exp(pnew.data$Prediction)
ggplot(pnew.data,aes(x=POSE.bmass2$Root_Biomass_g,y=Prediction_back,col=Treatment))+
  geom_point()+geom_abline(slope=1)+theme_minimal()
# plotly_day<-matrix(pnew.data$Day_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
# plotly_night<-matrix(pnew.data$Night_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
# plotly_zH<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="H")],byrow = T,ncol=7)
# plotly_zL<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="L")],byrow = T,ncol=7)
# plotly_zW<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="W")],byrow = T,ncol=7)
# plotly_zC<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="C")],byrow = T,ncol=7)
# fig<-plot_ly(showscale=F,type="surface", x=~plotly_day,y=~plotly_night)
# fig <- fig %>% add_surface(z = ~plotly_zC, opacity = 0.5, colorscale = list(c(0,1),c("rgb(107,255,184)","rgb(0,124,90)"))) #green
# fig <- fig %>% add_surface(z = ~plotly_zW, opacity = 0.5, colorscale = list(c(0,1),c("rgb(129,212,247)","rgb(12,177,247)"))) #blue
# fig <- fig %>% add_surface(z = ~plotly_zL, opacity = 0.5, colorscale = list(c(0,1),c("rgb(247,226,157)","rgb(255,192,7)"))) #yellow
# fig <- fig %>% add_surface(z = ~plotly_zH, opacity = 0.5, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) #pink
# fig

summary(pnew.data$Prediction_back)

ggplot(pnew.data,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
  geom_tile(aes(fill=Prediction_back))+
  facet_grid(.~Treatment)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.0005, space = "Lab",
                       name="Total Root Biomass") +
  theme_minimal()+
  theme(legend.position = "bottom")

#### total individual biomass glmmTMB ####
#rootratefit<-glmmTMB(Avg_Root_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment,
#                     family=Gamma(link = "log"),data=POSE_agg)

POSE.bmass2$Indiv_Bmass_g<-POSE.bmass2$AVG_RT_BMSS_g+POSE.bmass2$AVG_SH_BMSS_g

indiv.bmass.fit.full<-glmmTMB(Indiv_Bmass_g~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                                Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                                Day_Temp_Scale:Night_Temp_Scale+
                                Day_Temp_Scale:Night_Temp_Scale:Treatment,
                              family=Gamma(link="log"),data=POSE.bmass2)
summary(indiv.bmass.fit.full)


#### prediction ####
pnew.data<-POSE.bmass2[,c(1,5,26,24,25)]

pnew.data$Prediction<-predict(indiv.bmass.fit.full,pnew.data,allow.new.levels=T)
pnew.data$Prediction_back<-exp(pnew.data$Prediction)
ggplot(pnew.data,aes(x=POSE.bmass2$Indiv_Bmass_g,y=Prediction_back,col=Treatment))+
  geom_point()+geom_abline(slope=1)+theme_minimal()
# plotly_day<-matrix(pnew.data$Day_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
# plotly_night<-matrix(pnew.data$Night_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
# plotly_zH<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="H")],byrow = T,ncol=7)
# plotly_zL<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="L")],byrow = T,ncol=7)
# plotly_zW<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="W")],byrow = T,ncol=7)
# plotly_zC<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="C")],byrow = T,ncol=7)
# fig<-plot_ly(showscale=F,type="surface", x=~plotly_day,y=~plotly_night)
# fig <- fig %>% add_surface(z = ~plotly_zC, opacity = 0.5, colorscale = list(c(0,1),c("rgb(107,255,184)","rgb(0,124,90)"))) #green
# fig <- fig %>% add_surface(z = ~plotly_zW, opacity = 0.5, colorscale = list(c(0,1),c("rgb(129,212,247)","rgb(12,177,247)"))) #blue
# fig <- fig %>% add_surface(z = ~plotly_zL, opacity = 0.5, colorscale = list(c(0,1),c("rgb(247,226,157)","rgb(255,192,7)"))) #yellow
# fig <- fig %>% add_surface(z = ~plotly_zH, opacity = 0.5, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) #pink
# fig

summary(pnew.data$Prediction_back)

ggplot(pnew.data,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
  geom_tile(aes(fill=Prediction_back))+
  facet_grid(.~Treatment)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.0005, space = "Lab",
                       name="Total Individual Biomass") +
  theme_minimal()+
  theme(legend.position = "bottom")



#### total cuvette biomass glmmTMB ####
#rootratefit<-glmmTMB(Avg_Root_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment,
#                     family=Gamma(link = "log"),data=POSE_agg)
POSE.bmass2$Total_Bmass_g<-POSE.bmass2$Shoot_Biomass_g+POSE.bmass2$Root_Biomass_g



total.bmass.fit.full<-glmmTMB(Total_Bmass_g~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                               Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                               Day_Temp_Scale:Night_Temp_Scale+
                               Day_Temp_Scale:Night_Temp_Scale:Treatment,
                             family=Gamma(link="log"),data=POSE.bmass2)
summary(total.bmass.fit.full)


#### prediction ####
pnew.data<-POSE.bmass2[,c(1,5,27,24,25)]

pnew.data$Prediction<-predict(total.bmass.fit.full,pnew.data,allow.new.levels=T)
pnew.data$Prediction_back<-exp(pnew.data$Prediction)
ggplot(pnew.data,aes(x=POSE.bmass2$Total_Bmass_g,y=Prediction_back,col=Treatment))+
  geom_point()+geom_abline(slope=1)+theme_minimal()
# plotly_day<-matrix(pnew.data$Day_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
# plotly_night<-matrix(pnew.data$Night_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
# plotly_zH<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="H")],byrow = T,ncol=7)
# plotly_zL<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="L")],byrow = T,ncol=7)
# plotly_zW<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="W")],byrow = T,ncol=7)
# plotly_zC<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="C")],byrow = T,ncol=7)
# fig<-plot_ly(showscale=F,type="surface", x=~plotly_day,y=~plotly_night)
# fig <- fig %>% add_surface(z = ~plotly_zC, opacity = 0.5, colorscale = list(c(0,1),c("rgb(107,255,184)","rgb(0,124,90)"))) #green
# fig <- fig %>% add_surface(z = ~plotly_zW, opacity = 0.5, colorscale = list(c(0,1),c("rgb(129,212,247)","rgb(12,177,247)"))) #blue
# fig <- fig %>% add_surface(z = ~plotly_zL, opacity = 0.5, colorscale = list(c(0,1),c("rgb(247,226,157)","rgb(255,192,7)"))) #yellow
# fig <- fig %>% add_surface(z = ~plotly_zH, opacity = 0.5, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) #pink
# fig

summary(pnew.data$Prediction_back)

ggplot(pnew.data,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
  geom_tile(aes(fill=Prediction_back))+
  facet_grid(.~Treatment)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.0016, space = "Lab",
                       name="Total Cuvette Biomass") +
  theme_minimal()+
  theme(legend.position = "bottom")









