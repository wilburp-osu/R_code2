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
POSE<-read.csv("Data/POSE_Data/SA_POSE_Data_TEM.csv")
days<-unique(POSE$Day_Collected)


### clean data ####
POSE$Root_Pheno_Code<-as.numeric(paste0(POSE$Root_Pheno_Code))
POSE$Shoot_Pheno_Code<-as.numeric(paste0(POSE$Shoot_Pheno_Code))

POSE$Root_Pheno_Code[which(POSE$Root_Pheno_Code==999)]<-NA
POSE$Shoot_Pheno_Code[which(POSE$Shoot_Pheno_Code==999)]<-NA

POSE$Avg_Root_rate<-NA
POSE$Avg_Shoot_rate<-NA

Root_Code<-c(1,2,3,4)
Root_Length_mid<-c(0.5,5,10,20)
Root_match<-data.frame(Code=Root_Code,Length=Root_Length_mid)
POSE$Root_Length<-Root_match$Length[match(POSE$Root_Pheno_Code,Root_match$Code)]

Shoot_Code<-c(6,7,8)
Shoot_Length_mid<-c(0.5,10,20)
Shoot_match<-data.frame(Code=Shoot_Code,Length=Shoot_Length_mid)
POSE$Shoot_Length<-Shoot_match$Length[match(POSE$Shoot_Pheno_Code,Shoot_match$Code)]




for (i in 1:max(POSE$Group)){
  for (t in unique(POSE$Treatment)) {
    for (c in 1:4) {
      POSE_seed<-POSE[which(POSE$Group==i & POSE$Treatment==t & POSE$Cell==c),]
      if (any(POSE_seed$Root_Pheno_Code==1)){
      root_days<-POSE_seed$Root_Length[which(POSE_seed$Root_Pheno_Code>0&POSE_seed$Root_Pheno_Code<5)]
      #root_days_dup<-unique(root_days)
      #root_rate<-sum(root_days_dup)/length(root_days)
      root_rate<-max(root_days)/length(root_days)
      } else {root_rate<-0}
      if (any(POSE_seed$Shoot_Pheno_Code==6)){
      shoot_days<-(POSE_seed$Shoot_Length[which(POSE_seed$Shoot_Pheno_Code>5&POSE_seed$Shoot_Pheno_Code<9)])
      #shoot_days_dup<-unique(shoot_days)
      #shoot_rate<-sum(shoot_days_dup)/length(shoot_days)
      shoot_rate<-max(shoot_days)/length(shoot_days)
      } else {shoot_rate<-0}
      POSE$Avg_Root_rate[which(POSE$Group==i & POSE$Treatment==t & POSE$Cell==c)]<-root_rate
      POSE$Avg_Shoot_rate[which(POSE$Group==i & POSE$Treatment==t & POSE$Cell==c)]<-shoot_rate
    }
  }
}





### aggregate data for each cuvette for each day
POSE_agg<-aggregate(cbind(POSE$Avg_Root_rate,
                          POSE$Avg_Shoot_rate),
                    by=list(POSE$Day_Scale_7,POSE$Night_Scale_7,
                            POSE$Treatment,
                            POSE$Group),mean,na.rm=T)

POSE_agg[POSE_agg =="NaN"]<-NA
names(POSE_agg)<-c("Day_Temp","Night_Temp","Treatment",
                   "Group","Avg_Root_rate","Avg_Shoot_rate")


POSE_agg$Temp_Ratio<-POSE_agg$Day_Temp/POSE_agg$Night_Temp
POSE_agg$Treatment<-as.factor(paste0(POSE_agg$Treatment))
POSE_agg$Day_Temp_Scale<-scale(POSE_agg$Day_Temp)
POSE_agg$Night_Temp_Scale<-scale(POSE_agg$Night_Temp)

#### root rate glmmTMB ####
#rootratefit<-glmmTMB(Avg_Root_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment,
#                     family=Gamma(link = "log"),data=POSE_agg)

rootratefit.full<-glmmTMB(Avg_Root_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                       Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                       Day_Temp_Scale:Night_Temp_Scale+
                       Day_Temp_Scale:Night_Temp_Scale:Treatment,
                     family=Gamma(link="log"),data=POSE_agg)
summary(rootratefit.full)

rootratefit<-glmmTMB(Avg_Root_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                       Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                       Day_Temp_Scale:Night_Temp_Scale,
                     family=Gamma(link="log"),data=POSE_agg)
summary(rootratefit)



#### prediction ####
pnew.data<-POSE_agg[,c(1,2,3,4,8,9)]

pnew.data$Prediction<-predict(rootratefit,pnew.data,allow.new.levels=T)
pnew.data$Prediction_back<-exp(pnew.data$Prediction)
ggplot(pnew.data,aes(x=POSE_agg$Avg_Root_rate,y=Prediction_back,col=Treatment))+
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

summary(pnew.data$Prediction_back)

ggplot(pnew.data,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
  geom_raster(aes(fill=Prediction_back),interpolate=T)+
  facet_grid(.~Treatment)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 2.4, space = "Lab",
                       name="Root growth rate (final root length/days of growth)") +
  theme_minimal()+
  theme(legend.position = "bottom")


#### shoot rate glmmTMB ####
shootratefit.full<-glmmTMB(Avg_Shoot_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                        Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                        Day_Temp_Scale:Night_Temp_Scale+
                        Day_Temp_Scale:Night_Temp_Scale:Treatment,
                    family=ziGamma(link="log"),data=POSE_agg,ziformula=~1)
summary(shootratefit.full)

shootratefit<-glmmTMB(Avg_Shoot_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                             Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                             Day_Temp_Scale:Night_Temp_Scale,
                           family=ziGamma(link="log"),data=POSE_agg,ziformula=~1)
summary(shootratefit)

#### prediction ####
pnew.data.sh<-POSE_agg[,c(1,2,3,4,8,9)]

pnew.data.sh$Prediction<-predict(shootratefit.full,pnew.data.sh)
pnew.data.sh$Prediction_back<-exp(pnew.data.sh$Prediction)
ggplot(pnew.data.sh,aes(x=POSE_agg$Avg_Shoot_rate,y=Prediction_back,col=Treatment))+
  geom_point()+geom_abline(slope=1)+theme_minimal()
plotly_day<-matrix(pnew.data$Day_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
plotly_night<-matrix(pnew.data$Night_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
plotly_zH<-matrix(pnew.data.sh$Prediction_back[which(pnew.data.sh$Treatment=="H")],byrow = T,ncol=7)
plotly_zL<-matrix(pnew.data.sh$Prediction_back[which(pnew.data.sh$Treatment=="L")],byrow = T,ncol=7)
plotly_zW<-matrix(pnew.data.sh$Prediction_back[which(pnew.data.sh$Treatment=="W")],byrow = T,ncol=7)
plotly_zC<-matrix(pnew.data.sh$Prediction_back[which(pnew.data.sh$Treatment=="C")],byrow = T,ncol=7)
fig<-plot_ly(showscale=F,type="surface", x=~plotly_day,y=~plotly_night)
fig <- fig %>% add_surface(z = ~plotly_zC, opacity = 0.5, colorscale = list(c(0,1),c("rgb(107,255,184)","rgb(0,124,90)"))) #green
fig <- fig %>% add_surface(z = ~plotly_zW, opacity = 0.5, colorscale = list(c(0,1),c("rgb(129,212,247)","rgb(12,177,247)"))) #blue
fig <- fig %>% add_surface(z = ~plotly_zL, opacity = 0.5, colorscale = list(c(0,1),c("rgb(247,226,157)","rgb(255,192,7)"))) #yellow
fig <- fig %>% add_surface(z = ~plotly_zH, opacity = 0.5, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) #pink
fig

summary(pnew.data$Prediction_back)

ggplot(pnew.data,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
  geom_raster(aes(fill=Prediction_back),interpolate=T)+
  facet_grid(.~Treatment)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 2.4, space = "Lab",
                       name="Shoot growth rate (final root length/days of growth)") +
  theme_minimal()+
  theme(legend.position = "bottom")


#### plot raster root growth rate ####
POSE_root_H<-matrix(POSE_agg$Avg_Root_rate[which(POSE_agg$Treatment=="H")],ncol=7,byrow=T)
POSE_root_HR<-raster(POSE_root_H,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_root_HR,main="Root Growth Rate - High",
     breaks=round(seq(0,9,length.out=10),digits=2),col=terrain.colors(10,rev=T))
POSE_root_W<-matrix(POSE_agg$Avg_Root_rate[which(POSE_agg$Treatment=="W")],ncol=7,byrow=T)
POSE_root_WR<-raster(POSE_root_W,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_root_WR,main="Root Growth Rate - Water",
     breaks=round(seq(0,9,length.out=10),digits=2),col=terrain.colors(10,rev=T))
POSE_root_L<-matrix(POSE_agg$Avg_Root_rate[which(POSE_agg$Treatment=="L")],ncol=7,byrow=T)
POSE_root_LR<-raster(POSE_root_L,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_root_LR,main="Root Growth Rate - Low",
     breaks=round(seq(0,9,length.out=10),digits=2),col=terrain.colors(10,rev=T))
POSE_root_C<-matrix(POSE_agg$Avg_Root_rate[which(POSE_agg$Treatment=="C")],ncol=7,byrow=T)
POSE_root_CR<-raster(POSE_root_C,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_root_CR,main="Root Growth Rate - Control",
     breaks=round(seq(0,9,length.out=10),digits=2),col=terrain.colors(10,rev=T))

#### plot raster shoot growth rate ####
POSE_shoot_H<-matrix(POSE_agg$Avg_Shoot_rate[which(POSE_agg$Treatment=="H")],ncol=7,byrow=T)
POSE_shoot_HR<-raster(POSE_shoot_H,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_shoot_HR,main="Shoot Growth Rate - High",
     breaks=round(seq(0,12,length.out=10),digits=2),col=terrain.colors(10,rev=T))
POSE_shoot_W<-matrix(POSE_agg$Avg_Shoot_rate[which(POSE_agg$Treatment=="W")],ncol=7,byrow=T)
POSE_shoot_WR<-raster(POSE_shoot_W,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_shoot_WR,main="Shoot Growth Rate - Water",
     breaks=round(seq(0,12,length.out=10),digits=2),col=terrain.colors(10,rev=T))
POSE_shoot_L<-matrix(POSE_agg$Avg_Shoot_rate[which(POSE_agg$Treatment=="L")],ncol=7,byrow=T)
POSE_shoot_LR<-raster(POSE_shoot_L,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_shoot_LR,main="Shoot Growth Rate - Low",
     breaks=round(seq(0,12,length.out=10),digits=2),col=terrain.colors(10,rev=T))
POSE_shoot_C<-matrix(POSE_agg$Avg_Shoot_rate[which(POSE_agg$Treatment=="C")],ncol=7,byrow=T)
POSE_shoot_CR<-raster(POSE_shoot_C,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_shoot_CR,main="Shoot Growth Rate - Control",
     breaks=round(seq(0,12,length.out=10),digits=2),col=terrain.colors(10,rev=T))


