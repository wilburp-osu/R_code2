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
ARTR<-read.csv("Data/ARTR_Data/ARTR_Data_FINAL.csv")
days<-unique(ARTR$Day_Collected)


### clean data ####
ARTR$Root_Pheno_Code<-as.numeric(paste0(ARTR$Root_Pheno_Code))
ARTR$Shoot_Pheno_Code<-as.numeric(paste0(ARTR$Shoot_Pheno_Code))

ARTR$Root_Pheno_Code[which(ARTR$Root_Pheno_Code==999)]<-NA
ARTR$Shoot_Pheno_Code[which(ARTR$Shoot_Pheno_Code==999)]<-NA

ARTR$Avg_Root_rate<-NA
ARTR$Avg_Shoot_rate<-NA

Root_Code<-c(1,2,3,4)
Root_Length_mid<-c(0.5,5,10,20)
Root_match<-data.frame(Code=Root_Code,Length=Root_Length_mid)
ARTR$Root_Length<-Root_match$Length[match(ARTR$Root_Pheno_Code,Root_match$Code)]

Shoot_Code<-c(6,7,8)
Shoot_Length_mid<-c(0.5,10,20)
Shoot_match<-data.frame(Code=Shoot_Code,Length=Shoot_Length_mid)
ARTR$Shoot_Length<-Shoot_match$Length[match(ARTR$Shoot_Pheno_Code,Shoot_match$Code)]


for (i in 1:max(ARTR$Group)){
  for (t in unique(ARTR$Treatment)) {
    for (c in 1:4) {
      ARTR_seed<-ARTR[which(ARTR$Group==i & ARTR$Treatment==t & ARTR$Cell==c),]
      if (any(ARTR_seed$Root_Pheno_Code==1)){
        root_days<-ARTR_seed$Root_Length[
          which(ARTR_seed$Root_Pheno_Code>0&ARTR_seed$Root_Pheno_Code<5)]
        #root_days_dup<-unique(root_days)
        #root_rate<-sum(root_days_dup)/length(root_days)
        root_rate<-max(root_days)/length(root_days)
      } else {root_rate<-0}
      if (any(ARTR_seed$Shoot_Pheno_Code==6)){
        shoot_days<-(ARTR_seed$Shoot_Length[
          which(ARTR_seed$Shoot_Pheno_Code>5&ARTR_seed$Shoot_Pheno_Code<9)])
        #shoot_days_dup<-unique(shoot_days)
        #shoot_rate<-sum(shoot_days_dup)/length(shoot_days)
        shoot_rate<-max(shoot_days)/length(shoot_days)
      } else {shoot_rate<-0}
      ARTR$Avg_Root_rate[which(ARTR$Group==i & ARTR$Treatment==t & ARTR$Cell==c)]<-root_rate
      ARTR$Avg_Shoot_rate[which(ARTR$Group==i & ARTR$Treatment==t & ARTR$Cell==c)]<-shoot_rate
    }
  }
}



### aggregate data for each cuvette for each day
ARTR_agg<-aggregate(cbind(ARTR$Avg_Root_rate,
                          ARTR$Avg_Shoot_rate),
                    by=list(ARTR$Day_Scale_7,ARTR$Night_Scale_7,
                            ARTR$Treatment,
                            ARTR$Group),mean,na.rm=T)

ARTR_agg[ARTR_agg =="NaN"]<-NA
names(ARTR_agg)<-c("Day_Temp","Night_Temp","Treatment",
                   "Group","Avg_Root_rate","Avg_Shoot_rate")


ARTR_agg$Temp_Ratio<-ARTR_agg$Day_Temp/ARTR_agg$Night_Temp
ARTR_agg$Treatment<-as.factor(paste0(ARTR_agg$Treatment))
ARTR_agg$Day_Temp_Scale<-scale(ARTR_agg$Day_Temp)
ARTR_agg$Night_Temp_Scale<-scale(ARTR_agg$Night_Temp)

#### root rate glmmTMB ####
#rootratefit<-glmmTMB(Avg_Root_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment,
#                     family=Gamma(link = "log"),data=ARTR_agg)

rootratefit.full<-glmmTMB(Avg_Root_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                            Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                            Day_Temp_Scale:Night_Temp_Scale+
                            Day_Temp_Scale:Night_Temp_Scale:Treatment,
                          family=Gamma(link="log"),data=ARTR_agg)
summary(rootratefit.full)

rootratefit<-glmmTMB(Avg_Root_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                       Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                       Day_Temp_Scale:Night_Temp_Scale,
                     family=Gamma(link="log"),data=ARTR_agg)
summary(rootratefit)



#### prediction ####
pnew.data<-ARTR_agg[,c(1,2,3,4,8,9)]

pnew.data$Prediction<-predict(rootratefit,pnew.data,allow.new.levels=T)
pnew.data$Prediction_back<-exp(pnew.data$Prediction)
ggplot(pnew.data,aes(x=ARTR_agg$Avg_Root_rate,y=Prediction_back,col=Treatment))+
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
                           family=ziGamma(link="log"),data=ARTR_agg,ziformula=~1)
summary(shootratefit.full)

shootratefit<-glmmTMB(Avg_Shoot_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                        Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                        Day_Temp_Scale:Night_Temp_Scale,
                      family=ziGamma(link="log"),data=ARTR_agg,ziformula=~1)
summary(shootratefit)

#### prediction ####
pnew.data.sh<-ARTR_agg[,c(1,2,3,4,8,9)]

pnew.data.sh$Prediction<-predict(shootratefit.full,pnew.data.sh)
pnew.data.sh$Prediction_back<-exp(pnew.data.sh$Prediction)
ggplot(pnew.data.sh,aes(x=ARTR_agg$Avg_Shoot_rate,y=Prediction_back,col=Treatment))+
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
ARTR_root_H<-matrix(ARTR_agg$Avg_Root_rate[which(ARTR_agg$Treatment=="H")],ncol=7,byrow=T)
ARTR_root_HR<-raster(ARTR_root_H,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(ARTR_root_HR,main="Root Growth Rate - High",
     breaks=round(seq(0,9,length.out=10),digits=2),col=terrain.colors(10,rev=T))
ARTR_root_W<-matrix(ARTR_agg$Avg_Root_rate[which(ARTR_agg$Treatment=="W")],ncol=7,byrow=T)
ARTR_root_WR<-raster(ARTR_root_W,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(ARTR_root_WR,main="Root Growth Rate - Water",
     breaks=round(seq(0,9,length.out=10),digits=2),col=terrain.colors(10,rev=T))
ARTR_root_L<-matrix(ARTR_agg$Avg_Root_rate[which(ARTR_agg$Treatment=="L")],ncol=7,byrow=T)
ARTR_root_LR<-raster(ARTR_root_L,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(ARTR_root_LR,main="Root Growth Rate - Low",
     breaks=round(seq(0,9,length.out=10),digits=2),col=terrain.colors(10,rev=T))
ARTR_root_C<-matrix(ARTR_agg$Avg_Root_rate[which(ARTR_agg$Treatment=="C")],ncol=7,byrow=T)
ARTR_root_CR<-raster(ARTR_root_C,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(ARTR_root_CR,main="Root Growth Rate - Control",
     breaks=round(seq(0,9,length.out=10),digits=2),col=terrain.colors(10,rev=T))

#### plot raster shoot growth rate ####
ARTR_shoot_H<-matrix(ARTR_agg$Avg_Shoot_rate[which(ARTR_agg$Treatment=="H")],ncol=7,byrow=T)
ARTR_shoot_HR<-raster(ARTR_shoot_H,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(ARTR_shoot_HR,main="Shoot Growth Rate - High",
     breaks=round(seq(0,12,length.out=10),digits=2),col=terrain.colors(10,rev=T))
ARTR_shoot_W<-matrix(ARTR_agg$Avg_Shoot_rate[which(ARTR_agg$Treatment=="W")],ncol=7,byrow=T)
ARTR_shoot_WR<-raster(ARTR_shoot_W,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(ARTR_shoot_WR,main="Shoot Growth Rate - Water",
     breaks=round(seq(0,12,length.out=10),digits=2),col=terrain.colors(10,rev=T))
ARTR_shoot_L<-matrix(ARTR_agg$Avg_Shoot_rate[which(ARTR_agg$Treatment=="L")],ncol=7,byrow=T)
ARTR_shoot_LR<-raster(ARTR_shoot_L,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(ARTR_shoot_LR,main="Shoot Growth Rate - Low",
     breaks=round(seq(0,12,length.out=10),digits=2),col=terrain.colors(10,rev=T))
ARTR_shoot_C<-matrix(ARTR_agg$Avg_Shoot_rate[which(ARTR_agg$Treatment=="C")],ncol=7,byrow=T)
ARTR_shoot_CR<-raster(ARTR_shoot_C,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(ARTR_shoot_CR,main="Shoot Growth Rate - Control",
     breaks=round(seq(0,12,length.out=10),digits=2),col=terrain.colors(10,rev=T))

