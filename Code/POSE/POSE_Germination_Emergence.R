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

### calculate days til emergence and germination
POSE$germ_rate<-NA
POSE$emerg_rate<-NA

for (i in 1:max(POSE$Group)){
  for (t in unique(POSE$Treatment)) {
    for (c in 1:4) {
      POSE_seed<-POSE[which(POSE$Group==i & POSE$Treatment==t & POSE$Cell==c),]
      if (any(POSE_seed$Root_Pheno_Code==1)){
      Days<-POSE_seed$Days_Since_Install[which(POSE_seed$Root_Pheno_Code==1)]
      POSE$germ_rate[which(POSE$Group==i & POSE$Treatment==t & POSE$Cell==c)]<-1/Days
      } else {
      POSE$germ_rate[which(POSE$Group==i & POSE$Treatment==t & POSE$Cell==c)]<-0
      }
      if (any(POSE_seed$Shoot_Pheno_Code==6)){
        Days2<-POSE_seed$Days_Since_Install[which(POSE_seed$Shoot_Pheno_Code==6)]
        POSE$emerg_rate[which(POSE$Group==i & POSE$Treatment==t & POSE$Cell==c)]<-1/Days2  
      } else {
        POSE$emerg_rate[which(POSE$Group==i & POSE$Treatment==t & POSE$Cell==c)]<-0
      }
    }
  }
}
        




### aggregate data for each cuvette for each day
POSE_agg<-aggregate(cbind(POSE$Shoot_Pheno_Code,
                POSE$Root_Pheno_Code,
                POSE$germ_rate,
                POSE$emerg_rate),
                by=list(POSE$Day_Scale_7,POSE$Night_Scale_7,
                        POSE$Treatment,
                        POSE$Days_Since_Install,
                        POSE$Group),mean,na.rm=T)

POSE_agg[POSE_agg =="NaN"]<-NA
names(POSE_agg)<-c("Day_Temp","Night_Temp","Treatment",
                   "Days_Since_Install","Group",
                   "Root_Pheno_Code","Shoot_Pheno_Code",
                   "germ_rate","emerg_rate")


#### models for germination and emergence ####
#aggregate across all days 
POSE_agg2<-aggregate(cbind(POSE_agg$germ_rate,
                         POSE_agg$emerg_rate),
                   by=list(POSE_agg$Day_Temp,POSE_agg$Night_Temp,
                           POSE_agg$Treatment,
                           POSE_agg$Group),mean,na.rm=T)

POSE_agg2[POSE_agg2 =="NaN"]<-NA
names(POSE_agg2)<-c("Day_Temp","Night_Temp","Treatment",
                   "Group","germ_rate","emerg_rate")

POSE_agg2$Day_Temp_Scale<-scale(POSE_agg2$Day_Temp)
POSE_agg2$Night_Temp_Scale<-scale(POSE_agg2$Night_Temp)

#### germination glmmTMB ####
# added 3=way interaction but was not significant so dropped it
rootDaysfit<-glmmTMB(germ_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                       Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+
                       Treatment+Day_Temp_Scale:Night_Temp_Scale,
                     family=Gamma(link="log"),data=POSE_agg2)
summary(rootDaysfit)
# interactions between treatment and temps are not significant so dropped them
rootDaysfit<-glmmTMB(germ_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                       Treatment+Day_Temp_Scale:Night_Temp_Scale,
                     family=Gamma(link="log"),data=POSE_agg2)
summary(rootDaysfit)

#### prediction ####
pnew.data<-POSE_agg2[,c(1,2,3,4,7,8)]

pnew.data$Prediction<-predict(rootDaysfit,pnew.data)
pnew.data$Prediction_back<-exp(pnew.data$Prediction)
ggplot(pnew.data,aes(x=POSE_agg2$germ_rate,y=Prediction_back,col=Treatment))+
  geom_point()+geom_abline(slope=1)+theme_minimal()

plotly_x<-matrix(pnew.data$Day_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
plotly_y<-matrix(pnew.data$Night_Temp[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
plotly_zH<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="H")],byrow = T,ncol=7)
plotly_zL<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="L")],byrow = T,ncol=7)
plotly_zW<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="W")],byrow = T,ncol=7)
plotly_zC<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="C")],byrow = T,ncol=7)
fig<-plot_ly(showscale=F,type="surface", x=~plotly_x,y=plotly_y)
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
                       midpoint = 0.17, space = "Lab",
                       name="Germination Rate (1/days til germ)") +
  theme_minimal()+
  theme(legend.position = "bottom")


####emergence glmmTMB####
shootDaysfit<-glmmTMB(emerg_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                        Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+
                        Treatment+Day_Temp_Scale:Night_Temp_Scale,
                      family=ziGamma(link="log"),data=POSE_agg2,ziformula=~1)
summary(shootDaysfit)
# interactions between treatment and temps are not significant so dropped them
shootDaysfit<-glmmTMB(emerg_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                        Treatment,
                      family=ziGamma(link="log"),data=POSE_agg2,ziformula=~1)
summary(shootDaysfit)


#### prediction ####
pnew.data.emerg<-POSE_agg2[,c(1,2,3,4,7,8)]

pnew.data.emerg$Prediction<-predict(shootDaysfit,pnew.data.emerg)
pnew.data.emerg$Prediction_back<-exp(pnew.data.emerg$Prediction)
ggplot(pnew.data.emerg,aes(x=POSE_agg2$emerg_rate,y=Prediction_back,col=Treatment))+
  geom_point()+geom_abline(slope=1)+theme_minimal()

plotly_x<-matrix(pnew.data.emerg$Day_Temp[which(pnew.data.emerg$Treatment=="H")],byrow=T,ncol=7)
plotly_y<-matrix(pnew.data.emerg$Night_Temp[which(pnew.data.emerg$Treatment=="H")],byrow=T,ncol=7)
plotly_zH<-matrix(pnew.data.emerg$Prediction_back[which(pnew.data.emerg$Treatment=="H")],byrow = T,ncol=7)
plotly_zL<-matrix(pnew.data.emerg$Prediction_back[which(pnew.data.emerg$Treatment=="L")],byrow = T,ncol=7)
plotly_zW<-matrix(pnew.data.emerg$Prediction_back[which(pnew.data.emerg$Treatment=="W")],byrow = T,ncol=7)
plotly_zC<-matrix(pnew.data.emerg$Prediction_back[which(pnew.data.emerg$Treatment=="C")],byrow = T,ncol=7)
fig.emerg<-plot_ly(showscale=F,type="surface", x=~plotly_x,y=plotly_y)
fig.emerg <- fig.emerg %>% add_surface(z = ~plotly_zC, opacity = 0.5, colorscale = list(c(0,1),c("rgb(107,255,184)","rgb(0,124,90)"))) #green
fig.emerg <- fig.emerg %>% add_surface(z = ~plotly_zW, opacity = 0.5, colorscale = list(c(0,1),c("rgb(129,212,247)","rgb(12,177,247)"))) #blue
fig.emerg <- fig.emerg %>% add_surface(z = ~plotly_zL, opacity = 0.5, colorscale = list(c(0,1),c("rgb(247,226,157)","rgb(255,192,7)"))) #yellow
fig.emerg <- fig.emerg %>% add_surface(z = ~plotly_zH, opacity = 0.5, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) #pink
fig.emerg

summary(pnew.data$Prediction_back)

ggplot(pnew.data,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
  geom_raster(aes(fill=Prediction_back),interpolate=T)+
  facet_grid(.~Treatment)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.17, space = "Lab",
                       name="Emergence Rate (1/days til emergence)") +
  theme_minimal()+
  theme(legend.position = "bottom")


##########################
### plot raster germination
POSE_germ_H<-matrix(POSE_agg2$germ_rate[which(POSE_agg2$Treatment=="H")],ncol=7,byrow=T)
POSE_germ_HR<-raster(POSE_germ_H,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_germ_HR,main="Germination Rate - High",
     breaks=round(seq(0,0.5,length.out=10),digits=2),col=terrain.colors(10,rev=T))
POSE_germ_W<-matrix(POSE_agg2$germ_rate[which(POSE_agg2$Treatment=="W")],ncol=7,byrow=T)
POSE_germ_WR<-raster(POSE_germ_W,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_germ_WR,main="Germination Rate - Water",
     breaks=round(seq(0,0.5,length.out=10),digits=2),col=terrain.colors(10,rev=T))
POSE_germ_L<-matrix(POSE_agg2$germ_rate[which(POSE_agg2$Treatment=="L")],ncol=7,byrow=T)
POSE_germ_LR<-raster(POSE_germ_L,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_germ_LR,main="Germination Rate - Low",
     breaks=round(seq(0,0.5,length.out=10),digits=2),col=terrain.colors(10,rev=T))
POSE_germ_C<-matrix(POSE_agg2$germ_rate[which(POSE_agg2$Treatment=="C")],ncol=7,byrow=T)
POSE_germ_CR<-raster(POSE_germ_C,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_germ_CR,main="Germination Rate - Control",
     breaks=round(seq(0,0.5,length.out=10),digits=2),col=terrain.colors(10,rev=T))

### plot raster emergence
POSE_emerg_H<-matrix(POSE_agg2$emerg_rate[which(POSE_agg2$Treatment=="H")],ncol=7,byrow=T)
POSE_emerg_HR<-raster(POSE_emerg_H,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_emerg_HR,main="Emergence Rate - High",
     breaks=round(seq(0,0.3,length.out=10),digits=2),col=terrain.colors(10,rev=T))
POSE_emerg_W<-matrix(POSE_agg2$emerg_rate[which(POSE_agg2$Treatment=="W")],ncol=7,byrow=T)
POSE_emerg_WR<-raster(POSE_emerg_W,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_emerg_WR,main="Emergence - Water",
     breaks=round(seq(0,0.3,length.out=10),digits=2),col=terrain.colors(10,rev=T))
POSE_emerg_L<-matrix(POSE_agg2$emerg_rate[which(POSE_agg2$Treatment=="L")],ncol=7,byrow=T)
POSE_emerg_LR<-raster(POSE_emerg_L,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_emerg_LR,main="Emergence - Low",
     breaks=round(seq(0,0.3,length.out=10),digits=2),col=terrain.colors(10,rev=T))
POSE_emerg_C<-matrix(POSE_agg2$emerg_rate[which(POSE_agg2$Treatment=="C")],ncol=7,byrow=T)
POSE_emerg_CR<-raster(POSE_emerg_C,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(POSE_emerg_CR,main="Emergence - Control",
     breaks=round(seq(0,0.3,length.out=10),digits=2),col=terrain.colors(10,rev=T))



