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

#Set WD
setwd("C:/Users/18034/Dropbox/PC/Desktop/R_code/R_code2")

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
   temp<-data.list[[s]]
   
### clean data ####
temp$Root_Pheno_Code<-as.numeric(paste0(temp$Root_Pheno_Code))
temp$Shoot_Pheno_Code<-as.numeric(paste0(temp$Shoot_Pheno_Code))

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

clean.data.list[[s]]<-temp

}
###debug

### aggregate data for each cuvette for each day

str(clean.data.list)
agg.list<- list()

for (s in species) {
  temp <-clean.data.list[[s]]
temp_agg<-aggregate(cbind(temp$Shoot_Pheno_Code,
                          temp$Root_Pheno_Code,
                          temp$germ_rate,
                          temp$emerg_rate),
                    by=list(temp$Day_Scale_7,temp$Night_Scale_7,
                            temp$Treatment,
                            temp$Days_Since_Install,
                            temp$Group),mean,na.rm=T)

temp_agg[temp_agg =="NaN"]<-NA
names(temp_agg)<-c("Day_Temp","Night_Temp","Treatment",
                   "Days_Since_Install","Group",
                   "Root_Pheno_Code","Shoot_Pheno_Code",
                   "germ_rate","emerg_rate")
 agg.list[[s]] <- temp_agg
}

#### models for germination and emergence ####
#aggregate across all days 

agg.list2<-list()

for (s in species){
  temp_agg <- agg.list[[s]]
temp_agg2<-aggregate(cbind(temp_agg$germ_rate,
                           temp_agg$emerg_rate),
                     by=list(temp_agg$Day_Temp,temp_agg$Night_Temp,
                             temp_agg$Treatment,
                             temp_agg$Group),mean,na.rm=T)

temp_agg2[temp_agg2 =="NaN"]<-NA
names(temp_agg2)<-c("Day_Temp","Night_Temp","Treatment",
                    "Group","germ_rate","emerg_rate")

temp_agg2$Day_Temp_Scale<-scale(temp_agg2$Day_Temp)
temp_agg2$Night_Temp_Scale<-scale(temp_agg2$Night_Temp)

agg.list2[[s]] <- temp_agg2

}

#### germination glmmTMB ####

# should I put this in and change the first section: rootDaysfit.full.list <- list()?

rootDaysfit.list <-list()

for (s in species){
  temp_agg2 <- agg.list2[[s]]
  rootDaysfit<-glmmTMB(I(germ_rate+0.0001)~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                       Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+
                       Treatment+Day_Temp_Scale:Night_Temp_Scale,
                     family=Gamma(link="log"),data=temp_agg2)
  print(summary(rootDaysfit))
  
##why both named rootDaysfit? should I rename?##
  
  rootDaysfit<-glmmTMB(I(germ_rate+0.0001)~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                       Treatment+Day_Temp_Scale:Night_Temp_Scale,
                     family=Gamma(link="log"),data=temp_agg2)
  print(summary(rootDaysfit))
  
  rootDaysfit.list[[s]] <- rootDaysfit

}

### prediction

str(agg.list)
str(agg.list2)
prediction.list<- list()

for (s in species){
  
  rootDaysfit <- rootDaysfit.list[[s]]
  temp_agg<-agg.list2[[s]]
  pnew.data<-agg.list2[[s]]

  pnew.data$Prediction<-predict(rootDaysfit,pnew.data)
  pnew.data$Prediction_back<-exp(pnew.data$Prediction)
  print(ggplot(pnew.data,aes(x=temp_agg2$germ_rate,y=Prediction_back,col=Treatment))+
    geom_point()+geom_abline(slope=1)+theme_minimal()+ labs(title = paste0(s)))

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
  print(fig)

  print(summary(pnew.data$Prediction_back))

  prediction.list[[s]] <- pnew.data

}

for (s in species){
  pnew.data <- prediction.list[[s]]
germplot <- ggplot(pnew.data,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
  geom_raster(aes(fill=Prediction_back))+ #, interpolate=T
  facet_grid(.~Treatment)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.17, space = "Lab",
                       name="Germination Rate (1/days til germ)") +
  theme_minimal()+
  theme(legend.position = "bottom")+
  labs(title = paste0(s))

print(germplot)
  
}

####emergence glmmTMB####

shootDaysfit.list <-list()

for (s in species){
shootDaysfit<-glmmTMB(I(emerg_rate+0.0001)~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                        Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+
                        Treatment+Day_Temp_Scale:Night_Temp_Scale,
                      family=ziGamma(link="log"),data=temp_agg2,ziformula=~1)
print(summary(shootDaysfit))
# interactions between treatment and temps are not significant so dropped them
shootDaysfit<-glmmTMB(I(emerg_rate+0.0001)~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                        Treatment,
                      family=ziGamma(link="log"),data=temp_agg2,ziformula=~1)
print(summary(shootDaysfit))

shootDaysfit.list[[s]] <- shootDaysfit

}

###prediction

str(agg.list)
str(agg.list2)
prediction.list.emerg<- list()

for (s in species){
  shootDaysfit <- shootDaysfit.list[[s]]
  pnew.data.emerg<-agg.list2[[s]]

  pnew.data.emerg$Prediction<-predict(shootDaysfit,pnew.data.emerg)
  pnew.data.emerg$Prediction_back<-exp(pnew.data.emerg$Prediction)
  print(ggplot(pnew.data.emerg,aes(x=temp_agg2$emerg_rate,y=Prediction_back,col=Treatment))+
    geom_point()+geom_abline(slope=1)+theme_minimal()+ labs(title = paste0(s)))

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
  print(fig.emerg)

  print(summary(pnew.data.emerg$Prediction_back))

  prediction.list.emerg[[s]] <- pnew.data.emerg
}

for (s in species){
  pnew.data.emerg <- prediction.list.emerg[[s]]
  emerg.plot<- ggplot(pnew.data,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
    geom_raster(aes(fill=Prediction_back),interpolate=T)+
    facet_grid(.~Treatment)+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.17, space = "Lab",
                       name="Emergence Rate (1/days til emergence)") +
    theme_minimal()+
    theme(legend.position = "bottom")+
    labs(title = paste0(s))

  print(emerg.plot)

}
