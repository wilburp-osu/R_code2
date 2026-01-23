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
setwd("C:/Users/18034/Dropbox/PC/Desktop/Martyn_Lab/Peter_Thermodatatable_entry/4_SA")

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
        if (any(temp_seed$Root_Pheno_Code==1)){ 
          root_days<-temp_seed$Root_Length[
            which(temp_seed$Root_Pheno_Code>0&temp_seed$Root_Pheno_Code<5)]
          #root_days_dup<-unique(root_days)
          #root_rate<-sum(root_days_dup)/length(root_days)
          root_rate<-max(root_days)/length(root_days)
        } else {root_rate<-0}
        if (any(temp_seed$Shoot_Pheno_Code==6)){
          shoot_days<-(temp_seed$Shoot_Length[
            which(temp_seed$Shoot_Pheno_Code>5&temp_seed$Shoot_Pheno_Code<9)])
          #shoot_days_dup<-unique(shoot_days)
          #shoot_rate<-sum(shoot_days_dup)/length(shoot_days)
          shoot_rate<-max(shoot_days)/length(shoot_days)
        } else {shoot_rate<-0}
        temp$Avg_Root_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-root_rate
        temp$Avg_Shoot_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-shoot_rate
      }
    }
  }
  temp$Root_Pheno_Code[which(temp$Root_Pheno_Code==999)]<-NA
  temp$Shoot_Pheno_Code[which(temp$Shoot_Pheno_Code==999)]<-NA
  
  clean.data.list[[s]]<-temp
}

### aggregate data for each cuvette for each day
str(clean.data.list)
agg.list<- list()

for (s in species) {
  temp <-clean.data.list[[s]]
  temp_agg<-aggregate(cbind(temp$Avg_Root_rate,
                          temp$Avg_Shoot_rate),
                    by=list(temp$Day_Scale_7,temp$Night_Scale_7,
                            temp$Treatment,
                            temp$Group),mean,na.rm=T)

  temp_agg[temp_agg =="NaN"]<-NA
  names(temp_agg)<-c("Day_Temp","Night_Temp","Treatment",
                  "Group","Avg_Root_rate","Avg_Shoot_rate")


  temp_agg$Temp_Ratio<-temp_agg$Day_Temp/temp_agg$Night_Temp
  temp_agg$Treatment<-as.factor(paste0(temp_agg$Treatment))
  temp_agg$Day_Temp_Scale<-scale(temp_agg$Day_Temp)
  temp_agg$Night_Temp_Scale<-scale(temp_agg$Night_Temp)
 agg.list[[s]] <- temp_agg
}

#### root rate glmmTMB ####
#rootratefit<-glmmTMB(Avg_Root_rate~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment,
#                     family=Gamma(link = "log"),data=POSE_agg)

root.rate.full.list <- list()
root.rate.fit.list <-list()

for (s in species) {
  temp_agg<-agg.list[[s]]
  rootratefit.full<-glmmTMB(I(Avg_Root_rate+0.0001)~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                            Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                            Day_Temp_Scale:Night_Temp_Scale+
                            Day_Temp_Scale:Night_Temp_Scale:Treatment,
                          family=Gamma(link="log"),data=temp_agg)

  print(summary(rootratefit.full))
  
  rootratefit<-glmmTMB(I(Avg_Root_rate+0.0001)~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                       Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                       Day_Temp_Scale:Night_Temp_Scale,
                     family=Gamma(link="log"),data=temp_agg)
  
  print(summary(rootratefit))
  root.rate.full.list[[s]] <- rootratefit.full
  root.rate.fit.list[[s]] <- rootratefit
}


### prediction
str(root.rate.full.list)
str(root.rate.fit.list)
str(agg.list)

predict.list <- list()
  
for (s in species){
  rootratefit <- root.rate.fit.list[[s]]
  temp_agg<-agg.list[[s]]
  pnew.data <- agg.list[[s]]
  #pnew.data<-temp_agg[[,c(1,2,3,4,8,9)]]
  
  pnew.data$Prediction<-predict(rootratefit,pnew.data,allow.new.levels=T)
  pnew.data$Prediction_back<-exp(pnew.data$Prediction)
  print(ggplot(pnew.data,aes(x=temp_agg$Avg_Root_rate,y=Prediction_back,col=Treatment))+
    geom_point()+geom_abline(slope=1)+theme_minimal()+ labs(title = paste0(s)))
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
  print(fig)

  print(summary(pnew.data$Prediction_back))
  
  predict.list[[s]] <- pnew.data
  }


for (s in species){
  pnew.data <- predict.list[[s]]
    plots <- ggplot(pnew.data,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
    geom_raster(aes(fill=Prediction_back),interpolate=T)+
    facet_grid(.~Treatment)+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 2.4, space = "Lab",
                         name="Root growth rate (final root length/days of growth)") +
    theme_minimal()+
    theme(legend.position = "bottom")+
      labs(title = paste0(s))
    
 print(plots)
}   

#### shoot rate glmmTMB ####

shoot.rate.full.list <- list()
shoot.rate.fit.list <-list()

for (s in species) {
  temp_agg<-agg.list[[s]]
  shootratefit.full<-glmmTMB(I(Avg_Shoot_rate+0.0001)~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                              Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                              Day_Temp_Scale:Night_Temp_Scale+
                              Day_Temp_Scale:Night_Temp_Scale:Treatment,
                            family=Gamma(link="log"),data=temp_agg)
  
  print(summary(shootratefit.full))
  
  shootratefit<-glmmTMB(I(Avg_Shoot_rate+0.0001)~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
                         Day_Temp_Scale:Treatment+Night_Temp_Scale:Treatment+Treatment+
                         Day_Temp_Scale:Night_Temp_Scale,
                       family=Gamma(link="log"),data=temp_agg)
  
  print(summary(shootratefit))
  shoot.rate.full.list[[s]] <- shootratefit.full
  shoot.rate.fit.list[[s]] <- shootratefit
}

#### prediction ####

predict.list.sh <- list()

for (s in species){
  shootratefit <- shoot.rate.fit.list[[s]]
  shootratefit.full <- shoot.rate.full.list[[s]]
  temp_agg<-agg.list[[s]]
  pnew.data.sh <- agg.list[[s]]

  pnew.data.sh$Prediction<-predict(shootratefit.full,pnew.data.sh)
  pnew.data.sh$Prediction_back<-exp(pnew.data.sh$Prediction)
  print(ggplot(pnew.data.sh,aes(x=temp_agg$Avg_Shoot_rate,y=Prediction_back,col=Treatment))+
    geom_point()+geom_abline(slope=1)+theme_minimal()+ labs(title = paste0(s)))
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
  print(fig)

  summary(pnew.data$Prediction_back)
 predict.list.sh[[s]] <- pnew.data.sh

}

for (s in species){
  pnew.data.sh <- predict.list.sh[[s]]
plot.sh <- ggplot(pnew.data.sh,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
  geom_raster(aes(fill=Prediction_back),interpolate=T)+
  facet_grid(.~Treatment)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 2.4, space = "Lab",
                       name="Shoot growth rate (final root length/days of growth)") +
  theme_minimal()+
  theme(legend.position = "bottom")+
  labs(title = paste0(s))

print(plot.sh)
}

#### plot raster root growth rate ####

temp_root_H.list <- list()
temp_root_HR.list <- list()

temp_root_W.list <- list()
temp_root_WR.list <-list()

temp_root_L.list <- list()
temp_root_LR.list <- list()

temp_root_C.list <- list()
temp_root_CR.list <- list()

for (s in species){
  
  temp_agg<-agg.list[[s]]

  temp_root_H<-matrix(temp_agg$Avg_Root_rate[which(temp_agg$Treatment=="H")],ncol=7,byrow=T)
  temp_root_HR<-raster(temp_root_H,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
  plot(temp_root_HR,main=paste0("Root Growth Rate - High"," (",s,")"),
      breaks=round(seq(0,9,length.out=10),digits=2),col=terrain.colors(10,rev=T))
  temp_root_W<-matrix(temp_agg$Avg_Root_rate[which(temp_agg$Treatment=="W")],ncol=7,byrow=T)
  temp_root_WR<-raster(temp_root_W,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
  plot(temp_root_WR,main= paste0("Root Growth Rate - Water"," (",s,")"),
     breaks=round(seq(0,9,length.out=10),digits=2),col=terrain.colors(10,rev=T))
  temp_root_L<-matrix(temp_agg$Avg_Root_rate[which(temp_agg$Treatment=="L")],ncol=7,byrow=T)
  temp_root_LR<-raster(temp_root_L,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
  plot(temp_root_LR,main= paste0("Root Growth Rate - Low"," (",s,")"),
     breaks=round(seq(0,9,length.out=10),digits=2),col=terrain.colors(10,rev=T))
  temp_root_C<-matrix(temp_agg$Avg_Root_rate[which(temp_agg$Treatment=="C")],ncol=7,byrow=T)
  temp_root_CR<-raster(temp_root_C,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
  plot(temp_root_CR,main= paste0("Root Growth Rate - Control"," (",s,")"),
     breaks=round(seq(0,9,length.out=10),digits=2),col=terrain.colors(10,rev=T))
  
  temp_root_H.list[[s]] <- temp_agg
  temp_root_HR.list[[s]] <- temp_agg
  
  temp_root_W.list[[s]] <- temp_agg
  temp_root_WR.list[[s]] <- temp_agg
  
  temp_root_L.list[[s]] <- temp_agg
  temp_root_LR.list[[s]] <- temp_agg
  
  temp_root_C.list[[s]] <- temp_agg
  temp_root_CR.list[[s]] <- temp_agg

}

#### plot raster shoot growth rate ####
temp_shoot_H.list <- list()
temp_shoot_HR.list <- list()

temp_shoot_W.list <- list()
temp_shoot_WR.list <-list()

temp_shoot_L.list <- list()
temp_shoot_LR.list <- list()

temp_shoot_C.list <- list()
temp_shoot_CR.list <- list()

for (s in species){
  
  temp_agg<-agg.list[[s]]

  temp_shoot_H<-matrix(temp_agg$Avg_Shoot_rate[which(temp_agg$Treatment=="H")],ncol=7,byrow=T)
  temp_shoot_HR<-raster(temp_shoot_H,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
  plot(temp_shoot_HR,main=paste0("Shoot Growth Rate - High"," (",s,")"),
     breaks=round(seq(0,12,length.out=10),digits=2),col=terrain.colors(10,rev=T))
  temp_shoot_W<-matrix(temp_agg$Avg_Shoot_rate[which(temp_agg$Treatment=="W")],ncol=7,byrow=T)
  temp_shoot_WR<-raster(temp_shoot_W,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
  plot(temp_shoot_WR,main=paste0("Shoot Growth Rate - Water"," (",s,")"),
     breaks=round(seq(0,12,length.out=10),digits=2),col=terrain.colors(10,rev=T))
  temp_shoot_L<-matrix(temp_agg$Avg_Shoot_rate[which(temp_agg$Treatment=="L")],ncol=7,byrow=T)
  temp_shoot_LR<-raster(temp_shoot_L,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
  plot(temp_shoot_LR,main=paste0("Shoot Growth Rate - Low"," (",s,")"),
     breaks=round(seq(0,12,length.out=10),digits=2),col=terrain.colors(10,rev=T))
  temp_shoot_C<-matrix(temp_agg$Avg_Shoot_rate[which(temp_agg$Treatment=="C")],ncol=7,byrow=T)
  temp_shoot_CR<-raster(temp_shoot_C,crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
  plot(temp_shoot_CR,main=paste0("Shoot Growth Rate - Control"," (",s,")"),
     breaks=round(seq(0,12,length.out=10),digits=2),col=terrain.colors(10,rev=T))

temp_root_H.list[[s]] <- temp_agg
temp_root_HR.list[[s]] <- temp_agg

temp_root_W.list[[s]] <- temp_agg
temp_root_WR.list[[s]] <- temp_agg

temp_root_L.list[[s]] <- temp_agg
temp_root_LR.list[[s]] <- temp_agg

temp_root_C.list[[s]] <- temp_agg
temp_root_CR.list[[s]] <- temp_agg

}
