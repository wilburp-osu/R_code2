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


#create open list
data.list<-list()
#create species names
species<-c("POSE","ARTR","ACMI","ELEL")

# Read in data ####
for (s in species ) {
  raw.read<-read.csv(paste0("Data/",s,"_Data/",s,"_Data_FINAL.csv"))
  if (any(is.na(raw.read$Group))){
    raw.read2<-raw.read[-which(is.na(raw.read$Group)),]
  } else {raw.read2<-raw.read}
  data.list[[s]]<-raw.read2
}

#days<-unique(POSE$Day_Collected)

str(data.list)
clean.data.list<-list()

for (s in species){
  temp<-data.list[[s]]
  
  ## clean data ####
  temp$Root_Pheno_Code<-as.numeric(paste0(temp$Root_Pheno_Code))
  temp$Shoot_Pheno_Code<-as.numeric(paste0(temp$Shoot_Pheno_Code))
  
  temp$germ_rate<-NA
  temp$emerg_rate<-NA
  temp$days_til_germ<-NA
  temp$days_til_emerg<-NA
  
  ###  calculate days til emergence and germination ####
  
  for (i in 1:max(temp$Group)){
    for (t in unique(temp$Treatment)) {
      for (c in 1:4) {
        temp_seed<-temp[which(temp$Group==i & temp$Treatment==t & temp$Cell==c),]
        if (any(temp_seed$Root_Pheno_Code==1)){
          Days<-temp_seed$Days_Since_Install[which(temp_seed$Root_Pheno_Code==1)]
          if (length(Days)>1){
            print(paste0(c("1s",s,Days,i,t,c),collapse="_"))} #print errors for doubling up 1's
          temp$germ_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-1/Days
          temp$days_til_germ[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-Days
        } else {
          temp$germ_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-0
          temp$days_til_germ[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-0
        }
        if (any(temp_seed$Shoot_Pheno_Code==6)){
          Days2<-temp_seed$Days_Since_Install[which(temp_seed$Shoot_Pheno_Code==6)]
          if (length(Days2)>1){
            print(paste0(c("6s",s,Days2,i,t,c),collapse="_"))} #print errors for doubling up 6's
          temp$emerg_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-1/Days2  
          temp$days_til_emerg[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-Days2  
        } else {
          temp$emerg_rate[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-0
          temp$days_til_emerg[which(temp$Group==i & temp$Treatment==t & temp$Cell==c)]<-0
        }
      }
    }
  }
  
  clean.data.list[[s]]<-temp
  
}


### aggregate data for each cuvette for each day ######

str(clean.data.list)
agg.list<- list()

for (s in species) {
  temp <-clean.data.list[[s]]
  temp_agg<-aggregate(cbind(temp$Shoot_Pheno_Code,
                            temp$Root_Pheno_Code,
                            temp$germ_rate,
                            temp$emerg_rate,
                            temp$days_til_emerg,
                            temp$days_til_germ),
                      by=list(temp$Day_Scale_7,temp$Night_Scale_7,
                              temp$Treatment,
                              temp$Days_Since_Install,
                              temp$Group),mean,na.rm=T)
  
  temp_agg[temp_agg =="NaN"]<-NA
  names(temp_agg)<-c("Day_Scale_7","Night_Scale_7","Treatment",
                     "Days_Since_Install","Group",
                     "Root_Pheno_Code","Shoot_Pheno_Code",
                     "germ_rate","emerg_rate",
                     "days_til_emerg","days_til_germ")
  agg.list[[s]] <- temp_agg
}


### aggregate data across all days #####

agg.list2<-list()

for (s in species){
  temp_agg <- agg.list[[s]]
  temp_agg2<-aggregate(cbind(temp_agg$germ_rate,
                             temp_agg$emerg_rate,
                             temp_agg$days_til_germ,
                             temp_agg$days_til_emerg),
                       by=list(temp_agg$Day_Scale_7,temp_agg$Night_Scale_7,
                               temp_agg$Treatment,
                               temp_agg$Group),mean,na.rm=T)
  
  temp_agg2[temp_agg2 =="NaN"]<-NA
  names(temp_agg2)<-c("Day_Scale_7","Night_Scale_7","Treatment",
                      "Group","germ_rate","emerg_rate",
                      "days_til_germ","days_til_emerg")
  
  #temp_agg2$Day_Temp_Scale<-scale(temp_agg2$Day_Temp)
  #temp_agg2$Night_Temp_Scale<-scale(temp_agg2$Night_Temp)
  
  agg.list2[[s]] <- temp_agg2
  
}

# germination glmmTMB ####

# should I put this in and change the first section: rootDaysfit.full.list <- list()?

rootDaysfit.list <-list()
rootDaysfit.full.list <- list()
for (s in species){
  print(s)
  temp_agg2 <- clean.data.list[[s]]
  rootDaysfull<-glmmTMB(days_til_germ~(1|Group)+Day_Scale_7+Night_Scale_7+
                          Day_Scale_7:Treatment+Night_Scale_7:Treatment+
                          Treatment+Day_Scale_7:Night_Scale_7+
                          Day_Scale_7:Night_Scale_7:Treatment,
                        family=nbinom2(link = "log"),ziformula=~1,data=temp_agg2)
  print(summary(rootDaysfull))
  # 
  # ##why both named rootDaysfit? should I rename?## ## I corrected this (good catch)
  # 
  # rootDaysfit<-glmmTMB(I(germ_rate+0.0001)~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
  #                        Treatment+Day_Temp_Scale:Night_Temp_Scale,
  #                      family=Gamma(link="log"),data=temp_agg2)
  # print(summary(rootDaysfit))
  # 
  # if (s=="ACMI") {
  #   rootDaysfit<-glmmTMB(I(germ_rate+0.0001)~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
  #                          Day_Temp_Scale:Night_Temp_Scale,
  #                        family=Gamma(link="log"),data=temp_agg2)
  #   print(summary(rootDaysfit))
  # }
  # 
  # if (s=="ARTR") {
  #   rootDaysfit<-glmmTMB(I(germ_rate+0.0001)~(1|Group)+Day_Temp_Scale+Night_Temp_Scale,
  #                        family=Gamma(link="log"),data=temp_agg2)
  #   print(summary(rootDaysfit))
  # }
  # 
  # rootDaysfit.list[[s]] <- rootDaysfit
   rootDaysfit.full.list[[s]] <- rootDaysfull
  
}

## make plots for germiation #####

#str(agg.list)
#str(agg.list2)
prediction.list<- list()

for (s in species){
  
  rootDaysfit <- rootDaysfit.full.list[[s]]
  #temp_agg<-agg.list2[[s]]
  pnew.data<-agg.list2[[s]]
  
  pnew.data$Prediction<-predict(rootDaysfit,pnew.data)
  pnew.data$Prediction_back<-exp(pnew.data$Prediction)
  #print(ggplot(pnew.data,aes(x=days_til_germ,y=Prediction_back,col=Treatment))+
  #        geom_point()+geom_abline(slope=1)+theme_minimal()+ labs(title = paste0(s),"_Germination"))
  
  plotly_x<-matrix(pnew.data$Day_Scale_7[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
  plotly_y<-matrix(pnew.data$Night_Scale_7[which(pnew.data$Treatment=="H")],byrow=T,ncol=7)
  plotly_zH<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="H")],byrow = T,ncol=7)
  plotly_zL<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="L")],byrow = T,ncol=7)
  plotly_zW<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="W")],byrow = T,ncol=7)
  plotly_zC<-matrix(pnew.data$Prediction_back[which(pnew.data$Treatment=="C")],byrow = T,ncol=7)
  fig<-plot_ly(showscale=F,type="surface", x=~plotly_x,y=plotly_y) %>% layout(title=paste0(s,"_Germination"))
  fig <- fig %>% add_surface(z = ~plotly_zC, opacity = 0.5, colorscale = list(c(0,1),c("rgb(107,255,184)","rgb(0,124,90)"))) #green
  fig <- fig %>% add_surface(z = ~plotly_zW, opacity = 0.5, colorscale = list(c(0,1),c("rgb(129,212,247)","rgb(12,177,247)"))) #blue
  fig <- fig %>% add_surface(z = ~plotly_zL, opacity = 0.5, colorscale = list(c(0,1),c("rgb(247,226,157)","rgb(255,192,7)"))) #yellow
  fig <- fig %>% add_surface(z = ~plotly_zH, opacity = 0.5, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) #pink
  print(fig)
  
  print(summary(pnew.data$Prediction_back))
  
  prediction.list[[s]] <- pnew.data
  
}

for (s in species){
  print (s)
  pnew.data <- prediction.list[[s]]
  temp_agg2 <- agg.list2[[s]]
  
  for (i in unique(temp_agg2$Group)){
    tempp<-temp_agg2[which(temp_agg2$Group==i),]
    if (all(tempp$days_til_germ==0)) {
      print (i)
    }
  }
  
  mid.pointss<-mean(temp_agg2$days_til_germ[which(temp_agg2$Treatment=="C")])
  print(mid.pointss)
  germplot <- ggplot(pnew.data,aes(x=Day_Scale_7,y=Night_Scale_7))+
    geom_raster(aes(fill=Prediction_back))+ #, interpolate=T
    facet_grid(.~Treatment)+
    scale_fill_gradient2(low = "red", high = "blue", mid = "white",
                         midpoint = mid.pointss, space = "Lab",
                         name="Germination Rate (1/days til germ)") +
    theme_minimal()+
    theme(legend.position = "bottom")+
    labs(title = paste0(s))
  
  print(germplot)
  
}

# emergence glmmTMB models #####

#shootDaysfit.list <-list()
shootDaysfull.list <-list()
for (s in species){
  print(s)
  temp_agg2 <- clean.data.list[[s]]
  #hist(temp_agg2$days_til_emerg)
  shootDaysfull<-glmmTMB(days_til_emerg~(1|Group)+Day_Scale_7+Night_Scale_7+
                          Day_Scale_7:Treatment+Night_Scale_7:Treatment+
                          Treatment+Day_Scale_7:Night_Scale_7+
                          Day_Scale_7:Night_Scale_7:Treatment,
                        family=nbinom2(link = "log"),ziformula=~1,data=temp_agg2,
                        control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))
  print(summary(shootDaysfull))
  # interactions between treatment and temps are not significant so dropped them
  # if (s == "POSE"){
  #   shootDaysfit<-glmmTMB(I(emerg_rate)~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
  #                           Treatment,
  #                         family=ziGamma(link="log"),data=temp_agg3,ziformula=~1)
  #   print(summary(shootDaysfit))
  # }
  # 
  # if (s=="ARTR"){
  #   shootDaysfit<-glmmTMB(I(emerg_rate)~(1|Group)+Day_Temp_Scale+Night_Temp_Scale+
  #                           Treatment,
  #                         family=ziGamma(link="log"),data=temp_agg3,ziformula=~1)
  #   print(summary(shootDaysfit))
  # }
  
  if (s == "ACMI" | s=="ELEL"){
    shootDaysfull<-glmmTMB(days_til_emerg~(1|Group)+Day_Scale_7+Night_Scale_7+
                            Day_Scale_7:Treatment+Night_Scale_7:Treatment+
                            Treatment+Day_Scale_7:Night_Scale_7,
                          family=nbinom2(link = "log"),ziformula=~1,data=temp_agg2,
                          control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))
    print(summary(rootDaysfull))
  }
  
  # if (s == "ELEL") {
  #   shootDaysfit<-glmmTMB(I(emerg_rate)~(1|Group)+Day_Temp_Scale+Night_Temp_Scale +
  #                           Treatment,
  #                         family=ziGamma(link="log"),data=temp_agg3,ziformula=~1)
  #   print(summary(shootDaysfit))
  # }
  
  #shootDaysfit.list[[s]] <- shootDaysfit
  shootDaysfull.list[[s]] <- shootDaysfull
}

###prediction

str(agg.list)
str(agg.list2)
prediction.list.emerg<- list()

## make plots for emergence ######
for (s in species){
  shootDaysfit <- shootDaysfit.list[[s]]
  pnew.data.emerg<-agg.list2[[s]]
  #pnew.data.emerg<-temp_agg2[which(temp_agg2$emerg_rate>0),]
  pnew.data.emerg$Prediction<-predict(shootDaysfit,pnew.data.emerg,allow.new.levels=TRUE)
  pnew.data.emerg$Prediction_back<-exp(pnew.data.emerg$Prediction)
  print(ggplot(pnew.data.emerg,aes(x=temp_agg2$emerg_rate,y=Prediction_back,col=Treatment))+
          geom_point()+geom_abline(slope=1)+theme_minimal()+ labs(title = paste0(s,"_Emergence")))
  
  plotly_x<-matrix(pnew.data.emerg$Day_Temp[which(pnew.data.emerg$Treatment=="H")],byrow=T,ncol=7)
  plotly_y<-matrix(pnew.data.emerg$Night_Temp[which(pnew.data.emerg$Treatment=="H")],byrow=T,ncol=7)
  plotly_zH<-matrix(pnew.data.emerg$Prediction_back[which(pnew.data.emerg$Treatment=="H")],byrow = T,ncol=7)
  plotly_zL<-matrix(pnew.data.emerg$Prediction_back[which(pnew.data.emerg$Treatment=="L")],byrow = T,ncol=7)
  plotly_zW<-matrix(pnew.data.emerg$Prediction_back[which(pnew.data.emerg$Treatment=="W")],byrow = T,ncol=7)
  plotly_zC<-matrix(pnew.data.emerg$Prediction_back[which(pnew.data.emerg$Treatment=="C")],byrow = T,ncol=7)
  fig.emerg<-plot_ly(showscale=F,type="surface", x=~plotly_x,y=plotly_y) %>% layout(title=paste0(s,"_Emergence"))
  fig.emerg <- fig.emerg %>% add_surface(z = ~plotly_zC, opacity = 0.5, colorscale = list(c(0,1),c("rgb(107,255,184)","rgb(0,124,90)"))) #green
  fig.emerg <- fig.emerg %>% add_surface(z = ~plotly_zW, opacity = 0.5, colorscale = list(c(0,1),c("rgb(129,212,247)","rgb(12,177,247)"))) #blue
  fig.emerg <- fig.emerg %>% add_surface(z = ~plotly_zL, opacity = 0.5, colorscale = list(c(0,1),c("rgb(247,226,157)","rgb(255,192,7)"))) #yellow
  fig.emerg <- fig.emerg %>% add_surface(z = ~plotly_zH, opacity = 0.5, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) #pink
  print(fig.emerg)
  
  #print(summary(pnew.data.emerg$Prediction_back))
  
  prediction.list.emerg[[s]] <- pnew.data.emerg
}

for (s in species){
  pnew.data.emerg <- prediction.list.emerg[[s]]
  #temp_agg2 <- agg.list2[[s]]
  mid.pointss<-median(pnew.data.emerg$Prediction_back)
  emerg.plot<- ggplot(pnew.data,aes(x=Day_Temp_Scale,y=Night_Temp_Scale))+
    geom_raster(aes(fill=Prediction_back))+
    facet_grid(.~Treatment)+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = mid.pointss, space = "Lab",
                         name="Emergence Rate (1/days til emergence)") +
    theme_minimal()+
    theme(legend.position = "bottom")+
    labs(title = paste0(s))
  
  print(emerg.plot)
  
}
