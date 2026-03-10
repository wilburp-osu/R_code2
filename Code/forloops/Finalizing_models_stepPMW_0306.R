rm(list=ls()) # clear workspace
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

######################################
#create open list
data.list<-list()
#create species names
species<-c("POSE","ARTR","ACMI","ELEL")

##set wd
setwd("C:/Users/18034/Dropbox/PC/Desktop/R_code/R_code2")

# Read in data ####
for (s in species ) {
  raw.read<-read.csv(paste0("Data/",s,"_Data/",s,"_Data_FINAL.csv"))
  if (any(is.na(raw.read$Group))){
    raw.read2<-raw.read[-which(is.na(raw.read$Group)),]
  } else {raw.read2<-raw.read}
  data.list[[s]]<-raw.read2
}

#str(data.list)
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
#view(clean.data.list)
#str(clean.data.list)

clean.data.list.short<-list()

for (s in species){
  
  temp2 <-clean.data.list[[s]]
  
  # add unqiue ID into the clean data dataframes. 
  temp2$unique_ID<-paste0(temp2$Group,"_",temp2$Treatment,"_",temp2$Cell)
  
  # remove all columns not needed
  # what columns we need are... c("Group","days_til_emerg",
  #         "days_til_germ","Treatment","Day_Scale_7","Night_Scale_7","unique_ID)
  
  short.temp <- temp2[, c("Group",
                          "days_til_emerg",
                          "days_til_germ", 
                          "Treatment",
                          "Cell",
                          "Day_Scale_7", 
                          "Night_Scale_7", 
                          "unique_ID",
                          "germ_rate",
                          "emerg_rate")]
  
  #instead of unique() which would create a seed for each day
  
  
  indiv <- aggregate(cbind(days_til_emerg, 
                           days_til_germ,
                           germ_rate,
                           emerg_rate)
                     ~unique_ID + 
                       Group + 
                       Treatment +  
                       Cell + 
                       Day_Scale_7 +
                       Night_Scale_7, 
                     data = short.temp, 
                     mean, na.rm = TRUE)
  
  clean.data.list.short[[s]] <- indiv
  
}

short.data.list2 <- list()

for (s in species){
  temp3<-clean.data.list.short[[s]]
  
  ###  calculate days til emergence and germination ####
  
  for (i in 1:max(temp3$Group)){
    for (t in unique(temp3$Treatment)) {
      for (c in 1:4) {
        
        temp3$emerg_rate[temp3$emerg_rate == 0] <- NA
        temp3$germ_rate[temp3$germ_rate == 0] <- NA
        temp3$days_til_emerg[temp3$days_til_emerg == 0] <- NA
        temp3$days_til_germ[temp3$days_til_germ == 0] <- NA
        
      }
    }
  }
  
  short.data.list2[[s]]<-temp3
  
}

######################################
## germination glm
#should i be using (1|Group)?

#parameters from step

rootDaysfit.full.list <- list()
rootDays.output <- list()

for (s in species){
  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  germ_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","days_til_germ") #select the columns I want for the germ model
  temp_agg2_germ<-temp_agg2[,germ_cols] # make a new dataframe with just those values
  
  rootDaysfull<-glm(days_til_germ~Treatment+Day_Scale_7+Night_Scale_7,
                    family=poisson(link = "log"),data=temp_agg2_germ)
  print(summary(rootDaysfull))

  if (s=="ACMI") {
     rootDaysfull<-glm(days_til_germ~Day_Scale_7+Night_Scale_7+Treatment+
                        Day_Scale_7:Night_Scale_7,
                      family=poisson(link = "log"),data=temp_agg2_germ)
     print(summary(rootDaysfull))
  }
  
   if (s=="POSE") {
     rootDaysfull<-glm(days_til_germ~Day_Scale_7+Night_Scale_7+Treatment+
                        Day_Scale_7:Treatment,
                      family=poisson(link = "log"),data=temp_agg2_germ)
     print(summary(rootDaysfull))
   }
  
  if (s=="ARTR") {
    rootDaysfull<-glm(days_til_germ~Day_Scale_7+Night_Scale_7+Treatment+
                        Day_Scale_7:Treatment,
                      family=poisson(link = "log"),data=temp_agg2_germ)
    print(summary(rootDaysfull))
  }
  rootDaysfit.full.list[[s]] <- rootDaysfull
  
  rootDaysfull$coefficients
  
  summary(rootDaysfull)$coefficients
  
  rootDaysrast <- data.frame(summary(rootDaysfull)$coefficients)
  rootDaysrast$Species <- s
  rootDaysrast$coefficients <- rownames(rootDaysrast)
  rootDays.output[[s]]<-rootDaysrast

}

output.all.rootdays <- do.call(rbind, rootDays.output)

names(output.all.rootdays)

output.all.rootdays$sigvar <-
  ifelse(output.all.rootdays$Pr...z.. < 0.05, 1, 0)

rootDays_sig<-output.all.rootdays[which(output.all.rootdays$sigvar==1),]         

ggplot(output.all.rootdays[-which(output.all.rootdays$coefficients=="(Intercept)"),],
       aes(x=Species,y=coefficients))+
  geom_raster(aes(fill=Estimate))+
  geom_tile(data =rootDays_sig[-which(rootDays_sig$coefficients=="(Intercept)"),], 
            aes(x=Species,y=coefficients, fill = sigvar),fill="transparent",
            col="black",size=2)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, space = "Lab",
                       name="Effect")+
  geom_text(aes(label = paste0(round(Estimate,2))))+
  ylab("Variables & Interactions")+
  xlab("Species")+
  ggtitle("Temperature & Treatments Interactions with Days Until Germination")+
  theme_grey()

## emergence glm

#parameters from step

shootDaysfit.full.list <- list()
shootDays.output <- list()

for (s in species){
  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  emerg_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","days_til_emerg") #select the columns I want for the germ model
  temp_agg2_emerg<-temp_agg2[,emerg_cols] # make a new dataframe with just those values
  
  shootDaysfull<-glm(days_til_emerg~Treatment+Day_Scale_7+Night_Scale_7,
                    family=poisson(link = "log"),data=temp_agg2_emerg)
  print(summary(shootDaysfull))
  
  shootDaysfit.full.list[[s]] <- shootDaysfull
  
  shootDaysfull$coefficients
  
  summary(rootDaysfull)$coefficients
  
  shootDaysrast <- data.frame(summary(shootDaysfull)$coefficients)
  shootDaysrast$Species <- s
  shootDaysrast$coefficients <- rownames(shootDaysrast)
  shootDays.output[[s]]<-shootDaysrast
  
}

output.all.shootdays <- do.call(rbind, shootDays.output)

output.all.shootdays$sigvar <-
  ifelse(output.all.shootdays$Pr...z.. < 0.05, 1, 0)


shootDays_sig<-output.all.shootdays[which(output.all.shootdays$sigvar==1),]         

ggplot(output.all.shootdays[-which(output.all.shootdays$coefficients=="(Intercept)"),],
       aes(x=Species,y=coefficients))+
  geom_raster(aes(fill=Estimate))+
  geom_tile(data =shootDays_sig[-which(shootDays_sig$coefficients=="(Intercept)"),], 
            aes(x=Species,y=coefficients, fill = sigvar),fill="transparent",
            col="black",size=2)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, space = "Lab",
                       name="Effect")+
  geom_text(aes(label = paste0(round(Estimate,2))))+
  ylab("Variables & Interactions")+
  xlab("Species")+
  ggtitle("Temperature & Treatments Interactions with Days Until Emergence") +
  theme_grey()

##Rates
##germ

root.rate.full.list <- list()
root.rate.output <- list()

for (s in species) {

  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  temp_agg3<-temp_agg2[-which(temp_agg2$germ_rate==0),]
  grate_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","germ_rate") #select the columns I want for the germ model
  temp_agg3_grate<-temp_agg2[,grate_cols] # make a new dataframe with just those values
  
  rootratefit.full<-glm(germ_rate~Day_Scale_7+Night_Scale_7+Treatment,
                            family=Gamma(link="log"),data=temp_agg3_grate)
  print(summary(rootratefit.full))
  
  if (s=="POSE"){
    rootratefit.full<-glm(germ_rate~Day_Scale_7+Night_Scale_7+Treatment+
                                Day_Scale_7:Treatment,
                              family=Gamma(link="log"),data=temp_agg3_grate)
    print(summary(rootratefit.full)) 
  }
  
  if (s=="ACMI"){
    rootratefit.full<-glm(germ_rate~Day_Scale_7+Night_Scale_7+Treatment+
                                Day_Scale_7:Night_Scale_7,
                              family=Gamma(link="log"),data=temp_agg3_grate)
    print(summary(rootratefit.full)) 
  }
  
  if (s=="ELEL"){
    rootratefit.full<-glm(germ_rate~Day_Scale_7+Night_Scale_7+Treatment+
                                Day_Scale_7:Treatment+Day_Scale_7:Night_Scale_7,
                              family=Gamma(link="log"),data=temp_agg3_grate)
    print(summary(rootratefit.full)) 
  }
  
  root.rate.full.list[[s]] <- rootratefit.full
  
  rootratefit.full$coefficients
  
  summary(rootratefit.full)$coefficients
  
  rootraterast <- data.frame(summary(rootratefit.full)$coefficients)
  rootraterast$Species <- s
  rootraterast$coefficients <- rownames(rootraterast)
  root.rate.output[[s]]<-rootraterast

}

output.all.rootrate <- do.call(rbind, root.rate.output)

output.all.rootrate$sigvar <-
  ifelse(output.all.rootrate$Pr...t.. < 0.05, 1, 0)

rootrate_sig<-output.all.rootrate[which(output.all.rootrate$sigvar==1),]         

ggplot(output.all.rootrate[-which(output.all.rootrate$coefficients=="(Intercept)"),],
       aes(x=Species,y=coefficients))+
  geom_raster(aes(fill=Estimate))+
  geom_tile(data =rootrate_sig[-which(rootrate_sig$coefficients=="(Intercept)"),], 
            aes(x=Species,y=coefficients, fill = sigvar),fill="transparent",
            col="black",size=2)+
  scale_fill_gradient2(low = "red", high = "blue", mid = "white",
                       midpoint = 0, space = "Lab",
                       name="Effect")+
  geom_text(aes(label = paste0(round(Estimate,2))))+
  ylab("Variables & Interactions")+
  xlab("Species")+
  ggtitle("Temperature & Treatments Interactions with Germination Rate")+
  theme_grey()

##emerg

shoot.rate.full.list <- list()
shoot.rate.output <- list()

for (s in species) {
  
  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  temp_agg3<-temp_agg2[-which(temp_agg2$emerg_rate==0),]
  erate_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","emerg_rate") #select the columns I want for the germ model
  temp_agg3_erate<-temp_agg2[,erate_cols] # make a new dataframe with just those values
  
  shootratefit.full<-glm(emerg_rate~Day_Scale_7+Night_Scale_7+Treatment,
                        family=Gamma(link="log"),data=temp_agg3_erate)
  print(summary(shootratefit.full))
  
  if (s=="POSE"){
    shootratefit.full<-glm(emerg_rate~Day_Scale_7+Night_Scale_7+Treatment+
                            Night_Scale_7:Treatment,
                          family=Gamma(link="log"),data=temp_agg3_erate)
    print(summary(shootratefit.full)) 
  }
  
  if (s=="ELEL"){
    shootratefit.full<-glm(emerg_rate~Day_Scale_7+Night_Scale_7+Treatment+
                            Day_Scale_7:Treatment+Day_Scale_7:Night_Scale_7,
                          family=Gamma(link="log"),data=temp_agg3_erate)
    print(summary(shootratefit.full)) 
  }
  
  shoot.rate.full.list[[s]] <- shootratefit.full
  
  shootratefit.full$coefficients
  
  summary(shootratefit.full)$coefficients
  
  shootraterast <- data.frame(summary(shootratefit.full)$coefficients)
  shootraterast$Species <- s
  shootraterast$coefficients <- rownames(shootraterast)
  shoot.rate.output[[s]]<-shootraterast
  
  
}

output.all.shootrate <- do.call(rbind, shoot.rate.output)

output.all.shootrate$sigvar <-
  ifelse(output.all.shootrate$Pr...t.. < 0.05, 1, 0)

shootrate_sig<-output.all.shootrate[which(output.all.shootrate$sigvar==1),]         

ggplot(output.all.shootrate[-which(output.all.shootrate$coefficients=="(Intercept)"),],
       aes(x=Species,y=coefficients))+
  geom_raster(aes(fill=Estimate))+
  geom_tile(data =shootrate_sig[-which(shootrate_sig$coefficients=="(Intercept)"),], 
            aes(x=Species,y=coefficients, fill = sigvar),fill="transparent",
            col="black",size=2)+
  scale_fill_gradient2(low = "red", high = "blue", mid = "white",
                       midpoint = 0, space = "Lab",
                       name="Effect")+
  geom_text(aes(label = paste0(round(Estimate,2))))+
  ylab("Variables & Interactions")+
  xlab("Species")+
  ggtitle("Temperature & Treatments Interactions with Germination Rate")+
  theme_grey()


###########
##days til germ
#predict.list.root <- list()

#for (s in species){
#  rootDaysfit <- rootDaysfit.full.list[[s]]
#  #temp_agg<- short.data.list2[[s]]
#  pnew.data.root <- short.data.list2[[s]]
# #pnew.data<-temp_agg[[,c(1,2,3,4,8,9)]]
  
#  pnew.data.root$Prediction<-predict(rootDaysfull,pnew.data.root,allow.new.levels=T)
#  pnew.data.root$Prediction_back<-exp(pnew.data.root$Prediction)
#  print(ggplot(pnew.data.root,aes(x=days_til_germ,y=Prediction_back,col=Treatment))+
#          geom_point()+geom_abline(slope=1)+theme_minimal()+ labs(title = paste0(s)))
#  
#  predict.list.root[[s]] <- pnew.data.root
#

### raster figs ####
## is this right?
##for (s in species){
#  print(s)
#  pnew.data.root <- predict.list.root[[s]]
#  temp_agg2 <- short.data.list2[[s]]
#  
#  #temp_agg2<-temp_agg[-which(temp_agg$Avg_Root_rate==0),]
#  mid.pointss<-mean(temp_agg2$days_til_germ[which(temp_agg2$Treatment=="C")],na.rm=T)
#  
#  print(mid.pointss)
  
#  germplots <- ggplot(pnew.data.root,aes(x=Day_Scale_7,y=Night_Scale_7))+
#    geom_raster(aes(fill=Prediction_back))+
#    facet_grid(.~Treatment)+
#    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                         midpoint = mid.pointss, space = "Lab",
#                         name="Days until Germination") +
#    theme_minimal()+
#    ylab("Night Temperature (C)")+
#    xlab("Day Temperature (C)")+
#    theme(legend.position = "bottom")+
#    labs(title = paste0(s))
#  
#  #CairoPNG(file.path(paste0("Figures/",s,"_Root_growth_rate.png")), width = 12, height = 7, units = "in",
  #res = 300)
#  print(germplots)
  # dev.off()
#}   


##days til emerg
#predict.list.shoot <- list()

### 2D figs #####
#for (s in species){
#  shootDaysfit <- shootDaysfit.full.list[[s]]
  #temp_agg<- short.data.list2[[s]]
#  pnew.data.shoot <- short.data.list2[[s]]
  #pnew.data<-temp_agg[[,c(1,2,3,4,8,9)]]
  
#  pnew.data.shoot$Prediction<-predict(shootDaysfull,pnew.data.shoot,allow.new.levels=T)
#  pnew.data.shoot$Prediction_back<-exp(pnew.data.shoot$Prediction)
#  print(ggplot(pnew.data.shoot,aes(x=days_til_emerg,y=Prediction_back,col=Treatment))+
#          geom_point()+geom_abline(slope=1)+theme_minimal()+ labs(title = paste0(s)))
  
#  predict.list.shoot[[s]] <- pnew.data.shoot
#}

### raster figs ####
## is this right? look way different than original
#for (s in species){
#  print(s)
#  pnew.data.shoot <- predict.list.shoot[[s]]
#  temp_agg2 <- short.data.list2[[s]]
#  
  #temp_agg2<-temp_agg[-which(temp_agg$Avg_Root_rate==0),]
#  mid.pointss<-mean(temp_agg2$days_til_emerg[which(temp_agg2$Treatment=="C")],na.rm=T)
  
#  print(mid.pointss)
  
#  emergplots <- ggplot(pnew.data,aes(x=Day_Scale_7,y=Night_Scale_7))+
#    geom_raster(aes(fill=Prediction_back))+
#    facet_grid(.~Treatment)+
#    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                         midpoint = mid.pointss, space = "Lab",
#                         name="Days until Emergence") +
#    theme_minimal()+
#    ylab("Night Temperature (C)")+
#    xlab("Day Temperature (C)")+
#    theme(legend.position = "bottom")+
#    labs(title = paste0(s))
  
  #CairoPNG(file.path(paste0("Figures/",s,"_Root_growth_rate.png")), width = 12, height = 7, units = "in",
           #res = 300)
#  print(emergplots)
 # dev.off()
#}   

##########

#predict.list.grate <- list()
##Rates
#for (s in species){
#  rootratefit <- root.rate.full.list[[s]]
  #temp_agg<-agg.list[[s]]
#  pnew.data.grate <-  short.data.list2[[s]]
  #pnew.data<-temp_agg[[,c(1,2,3,4,8,9)]]
  
#  pnew.data.grate$Prediction<-predict(rootratefit,pnew.data.grate,allow.new.levels=T)
#  pnew.data.grate$Prediction_back<-exp(pnew.data.grate$Prediction)
#  print(ggplot(pnew.data.grate,aes(x=germ_rate,y=Prediction_back,col=Treatment))+
#          geom_point()+geom_abline(slope=1)+theme_minimal()+ labs(title = paste0(s)))
  ##no line?
#  predict.list.grate[[s]] <- pnew.data.grate
#}

### raster figs ####
#for (s in species){
#  print(s)
#  pnew.data.grate <- predict.list.grate[[s]]
#  temp_agg2 <- short.data.list2[[s]]
  
  #temp_agg2<-temp_agg[-which(temp_agg$Avg_Root_rate==0),]
#  mid.pointss<-mean(temp_agg2$germ_rate[which(temp_agg2$Treatment=="C")],na.rm=T)
  
#  print(mid.pointss)
  
#  plots <- ggplot(pnew.data.grate,aes(x=Day_Scale_7,y=Night_Scale_7))+
#    geom_raster(aes(fill=Prediction_back))+
#    facet_grid(.~Treatment)+
#    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                         midpoint = mid.pointss, space = "Lab",
##                         name="Root growth rate (final root length/days of growth)") +
#    theme_minimal()+
#    ylab("Night Temperature (C)")+
#    xlab("Day Temperature (C)")+
#    theme(legend.position = "bottom")+
#    labs(title = paste0(s))
  
  #CairoPNG(file.path(paste0("Figures/",s,"_Root_growth_rate.png")), width = 12, height = 7, units = "in",
           #res = 300)
#  print(plots)
  #dev.off()
#}   

###

#predict.list.grate <- list()
##Rates
#for (s in species){
#  rootratefit <- root.rate.full.list[[s]]
  #temp_agg<-agg.list[[s]]
#  pnew.data.grate <-  short.data.list2[[s]]
  #pnew.data<-temp_agg[[,c(1,2,3,4,8,9)]]
  
#  pnew.data.grate$Prediction<-predict(rootratefit,pnew.data.grate,allow.new.levels=T)
#  pnew.data.grate$Prediction_back<-exp(pnew.data.grate$Prediction)
#  print(ggplot(pnew.data.grate,aes(x=germ_rate,y=Prediction_back,col=Treatment))+
#          geom_point()+geom_abline(slope=1)+theme_minimal()+ labs(title = paste0(s)))
  
#  predict.list.grate[[s]] <- pnew.data.grate
#}

### raster figs ####
#for (s in species){
#  print(s)
#  pnew.data.grate <- predict.list.grate[[s]]
#  temp_agg2 <- short.data.list2[[s]]
#  
  #temp_agg2<-temp_agg[-which(temp_agg$Avg_Root_rate==0),]
#  mid.pointss<-mean(temp_agg2$germ_rate[which(temp_agg2$Treatment=="C")],na.rm=T)
  
#  print(mid.pointss)
  
#  plots <- ggplot(pnew.data.grate,aes(x=Day_Scale_7,y=Night_Scale_7))+
#    geom_raster(aes(fill=Prediction_back))+
#    facet_grid(.~Treatment)+
#    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                         midpoint = mid.pointss, space = "Lab",
#                         name="Root growth rate (final root length/days of growth)") +
#    theme_minimal()+
#    ylab("Night Temperature (C)")+
#    xlab("Day Temperature (C)")+
#    theme(legend.position = "bottom")+
#    labs(title = paste0(s))
  
  #CairoPNG(file.path(paste0("Figures/",s,"_Root_growth_rate.png")), width = 12, height = 7, units = "in",
  #res = 300)
#  print(plots)
  #dev.off()
#}   

###
#Emerg
#predict.list.erate <- list()
##Rates
#for (s in species){
#  shootratefit <- shoot.rate.full.list[[s]]
  #temp_agg<-agg.list[[s]]
#  pnew.data.erate <-  short.data.list2[[s]]
  #pnew.data<-temp_agg[[,c(1,2,3,4,8,9)]]
  
#  pnew.data.erate$Prediction<-predict(shootratefit,pnew.data.erate,allow.new.levels=T)
# pnew.data.erate$Prediction_back<-exp(pnew.data.erate$Prediction)
#  print(ggplot(pnew.data.erate,aes(x=emerg_rate,y=Prediction_back,col=Treatment))+
#          geom_point()+geom_abline(slope=1)+theme_minimal()+ labs(title = paste0(s)))
  
#  predict.list.erate[[s]] <- pnew.data.erate
#}

### raster figs ####
#for (s in species){
#  print(s)
#  pnew.data.erate <- predict.list.erate[[s]]
#  temp_agg2 <- short.data.list2[[s]]
  
  #temp_agg2<-temp_agg[-which(temp_agg$Avg_Root_rate==0),]
#  mid.pointss<-mean(temp_agg2$emerg_rate[which(temp_agg2$Treatment=="C")],na.rm=T)
  
#  print(mid.pointss)
  
#  plots <- ggplot(pnew.data.erate,aes(x=Day_Scale_7,y=Night_Scale_7))+
#    geom_raster(aes(fill=Prediction_back))+
#    facet_grid(.~Treatment)+
#    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                         midpoint = mid.pointss, space = "Lab",
#                         name="Shoot growth rate (final shoot length/days of growth)") +
#    theme_minimal()+
#    ylab("Night Temperature (C)")+
#    xlab("Day Temperature (C)")+
#    theme(legend.position = "bottom")+
#    labs(title = paste0(s))
  
  #CairoPNG(file.path(paste0("Figures/",s,"_Root_growth_rate.png")), width = 12, height = 7, units = "in",
  #res = 300)
#  print(plots)
  #dev.off()
#}   


