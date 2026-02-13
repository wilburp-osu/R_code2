######################################
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
library(MuMIn)


######################################
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

############################### TEM comments
# temp$unique_ID<-paste0(temp$Group,"_",temp$Treatment,"_",temp$Cell)
# add unqiue ID into the clean data dataframes.
# remove all columns not needed
# what columns we need are... c("Group","days_til_emerg",
#         "days_til_germ","Treatment","Day_Scale_7","Night_Scale_7","unique_ID)
# then unique(dataframe)
# should have a dataframe with one row for each seed in the end 784 rows for each species
# add back in (1|Group)

### aggregate data for each cuvette for each day ######
view(clean.data.list)
str(clean.data.list)

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
                       Night_Scale_7 , 
                     data = short.temp, 
                     mean, na.rm = TRUE)
  
  clean.data.list.short[[s]] <- indiv
  
}

##can we do data analysis by taking each individual out
# that is the idea behind the "for those that did emerge..."
#but then for loop wouldnt work because unequal lengths

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

##Run stats to prep for dredge

rootDaysfit.list <-list()
rootDaysfit.full.list <- list()
for (s in species){
  print(s)
  temp_agg2 <- short.data.list2[[s]]
  
  #pp<-ggplot(temp_agg2,aes(y=days_til_germ,x=as.factor(Night_Scale_7),fill=Treatment))+
   # geom_boxplot()+
    #ggtitle(paste0(s,"_germ_Night"))
  #print(pp)
  
  #pp<-ggplot(temp_agg2,aes(y=days_til_germ,x=as.factor(Day_Scale_7),fill=Treatment))+
   # geom_boxplot()+
    #ggtitle(paste0(s,"_germ_Day"))
  #print(pp)
  
  rootDaysfull<-glmmTMB(days_til_germ~species*(Day_Scale_7+Night_Scale_7+
                          Treatment),
                        family=nbinom2(link = "log"),ziformula=~1,data=temp_agg2,
                        control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))
  print(summary(rootDaysfull))
  # 
  # ##why both named rootDaysfit? should I rename?## ## I corrected this (good catch)
  # 
  rootDaysfit<-glmmTMB(days_til_germ~Day_Scale_7+Night_Scale_7+
                         Day_Scale_7:Night_Scale_7,
                       family=nbinom2(link = "log"),ziformula=~1,data=temp_agg2[which(temp_agg2$Treatment=="C"),],
                       control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))
  print(summary(rootDaysfit))
  
  # rootDaysfit<-glmmTMB(days_til_germ~Night_Scale_7+Day_Scale_7,
  #                      family=nbinom2(link = "log"),ziformula=~1,data=temp_agg2[which(temp_agg2$Treatment=="C"),],
  #                      control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))
  # 
  #                      
  #                      
  # print(summary(rootDaysfit))
  # 
  if (s=="ARTR") {
    rootDaysfit<-glmmTMB(days_til_germ~Day_Scale_7+Night_Scale_7,
                         family=nbinom2(link = "log"),ziformula=~1,data=temp_agg2[which(temp_agg2$Treatment=="C"),],
                         control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))
    print(summary(rootDaysfit))
    
  }
  # 
  # if (s=="ARTR") {
  #   rootDaysfit<-glmmTMB(I(germ_rate+0.0001)~(1|Group)+Day_Temp_Scale+Night_Temp_Scale,
  #                        family=Gamma(link="log"),data=temp_agg2)
  #   print(summary(rootDaysfit))
  # }
  # 
  rootDaysfit.list[[s]] <- rootDaysfit
  rootDaysfit.full.list[[s]] <- rootDaysfull
  
}

###DREDGE TIME!
#Germ
options(na.action = na.omit)
germ.dredge.list <- list()

for (s in species){
  print(s)
  global.model <- short.data.list2[[s]]
  
  dredge.models <- dredge(global.model)
  
  print(dredge.models)
  
  germ.dredge.list[[s]] <- dredge.models
  
}
warnings()


