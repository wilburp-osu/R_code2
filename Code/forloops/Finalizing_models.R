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

#parameters from dredge

rootDaysfit.full.list <- list()

for (s in species){
  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  germ_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","days_til_germ") #select the columns I want for the germ model
  temp_agg2_germ<-temp_agg2[,germ_cols] # make a new dataframe with just those values
  
  rootDaysfull<-glm(days_til_germ~Day_Scale_7+Night_Scale_7,
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
  
  rootDaysfit.full.list[[s]] <- rootDaysfull
  
}

####
#parameters from Step
####

rootDaysfit.full.list <- list()

for (s in species){
  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  germ_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","days_til_germ") #select the columns I want for the germ model
  temp_agg2_germ<-temp_agg2[,germ_cols] # make a new dataframe with just those values
  
  rootDaysfull<-glm(days_til_germ~Day_Scale_7+Night_Scale_7+Treatment,
                    family=poisson(link = "log"),data=temp_agg2_germ)
  print(summary(rootDaysfull))
  
  if (s=="ACMI") {
    rootDaysfull<-glm(days_til_germ~Day_Scale_7+Night_Scale_7+Treatment+
                        Day_Scale_7:Night_Scale_7+Treatment:Day_Scale_7,
                      family=poisson(link = "log"),data=temp_agg2_germ)
    print(summary(rootDaysfull))
  }
  
  if (s=="POSE") {
    rootDaysfull<-glm(days_til_germ~Day_Scale_7+Night_Scale_7+Treatment+
                        Treatment:Day_Scale_7,
                      family=poisson(link = "log"),data=temp_agg2_germ)
    print(summary(rootDaysfull))
  }
  
  rootDaysfit.full.list[[s]] <- rootDaysfull
  
}
####

## emergence glm

#parameters from dredge

shootDaysfit.full.list <- list()

for (s in species){
  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  emerg_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","days_til_emerg") #select the columns I want for the germ model
  temp_agg2_emerg<-temp_agg2[,emerg_cols] # make a new dataframe with just those values
  
  shootDaysfull<-glm(days_til_emerg~Day_Scale_7+Night_Scale_7,
                    family=poisson(link = "log"),data=temp_agg2_emerg)
  print(summary(shootDaysfull))
  
  if (s=="ACMI") {
    shootDaysfull<-glm(days_til_emerg~Day_Scale_7+Night_Scale_7+Treatment,
                      family=poisson(link = "log"),data=temp_agg2_emerg)
    print(summary(shootDaysfull))
  }
  
  if (s=="POSE") {
    shootDaysfull<-glm(days_til_emerg~Day_Scale_7+Night_Scale_7+Treatment,
                      family=poisson(link = "log"),data=temp_agg2_emerg)
    print(summary(shootDaysfull))
  }
  
  rootDaysfit.full.list[[s]] <- rootDaysfull
  
}


##Rates
##germ

root.rate.full.list <- list()

for (s in species) {

  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  temp_agg3<-temp_agg2[-which(temp_agg2$germ_rate==0),]
  germ_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","germ_rate") #select the columns I want for the germ model
  temp_agg3_germ<-temp_agg2[,germ_cols] # make a new dataframe with just those values
  
  rootratefit.full<-glm(germ_rate~Day_Scale_7+Night_Scale_7,
                            family=Gamma(link="log"),data=temp_agg3_germ)
  print(summary(rootratefit.full))
  
  if (s=="POSE"){
    rootratefit.full<-glm(germ_rate~Day_Scale_7+Night_Scale_7+Treatment+
                                Day_Scale_7:Treatment,
                              family=Gamma(link="log"),data=temp_agg3_germ)
    print(summary(rootratefit.full)) 
  }
  
  if (s=="ARTR"){
    rootratefit.full<-glm(germ_rate~Day_Scale_7+Night_Scale_7+Treatment,
                              family=Gamma(link="log"),data=temp_agg3_germ)
    print(summary(rootratefit.full)) 
  }
  
  if (s=="ACMI"){
    rootratefit.full<-glm(germ_rate~Day_Scale_7+Night_Scale_7+
                                Day_Scale_7:Night_Scale_7,
                              family=Gamma(link="log"),data=temp_agg3_germ)
    print(summary(rootratefit.full)) 
  }
  
  if (s=="ELEL"){
    rootratefit.full<-glm(germ_rate~Day_Scale_7+Night_Scale_7+Treatment+
                                Day_Scale_7:Treatment+Day_Scale_7:Night_Scale_7,
                              family=Gamma(link="log"),data=temp_agg3_germ)
    print(summary(rootratefit.full)) 
  }
  
  root.rate.full.list[[s]] <- rootratefit.full

}

##emerg

shoot.rate.full.list <- list()

for (s in species) {
  
  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  temp_agg3<-temp_agg2[-which(temp_agg2$emerg_rate==0),]
  emerg_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","emerg_rate") #select the columns I want for the germ model
  temp_agg3_emerg<-temp_agg2[,emerg_cols] # make a new dataframe with just those values
  
  shootratefit.full<-glm(emerg_rate~Day_Scale_7+Night_Scale_7,
                        family=Gamma(link="log"),data=temp_agg3_emerg)
  print(summary(shootratefit.full))
  
  if (s=="POSE"){
    shootratefit.full<-glm(emerg_rate~Day_Scale_7+Night_Scale_7+Treatment+
                            Day_Scale_7:Treatment+Night_Scale_7:Treatment,
                          family=Gamma(link="log"),data=temp_agg3_emerg)
    print(summary(shootratefit.full)) 
  }
  
  if (s=="ARTR"){
    shootratefit.full<-glm(emerg_rate~Day_Scale_7+Night_Scale_7+Treatment,
                          family=Gamma(link="log"),data=temp_agg3_emerg)
    print(summary(shootratefit.full)) 
  }
  
  if (s=="ELEL"){
    shootratefit.full<-glm(emerg_rate~Day_Scale_7+Night_Scale_7+Treatment+
                            Day_Scale_7:Treatment+Day_Scale_7:Night_Scale_7,
                          family=Gamma(link="log"),data=temp_agg3_emerg)
    print(summary(shootratefit.full)) 
  }
  
  shoot.rate.full.list[[s]] <- shootratefit.full
  
}


