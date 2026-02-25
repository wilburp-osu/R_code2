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
#setwd("C:/Users/18034/Dropbox/PC/Desktop/R_code/R_code2")

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

##STEPS
#Germ

### for stepwise regression you usually start with a null model (for forward selection at least) and
##   include a dataset with just the columns that you potentially want to include
##   in the model

for (s in species){
  
  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  germ_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","days_til_germ") #select the columns I want for the germ model
  temp_agg2_germ<-temp_agg2[,germ_cols] # make a new dataframe with just those values

  
    
  #fit the null model
  null.model<-glm(days_til_germ~. ,
                  family=poisson(link = "log"),data=temp_agg2_germ)
  
  #fit the global models
  global.model<-glm(days_til_germ~Day_Scale_7+Night_Scale_7+
                      Day_Scale_7:Treatment+Night_Scale_7:Treatment+
                      Treatment+Day_Scale_7:Night_Scale_7+
                      Day_Scale_7:Night_Scale_7:Treatment,
                    family=poisson(link = "log"),data=temp_agg2_germ)
  
  #run the step function in both directions bounded by the null and global model
  models.step<- step(null.model,direction="both",scope = list(lower=null.model, upper=global.model))
  
  print(summary(models.step))
  #models.step$anova
  
}

## It seems like we are getting similar results with dredge and with step. In the end, 
#    we will have to choose one and go forward from there for the manuscript.

#Emerg

for (s in species){
  
  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  emerg_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","days_til_emerg") #select the columns I want for the germ model
  temp_agg2_emerg<-temp_agg2[,emerg_cols] # make a new dataframe with just those values
  
  
  
  #fit the null model
  null.model<-glm(days_til_emerg~. ,
                  family=poisson(link = "log"),data=temp_agg2_emerg)
  
  #fit the global models
  global.model<-glm(days_til_emerg~Day_Scale_7+Night_Scale_7+
                      Day_Scale_7:Treatment+Night_Scale_7:Treatment+
                      Treatment+Day_Scale_7:Night_Scale_7+
                      Day_Scale_7:Night_Scale_7:Treatment,
                    family=poisson(link = "log"),data=temp_agg2_emerg)
  
  #run the step function in both directions bounded by the null and global model
  models.step<- step(null.model,direction="both",scope = list(lower=null.model, upper=global.model))
  
  print(summary(models.step))
  #models.step$anova
  
}

## Rates


#Germ
for (s in species){
  
  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  grate_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","germ_rate") #select the columns I want for the germ model
  temp_agg2_grate<-temp_agg2[,grate_cols] # make a new dataframe with just those values
  
  
  
  #fit the null model
  null.model<-glm(germ_rate~. ,
                  family=Gamma(link = "log"),data=temp_agg2_grate)
  
  #fit the global models
  global.model<-glm(germ_rate~Day_Scale_7+Night_Scale_7+
                      Day_Scale_7:Treatment+Night_Scale_7:Treatment+
                      Treatment+Day_Scale_7:Night_Scale_7+
                      Day_Scale_7:Night_Scale_7:Treatment,
                    family=Gamma(link = "log"),data=temp_agg2_grate)
  
  #run the step function in both directions bounded by the null and global model
  models.step<- step(null.model,direction="both",scope = list(lower=null.model, upper=global.model))
  
  print(summary(models.step))
  #models.step$anova
  
}


#Emerg
for (s in species){
  
  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  erate_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","emerg_rate") #select the columns I want for the germ model
  temp_agg2_erate<-temp_agg2[,erate_cols] # make a new dataframe with just those values
  
  
  
  #fit the null model
  null.model<-glm(emerg_rate~. ,
                  family=Gamma(link = "log"),data=temp_agg2_erate)
  
  #fit the global models
  global.model<-glm(emerg_rate~Day_Scale_7+Night_Scale_7+
                      Day_Scale_7:Treatment+Night_Scale_7:Treatment+
                      Treatment+Day_Scale_7:Night_Scale_7+
                      Day_Scale_7:Night_Scale_7:Treatment,
                    family=Gamma(link = "log"),data=temp_agg2_erate)
  
  #run the step function in both directions bounded by the null and global model
  models.step<- step(null.model,direction="both",scope = list(lower=null.model, upper=global.model))
  
  print(summary(models.step))
  #models.step$anova
  
}

### With BIC
## days til
#Germ

for (s in species){
  
  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  germ_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","days_til_germ") #select the columns I want for the germ model
  temp_agg2_germ<-temp_agg2[,germ_cols] # make a new dataframe with just those values
  
  
  
  #fit the null model
  null.model<-glm(days_til_germ~. ,
                  family=poisson(link = "log"),data=temp_agg2_germ)
  
  #fit the global models
  global.model<-glm(days_til_germ~Day_Scale_7+Night_Scale_7+
                      Day_Scale_7:Treatment+Night_Scale_7:Treatment+
                      Treatment+Day_Scale_7:Night_Scale_7+
                      Day_Scale_7:Night_Scale_7:Treatment,
                    family=poisson(link = "log"),data=temp_agg2_germ)
  
  #run the step function in both directions bounded by the null and global model
  models.step<- step(null.model,direction="both",
                     scope = list(lower=null.model, upper=global.model),
                     k = log(784)) ##784 is the number of observations
  ## I think this makes what is labelled as AIC turn to BIC
  
  print(summary(models.step))
  #models.step$anova
  
}

##Emerg
for (s in species){
  
  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  emerg_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","days_til_emerg") #select the columns I want for the germ model
  temp_agg2_emerg<-temp_agg2[,emerg_cols] # make a new dataframe with just those values
  
  
  
  #fit the null model
  null.model<-glm(days_til_emerg~. ,
                  family=poisson(link = "log"),data=temp_agg2_emerg)
  
  #fit the global models
  global.model<-glm(days_til_emerg~Day_Scale_7+Night_Scale_7+
                      Day_Scale_7:Treatment+Night_Scale_7:Treatment+
                      Treatment+Day_Scale_7:Night_Scale_7+
                      Day_Scale_7:Night_Scale_7:Treatment,
                    family=poisson(link = "log"),data=temp_agg2_emerg)
  
  #run the step function in both directions bounded by the null and global model
  models.step<- step(null.model,direction="both",
                     scope = list(lower=null.model, upper=global.model),
                     k = log(784)) ##784 is the number of observations
  ## I think this makes what is labelled as AIC turn to BIC
  
  print(summary(models.step))
  #models.step$anova
  
}

##Germ rate
for (s in species){
  
  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  grate_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","germ_rate") #select the columns I want for the germ model
  temp_agg2_grate<-temp_agg2[,grate_cols] # make a new dataframe with just those values
  
  
  
  #fit the null model
  null.model<-glm(germ_rate~. ,
                  family=Gamma(link = "log"),data=temp_agg2_grate)
  
  #fit the global models
  global.model<-glm(germ_rate~Day_Scale_7+Night_Scale_7+
                      Day_Scale_7:Treatment+Night_Scale_7:Treatment+
                      Treatment+Day_Scale_7:Night_Scale_7+
                      Day_Scale_7:Night_Scale_7:Treatment,
                    family=Gamma(link = "log"),data=temp_agg2_grate)
  
  #run the step function in both directions bounded by the null and global model
  models.step<- step(null.model,direction="both",
                     scope = list(lower=null.model, upper=global.model),
                     k = log(784)) ##784 is the number of observations
  ## I think this makes what is labelled as AIC turn to BIC
  
  print(summary(models.step))
  #models.step$anova
  
}

##Emerg rate
for (s in species){
  
  print(s)
  temp_agg2 <- short.data.list2[[s]] # pull out the data
  erate_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","emerg_rate") #select the columns I want for the germ model
  temp_agg2_erate<-temp_agg2[,erate_cols] # make a new dataframe with just those values
  
  
  
  #fit the null model
  null.model<-glm(emerg_rate~. ,
                  family=Gamma(link = "log"),data=temp_agg2_erate)
  
  #fit the global models
  global.model<-glm(emerg_rate~Day_Scale_7+Night_Scale_7+
                      Day_Scale_7:Treatment+Night_Scale_7:Treatment+
                      Treatment+Day_Scale_7:Night_Scale_7+
                      Day_Scale_7:Night_Scale_7:Treatment,
                    family=Gamma(link = "log"),data=temp_agg2_erate)
  
  #run the step function in both directions bounded by the null and global model
  models.step<- step(null.model,direction="both",
                     scope = list(lower=null.model, upper=global.model),
                     k = log(784)) ##784 is the number of observations
  ## I think this makes what is labelled as AIC turn to BIC
  
  print(summary(models.step))
  #models.step$anova
  
}

