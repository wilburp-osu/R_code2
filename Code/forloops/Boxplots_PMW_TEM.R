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
library(ggplot2)

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


short.all<-do.call(rbind,short.data.list2)  

##PLOTS

short.data.list2[["POSE"]]$Temp_bin <- cut(short.data.list2[["POSE"]]$Day_Scale_7, breaks = 7)

ggplot(short.data.list2[["POSE"]], aes(x = Temp_bin, y = days_til_germ, 
                                       fill = Night_Scale_7)) +
  geom_boxplot() +
  geom_line(aes(mean(days_til_germ)))



# Create a scatter plot and facet by 'cut'
ggplot(short.data.list2[["POSE"]], aes(x = Night_Scale_7, y = days_til_germ)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ Day_Scale_7, nrow= 1) + # This creates a separate plot for each 'cut' level
  labs(title = "Days Til germ vs Night temp faceted by Day temp",
       x = "Night Temp",
       y = "Days til") #+
theme_minimal()

## base r

for  (s in species){
  
  boxplot(short.data.list2[[s]]$days_til_germ~short.data.list2[[s]]$Night_Scale_7,
          col="LightBlue",
          main=paste0("germ_night (", s, ")"),
          ylab="Days Until Germination",
          xlab="Night Temperature")
  
  global.model<-lm(short.data.list2[[s]]$days_til_germ~
                     short.data.list2[[s]]$Night_Scale_7)
  
  abline(global.model,lwd=2)
  
}

## ggplot

#night scale

# 
#   night.t.boxplots <- ggplot(short.all, aes(x=factor(Night_Scale_7), y=days_til_germ,fill=Treatment,col=Treatment)) +
#     geom_boxplot() +
#     stat_boxplot(geom ='errorbar') +
#     geom_smooth(se=FALSE, color="black", aes(group=1))+
#     stat_summary(fun = "mean", geom = "point", size = 3, shape = 15)+
#     labs( title = paste0("Days Until Germ x Night Scale (", s, ")"),
#           x = "Night Temperature",
#           y = "Days Until Germination")


### TEM code with just a linear model #########  
## I could see this being useful for just a general linear trend
## but we cannot see any variation or change in the trend

for (s in species){
  
  temp_agg3 <- short.data.list2[[s]] 
  
  night.t.boxplots <- ggplot(temp_agg3, aes(factor(round(Night_Scale_7,1)), days_til_germ,fill=Treatment)) +
    geom_smooth(method="lm",se=F,aes(group=Treatment,col=Treatment),size=2)+ ##put the line int he background
    geom_boxplot() +
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0("Days Until Germ x Night Scale (", s, ")"),
          x = "Night Temperature",
          y = "Days Until Germination")+
    theme_minimal(base_size = 15)+
    scale_color_manual(values=c("grey60","green3","orchid3","dodgerblue3"))+
    scale_fill_manual(values=c("grey80","green","orchid1","skyblue"))
  
  print(night.t.boxplots)
  
}



### TEM code with just a gam model varying the flexibility (span) #########  

## I could see this being useful for overall 
## trends (especially where ELEL inflects back up in the end)
## but could be too busy


for (s in species){
  
  temp_agg3 <- short.data.list2[[s]] 
  
  night.t.boxplots <- ggplot(temp_agg3, aes(factor(round(Night_Scale_7,1)), days_til_germ,fill=Treatment)) +
    geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,size=2)+ ##put the line int he background
    #geom_smooth(se=F,aes(group=Treatment,col="black"),span = 0.9,size=1)+ ##put the line int he background
    geom_boxplot() +
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0("Days Until Germ x Night Scale (", s, ")"),
          x = "Night Temperature",
          y = "Days Until Germination")+
    theme_minimal(base_size = 15)+
    scale_color_manual(values=c("grey60","green3","orchid3","dodgerblue3"))+
    scale_fill_manual(values=c("grey80","green","orchid1","skyblue"))
  
  print(night.t.boxplots)
  
}  

### TEM edit for day scale - I used the GAM here
#day scale

for (s in species){
  
  temp_agg3 <- short.data.list2[[s]] 
  
  day.t.boxplots <- ggplot(temp_agg3, aes(factor(round(Day_Scale_7,2)), days_til_germ,fill=Treatment)) +
    geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,size=2)+ ##put the line int he background
    #geom_smooth(se=F,aes(group=Treatment,col="black"),span = 0.9,size=1)+ ##put the line int he background
    geom_boxplot() +
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0("Days Until Germ x Day Scale (", s, ")"),
          x = "Day Temperature",
          y = "Days Until Germination")+
    theme_minimal(base_size = 15)+
    scale_color_manual(values=c("grey60","green3","orchid3","dodgerblue3"))+
    scale_fill_manual(values=c("grey80","green","orchid1","skyblue"))
  
  print(day.t.boxplots)
  
}

#treatment

for (s in species){
  
  temp_agg3 <- short.data.list2[[s]] 
  
  treatment.boxplots <- ggplot(temp_agg3, aes(x = Treatment, y = days_til_germ)) +
    geom_boxplot() +
    stat_boxplot(geom ='errorbar') +
    geom_smooth(method = "glm", aes(y= days_til_germ, x = Treatment, group=1),
                formula = y~x, se=FALSE, color="black")+
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 15)+
    labs( title = paste0("Days Until Germ x Treatment (", s, ")"),
          x = "Treatment",
          y = "Days Until Germination")
  
  print(treatment.boxplots)
}


#test

for (s in species){
  
  temp_agg3 <- short.data.list2[[s]] 
  
  germ.rate.boxplots <- ggplot(temp_agg3, aes(factor(Night_Scale_7:Day_Scale_7), germ_rate)) +
    geom_boxplot() +
    stat_boxplot(geom ='errorbar') +
    geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))+
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 15)+
    labs( title = paste0("Germ Rate x Night Temp (", s, ")"),
          # x = "Night Temp",
          y = "Germination Rate")
  
  print(germ.rate.boxplots)
  
}

#next steps
# am I doing this right? 
# Figure out multiple variables (i.e. day:night vs rate)
# is that what I am supposed to do?

## Test

for (s in species){
  
  temp_agg3 <- short.data.list2[[s]] 
  
  germ.rate.boxplots <- ggplot(temp_agg3, aes(x = factor(Night_Scale_7), 
                                              y = germ_rate)) +
    geom_boxplot() +
    stat_boxplot(geom ='errorbar') +
    #geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))+
    geom_smooth(method = "glm", aes(x = factor(Night_Scale_7), y = germ_rate, group = 1),
                formula = y ~ x , se= FALSE, color="black")+
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 15)+
    labs( title = paste0("Germ Rate x Night Temp (", s, ")"),
          x = "Night Temp",
          y = "Germination Rate")
  
  print(germ.rate.boxplots)
}


for (s in species){
  
  temp_agg3 <- short.data.list2[[s]] 
  
  treatment.boxplots <- ggplot(temp_agg3, aes(x = Treatment, y = days_til_germ)) +
    geom_boxplot() +
    stat_boxplot(geom ='errorbar') +
    geom_smooth(method = "glm", aes(y= days_til_germ, x = Treatment, group=1),
                formula = y~x, se=FALSE, color="black")+
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 15)+
    labs( title = paste0("Days Until Germ x Treatment (", s, ")"),
          x = "Treatment",
          y = "Days Until Germination")
  
  print(treatment.boxplots)
}


model_shift <- glm(germ_rate~Night_Scale_7,
                   #use gamma for non integer
                   family=Gamma(link = "log"),data=temp_agg3)
x <- seq(min(temp_agg3$Night_Scale_7), max(temp_agg3$Night_Scale_7), .1)
y <- predict(model_shift, list(Night_Scale_7 = x), type = 'response')
plot_data <- data.frame(Night_Scale_7 = x, germ_rate = y)

ggplot(aes(x = Night_Scale_7, y = germ_rate), data = plot_data) +
  geom_point()
