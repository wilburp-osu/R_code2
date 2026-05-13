# Create Open List For Clean Data
clean.data.list<-list()
# Run a forloop to clean data and calculate emergence rate & days until emergence
for (s in species){
  temp<-data.list[[s]]
  
  ## clean data ####
  temp$Root_Pheno_Code<-as.numeric(paste0(temp$Root_Pheno_Code))
  temp$Shoot_Pheno_Code<-as.numeric(paste0(temp$Shoot_Pheno_Code))
  # Create Empty Columns To Fill
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

# Open list for aggregated data
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
  
  # calculate means of emergence and germination by individual seed
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

# Create second open list 
# This one labelled so that it can be used in boxplots as well
short.data.list2.Emerg <- list()

for (s in species){
  temp3<-clean.data.list.short[[s]]
  
  # Changing Zeroes to NA for seeds that never germinated or emerged
  
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
  
  short.data.list2.Emerg[[s]]<-temp3
  
}

# Run Step Function Forward With Backwards Look To Get Final Models

for (s in species){
  
  print(s)
  temp_agg2 <- short.data.list2.Emerg[[s]] # pull out the data
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
  
}

shootDaysfit.full.list <- list()
shootDays.output <- list()

## Fitting and Getting Model Output for Each Species ####

for (s in species){
  print(s)
  temp_agg2 <- short.data.list2.Emerg[[s]] # pull out the data
  emerg_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","days_til_emerg") #select the columns I want for the germ model
  temp_agg2_emerg<-temp_agg2[,emerg_cols] # make a new dataframe with just those values
  
  shootDaysfull<-glm(days_til_emerg~Treatment+Day_Scale_7+Night_Scale_7,
                     family=poisson(link = "log"),data=temp_agg2_emerg)
  print(summary(shootDaysfull))
  
  shootDaysfit.full.list[[s]] <- shootDaysfull
  
  shootDaysrast <- data.frame(summary(shootDaysfull)$coefficients)
  shootDaysrast$Species <- s
  shootDaysrast$coefficients <- rownames(shootDaysrast)
  shootDays.output[[s]]<-shootDaysrast
  
}

## Create Raster Figures of Containing Each Species

output.all.shootdays <- do.call(rbind, shootDays.output)

# add in dummy coefficients 

output.all.shootdays[25,] <- c(NA,NA,NA,100,"POSE","Day_Scale_7:TreatmentW")
output.all.shootdays[26,] <- c(NA,NA,NA,100,"POSE","Day_Scale_7:TreatmentL")
output.all.shootdays[27,] <- c(NA,NA,NA,100,"POSE","Day_Scale_7:TreatmentH")
output.all.shootdays[28,] <- c(NA,NA,NA,100,"POSE","Day_Scale_7:Night_Scale_7")

output.all.shootdays$coefficients<-factor(output.all.shootdays$coefficients,
                                          levels = c("Day_Scale_7:TreatmentW",
                                                     "Day_Scale_7:TreatmentL", 
                                                     "Day_Scale_7:TreatmentH",
                                                     "Day_Scale_7:Night_Scale_7",
                                                     "(Intercept)", "TreatmentW",
                                                     "TreatmentL","TreatmentH",
                                                     "Day_Scale_7", "Night_Scale_7"),
                                          labels = c("Day Temperature: \n Water Treatment", 
                                                     "Day Temperature: \n Low SA Treatment", 
                                                     "Day Temperature: \n High SA Treatment", 
                                                     "Day Temperature: \n Night Temperature",
                                                     "(Intercept)", "Water Treatment", 
                                                     "Low SA Treatment", "High SA Treatment", 
                                                     "Day Temperature", "Night Temperature"))


output.all.shootdays$sigvar <-
  ifelse(as.numeric(output.all.shootdays$Pr...z..) < 0.05, 1, 0)


shootDays_sig<-output.all.shootdays[which(output.all.shootdays$sigvar==1),]         

Rounding<- function (x) {ifelse (is.numeric(x),round(x,4),x)} #do we need this any more?

output.all.shootdays$Estimate <- as.numeric(output.all.shootdays$Estimate)

Emergence_Rasters <- ggplot(output.all.shootdays[-which(output.all.shootdays$coefficients=="(Intercept)"),],
       aes(x=Species,y=coefficients))+
  geom_raster(aes(fill=as.numeric(Estimate)))+
  geom_tile(data =shootDays_sig[-which(shootDays_sig$coefficients=="(Intercept)"),], 
            aes(x=Species,y=coefficients, fill = sigvar),fill="transparent",
            col="black",linewidth=2)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",na.value = "transparent",
                       midpoint = 0, space = "Lab",
                       name="log(Effect)")+
  #geom_text(aes(label = (label = ifelse(is.na(Estimate), "", paste0(round(as.numeric((exp(Estimate)*sign(Estimate))),4))))))+
  ylab("")+
  xlab("")+
  theme_grey()+
  theme(legend.position="bottom") # using base_size = 22 will increase text size

ggsave("Figures/shootDays_sig.png", Emergence_Rasters)

