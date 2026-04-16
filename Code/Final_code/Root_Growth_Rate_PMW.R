######################################
clean.data.list<-list()

# GROWTH RATE Calculations #####

for (s in species) {
  temp<-data.list[[s]]
  ### clean data 
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

## aggregate data for each cuvette for each day ####
#str(clean.data.list)
agg.list<- list()

for (s in species) {
  temp <-clean.data.list[[s]]
  temp_agg<-aggregate(cbind(temp$Avg_Root_rate,
                            temp$Avg_Shoot_rate),
                      by=list(temp$Day_Scale_7,temp$Night_Scale_7,
                              temp$Treatment,
                              temp$Group),mean,na.rm=T)
  
  temp_agg[temp_agg =="NaN"]<-NA
  names(temp_agg)<-c("Day_Scale_7","Night_Scale_7","Treatment",
                     "Group","Avg_Root_rate","Avg_Shoot_rate")
  
  
  #temp_agg$Temp_Ratio<-temp_agg$Day_Temp/temp_agg$Night_Temp
  temp_agg$Treatment<-as.factor(paste0(temp_agg$Treatment))
  #temp_agg$Day_Scale_7<-scale(temp_agg$Day_Temp)
  #temp_agg$Night_Scale_7<-scale(temp_agg$Night_Temp)
  agg.list[[s]] <- temp_agg
}

## aggregate data for each cuvette for each day for only seeds that germinated####
#str(clean.data.list)
agg.list.sub<- list()

for (s in species) {
  temp <-clean.data.list[[s]]
  temp2<-temp[-which(temp$Avg_Root_rate==0),]
  temp_agg<-aggregate(cbind(temp2$Avg_Root_rate,
                            temp2$Avg_Shoot_rate),
                      by=list(temp2$Day_Scale_7,temp2$Night_Scale_7,
                              temp2$Treatment,
                              temp2$Group),mean,na.rm=T)
  
  temp_agg[temp_agg =="NaN"]<-NA
  names(temp_agg)<-c("Day_Scale_7","Night_Scale_7","Treatment",
                     "Group","Avg_Root_rate","Avg_Shoot_rate")
  
  
  #temp_agg$Temp_Ratio<-temp_agg$Day_Temp/temp_agg$Night_Temp
  temp_agg$Treatment<-as.factor(paste0(temp_agg$Treatment))
  #temp_agg$Day_Scale_7<-scale(temp_agg$Day_Temp)
  #temp_agg$Night_Scale_7<-scale(temp_agg$Night_Temp)
  agg.list.sub[[s]] <- temp_agg
}

#Germ

# for (s in species){
#   print(s)
#   temp_agg2 <- agg.list.sub[[s]] 
#   
#   global.model<-glm(Avg_Root_rate~Day_Scale_7+Night_Scale_7+
#                       Day_Scale_7:Treatment+Night_Scale_7:Treatment+
#                       Treatment+Day_Scale_7:Night_Scale_7+
#                       Day_Scale_7:Night_Scale_7:Treatment,
#                     family=Gamma(link = "log"),data=temp_agg2)
#   
#   models.step<- step(global.model)
#   
# }

for (s in species){
  
  print(s)
  temp_agg2 <- agg.list.sub[[s]] # pull out the data
  rgrowth_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","Avg_Root_rate") #select the columns I want for the germ model
  temp_agg2_rgrowth<-temp_agg2[,rgrowth_cols] # make a new dataframe with just those values
  
  #fit the null model
  null.model<-glm(Avg_Root_rate~. ,
                  #should link be inverse?
                  family=Gamma(link = "log"),data=temp_agg2_rgrowth)
  
  #fit the global models
  global.model<-glm(Avg_Root_rate~Day_Scale_7+Night_Scale_7+
                      Day_Scale_7:Treatment+Night_Scale_7:Treatment+
                      Treatment+Day_Scale_7:Night_Scale_7+
                      Day_Scale_7:Night_Scale_7:Treatment,
                    family=Gamma(link = "log"),data=temp_agg2_rgrowth)
  
  #run the step function in both directions bounded by the null and global model
  models.step<- step(null.model,direction="both",scope = list(lower=null.model, upper=global.model))
  
  print(summary(models.step))
  }
  
#models.step$anova
  
root.rate.full.list <- list()
root.rate.output <- list()

for (s in species) {

  print(s)
  temp_agg2 <- agg.list.sub[[s]] # pull out the data
  temp_agg3<-temp_agg2[-which(temp_agg2$Avg_Root_rate==0),]
  rgrowth_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","Avg_Root_rate") #select the columns I want for the germ model
  temp_agg3_rgrowth<-temp_agg2[,rgrowth_cols] # make a new dataframe with just those values
  
  rootratefit.full<-glm(Avg_Root_rate~Day_Scale_7+Night_Scale_7+Treatment,
                            family=Gamma(link="log"),data=temp_agg3_rgrowth)
  print(summary(rootratefit.full))
  
  if (s=="POSE"){
    rootratefit.full<-glm(Avg_Root_rate~Day_Scale_7+Night_Scale_7+Treatment+
                                Day_Scale_7:Night_Scale_7,
                              family=Gamma(link="log"),data=temp_agg3_rgrowth)
    print(summary(rootratefit.full)) 
  }
  
  if (s=="ARTR"){
    rootratefit.full<-glm(Avg_Root_rate~Day_Scale_7+Night_Scale_7+Treatment+
                                Treatment:Day_Scale_7,
                              family=Gamma(link="log"), data=temp_agg3_rgrowth)
    print(summary(rootratefit.full))
  }
  
  if (s=="ACMI"){
    rootratefit.full<-glm(Avg_Root_rate~Day_Scale_7+Night_Scale_7+Treatment+
                                Day_Scale_7:Night_Scale_7,
                              family=Gamma(link="log"),data=temp_agg3_rgrowth)
    print(summary(rootratefit.full)) 
  }
  
  if (s=="ELEL"){
    rootratefit.full<-glm(Avg_Root_rate~Day_Scale_7+Night_Scale_7+Treatment+
                                Day_Scale_7:Night_Scale_7,
                              family=Gamma(link="log"),data=temp_agg3_rgrowth)
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

output.all.rootrate$coefficients<-factor(output.all.rootrate$coefficients,
                                         levels = c("Day_Scale_7:TreatmentW",
                                                    "Day_Scale_7:TreatmentL",
                                                    "Day_Scale_7:TreatmentH",
                                                    "Day_Scale_7:Night_Scale_7",
                                                    "(Intercept)", "TreatmentW",
                                                    "TreatmentL","TreatmentH",
                                                    "Day_Scale_7","Night_Scale_7"),
                                         labels = c("Day Temperature:Water Treatment",
                                                    "Day Temperature: Low SA Treatment", 
                                                    "Day Temperature: High SA Treatment",
                                                    "Day Temperature: Night Temperature",
                                                    "(Intercept)", "Water Treatment", 
                                                    "Low SA Treatment", "High SA Treatment", 
                                                    "Day Temperature", "Night Temperature"))

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
                       #log(Effect) because of the Gamma distribution
                       name="log(Effect)")+
  geom_text(aes(label = paste0(round(Estimate,4))))+
  ylab("")+
  xlab("")+
  #ylab("Variables & Interactions")+
  #xlab("Species")+
  #ggtitle("Temperature & Treatments Interactions with Root Growth Rate")+
  theme_grey()+
  theme(plot.title = element_text(margin = margin(l = -20)))


