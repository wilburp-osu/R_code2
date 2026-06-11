# Create Open List For Clean Data
clean.data.list<-list()

# Run ForLoop for Growth Rates after Seed Germinates

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
# Open List For Aggregation
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
  
  temp_agg$Treatment<-as.factor(paste0(temp_agg$Treatment))
  agg.list[[s]] <- temp_agg
}

## aggregate data for each cuvette for each day for only seeds that germinated####
agg.list.sub.EGrate<- list()

for (s in species) {
  temp <-clean.data.list[[s]]
  ## is it right to do this with Shoot rate? its the only way to get it to work
  temp2<-temp[-which(temp$Avg_Shoot_rate==0),]
  temp_agg<-aggregate(cbind(temp2$Avg_Root_rate,
                            temp2$Avg_Shoot_rate),
                      by=list(temp2$Day_Scale_7,temp2$Night_Scale_7,
                              temp2$Treatment,
                              temp2$Group),mean,na.rm=T)
  
  temp_agg[temp_agg =="NaN"]<-NA
  names(temp_agg)<-c("Day_Scale_7","Night_Scale_7","Treatment",
                     "Group","Avg_Root_rate","Avg_Shoot_rate")
  
  temp_agg$Treatment<-as.factor(paste0(temp_agg$Treatment))
  agg.list.sub.EGrate[[s]] <- temp_agg
}

# Run Step Function Forward With Backwards Look To Get Final Models


for (s in species){
  
  print(s)
  temp_agg2 <- agg.list.sub.EGrate[[s]] # pull out the data
  sgrowth_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","Avg_Shoot_rate") #select the columns I want for the germ model
  temp_agg2_sgrowth<-temp_agg2[,sgrowth_cols] # make a new dataframe with just those values
  
  #fit the null model
  null.model<-glm(Avg_Shoot_rate~. ,
                  family=Gamma(link="log"),data=temp_agg2_sgrowth)
  
  #fit the global models
  global.model<-glm(Avg_Shoot_rate~Day_Scale_7+Night_Scale_7+
                      Day_Scale_7:Treatment+Night_Scale_7:Treatment+
                      Treatment+Day_Scale_7:Night_Scale_7+
                      Day_Scale_7:Night_Scale_7:Treatment,
                    family=Gamma(link="log"),data=temp_agg2_sgrowth)
  
  #run the step function in both directions bounded by the null and global model
  models.step<- step(null.model,direction="both",scope = list(lower=null.model, upper=global.model))
  
  print(summary(models.step))
}

shoot.rate.full.list <- list()
shoot.rate.output <- list()

## Fitting and Getting Model Output for Each Species ####

for (s in species) {
  
  print(s)
  temp_agg2 <- agg.list.sub.EGrate[[s]] # pull out the data
  temp_agg3<-temp_agg2[-which(temp_agg2$Avg_Shoot_rate==0),]
  sgrowth_cols<-c("Treatment","Day_Scale_7","Night_Scale_7","Avg_Shoot_rate") #select the columns I want for the germ model
  temp_agg3_sgrowth<-temp_agg2[,sgrowth_cols] # make a new dataframe with just those values
  
  # create glm of root growth rate
  shootratefit.full<-glm(Avg_Shoot_rate~Day_Scale_7+Night_Scale_7+Treatment,
                        family=Gamma(link="log"),data=temp_agg3_sgrowth)
  print(summary(shootratefit.full))
  
  if (s=="POSE"){
    shootratefit.full<-glm(Avg_Shoot_rate~Day_Scale_7+Night_Scale_7+Treatment+
                            Day_Scale_7:Night_Scale_7,
                          family=Gamma(link="log"),data=temp_agg3_sgrowth)
    print(summary(shootratefit.full)) 
  }
  
  if (s=="ACMI"){
    shootratefit.full<-glm(Avg_Shoot_rate~Day_Scale_7+Night_Scale_7+Treatment+
                            Day_Scale_7:Night_Scale_7,
                          family=Gamma(link="log"),data=temp_agg3_sgrowth)
    print(summary(shootratefit.full)) 
  }
  
  if (s=="ELEL"){
    shootratefit.full<-glm(Avg_Shoot_rate~Day_Scale_7+Night_Scale_7+Treatment+
                            Day_Scale_7:Night_Scale_7,
                          family=Gamma(link="log"),data=temp_agg3_sgrowth)
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

## Create Raster Figures of Containing Each Species

output.all.shootrate <- do.call(rbind, shoot.rate.output)

# Input dummy coefficients

output.all.shootrate[28,] <- c(NA,NA,NA,100,"POSE","Day_Scale_7:TreatmentW")
output.all.shootrate[29,] <- c(NA,NA,NA,100,"POSE","Day_Scale_7:TreatmentL")
output.all.shootrate[30,] <- c(NA,NA,NA,100,"POSE","Day_Scale_7:TreatmentH")

output.all.shootrate$coefficients<-factor(output.all.shootrate$coefficients,
                                         levels = c("Day_Scale_7:TreatmentW",
                                                    "Day_Scale_7:TreatmentL",
                                                    "Day_Scale_7:TreatmentH",
                                                    "Day_Scale_7:Night_Scale_7",
                                                    "(Intercept)", "TreatmentW",
                                                    "TreatmentL","TreatmentH",
                                                    "Day_Scale_7","Night_Scale_7"),
                                         labels = c("Day Temperature: \n Water Treatment", 
                                                    "Day Temperature: \n Low SA Treatment", 
                                                    "Day Temperature: \n High SA Treatment", 
                                                    "Day Temperature: \n Night Temperature",
                                                    "(Intercept)", "Water Treatment", 
                                                    "Low SA Treatment", "High SA Treatment", 
                                                    "Day Temperature", "Night Temperature"))

output.all.shootrate$sigvar <-
  ifelse(as.numeric(output.all.shootrate$Pr...t..) < 0.05, 1, 0)

shootrate_sig<-output.all.shootrate[which(output.all.shootrate$sigvar==1),]         

Rounding<- function (x) {ifelse (is.numeric(x),round(x,4),x)}

output.all.shootrate$Estimate <- as.numeric(output.all.shootrate$Estimate)

Shoot_Growth_Raster <- ggplot(output.all.shootrate[-which(output.all.shootrate$coefficients=="(Intercept)"),],
       aes(x=Species,y=coefficients))+
  geom_raster(aes(fill=as.numeric(Estimate)))+
  geom_tile(data =shootrate_sig[-which(shootrate_sig$coefficients=="(Intercept)"),], 
            aes(x=Species,y=coefficients, fill = sigvar),fill="transparent",
            col="black",size=2)+
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", na.value = "transparent",
                       midpoint = 0, space = "Lab",
                       name="log(Effect) (mm/day)")+
  geom_text(aes(label = (label = ifelse(is.na(Estimate), "", paste0(round((Estimate),3)))), lsize = 6))+
  ylab("")+
  xlab("")+
  theme_grey()+
  theme(legend.position="bottom", legend.key.width = unit(1, 'cm'))

print(Shoot_Growth_Raster)

ggsave("Figures/shootrate_sig.png", Shoot_Growth_Raster)
