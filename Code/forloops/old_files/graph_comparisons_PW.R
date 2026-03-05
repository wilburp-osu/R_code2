install.packages("patchwork")
library(patchwork)

####################################
###Arranging plots in 2x2 GRIDs

###########
#Shoot Growth Rate

plots <- list()

for (s in species){
  print(s)
  pnew.data.sh <- predict.list.sh[[s]]
  temp_agg2 <- agg.list.sub[[s]]
  
  temp.shoot.rates.C<-temp_agg2$Avg_Shoot_rate[which(temp_agg2$Treatment=="C")]
  #temp.shoot.rates<-temp_agg2$Avg_Shoot_rate
  #summary(temp.shoot.rates.C)
  #summary(temp.shoot.rates)
  
  mid.pointss<-mean(temp.shoot.rates.C,na.rm=T)
  
  print(mid.pointss)
  
  plot.sh <- ggplot(pnew.data.sh,aes(x=Day_Scale_7,y=Night_Scale_7))+
    geom_raster(aes(fill=Prediction_back))+
    facet_grid(.~Treatment)+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = mid.pointss, space = "Lab",
                         name="Shoot growth rate (final root length/days of growth)") +
    theme_minimal()+
    ylab("Night Temperature (C)")+
    xlab("Day Temperature (C)")+
    theme(legend.position = "bottom")+
    labs(title = paste0(s))+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  
  CairoPNG(file.path(paste0("Figures/",s,"_Shoot_growth_rate.png")), width = 12, height = 7, units = "in",
           res = 300)
  print(plot.sh)
  dev.off()
  plots[[s]] <- plot.sh
  
}

# Arrange in 2x2 grid
(plots[["POSE"]] | plots[["ACMI"]]) /
  (plots[["ARTR"]] | plots[["ELEL"]])

#############
##Root Growth Rate

plots1 <- list()

for (s in species){
  print(s)
  pnew.data <- predict.list[[s]]
  temp_agg2 <- agg.list.sub[[s]]
  
  #temp_agg2<-temp_agg[-which(temp_agg$Avg_Root_rate==0),]
  mid.pointss<-mean(temp_agg2$Avg_Root_rate[which(temp_agg2$Treatment=="C")],na.rm=T)
  
  print(mid.pointss)
  
  plots <- ggplot(pnew.data,aes(x=Day_Scale_7,y=Night_Scale_7))+
    geom_raster(aes(fill=Prediction_back))+
    facet_grid(.~Treatment)+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = mid.pointss, space = "Lab",
                         name="Root growth rate (final root length/days of growth)") +
    theme_minimal()+
    ylab("Night Temperature (C)")+
    xlab("Day Temperature (C)")+
    theme(legend.position = "bottom")+
    labs(title = paste0(s))
  
  
  
  
  CairoPNG(file.path(paste0("Figures/",s,"_Root_growth_rate.png")), width = 12, height = 7, units = "in",
           res = 300)
  print(plots)
  dev.off()
  plots1[[s]] <- plots
}   

# Arrange in 2x2 grid
(plots1[["POSE"]] | plots1[["ACMI"]]) /
  (plots1[["ARTR"]] | plots1[["ELEL"]])

#############
##Emerg plots

emergplotlist <- list()

for (s in species){
  pnew.data.emerg <- prediction.list.emerg[[s]]
  temp_agg2 <- agg.list2[[s]]
  
  # for (i in unique(temp_agg2$Group)){
  #   tempp<-temp_agg2[which(temp_agg2$Group==i),]
  #   if (all(tempp$days_til_germ==0)) {
  #     print (i)
  #   }
  # }
  
  mid.pointss<-mean(temp_agg2$days_til_emerg[which(temp_agg2$Treatment=="C")],na.rm=T)
  emerg.plot<- ggplot(pnew.data.emerg,aes(x=Day_Scale_7,y=Night_Scale_7))+
    geom_raster(aes(fill=Prediction_back))+
    facet_grid(.~Treatment)+
    scale_fill_gradient2(low = "red", high = "blue", mid = "white",
                         midpoint = mid.pointss, space = "Lab",
                         name="Emergence Rate (1/days til emergence)") +
    theme_minimal()+
    ylab("Night Temperature (C)")+
    xlab("Day Temperature (C)")+
    theme(legend.position = "bottom")+
    labs(title = paste0(s))
  
  CairoPNG(file.path(paste0("Figures/",s,"_Days_til_Emerg_raster.png")), width = 12, height = 7, units = "in",
           res = 300)
  print(emerg.plot)
  dev.off()
  emergplotlist[[s]]<-emerg.plot
  
}

# Arrange in 2x2 grid
(emergplotlist[["POSE"]] | emergplotlist[["ACMI"]]) /
  (emergplotlist[["ARTR"]] | emergplotlist[["ELEL"]])

#############
## Germ rates

germplots <- list()

for (s in species){
  print (s)
  pnew.data <- prediction.list[[s]]
  temp_agg2 <- agg.list2[[s]]
  
  # for (i in unique(temp_agg2$Group)){
  #   tempp<-temp_agg2[which(temp_agg2$Group==i),]
  #   if (all(tempp$days_til_germ==0)) {
  #     print (i)
  #   }
  # }
  
  mid.pointss<-mean(temp_agg2$days_til_germ[which(temp_agg2$Treatment=="C")],na.rm=T)
  print(mid.pointss)
  germplot <- ggplot(pnew.data,aes(x=Day_Scale_7,y=Night_Scale_7))+
    geom_raster(aes(fill=Prediction_back))+ #, interpolate=T
    facet_grid(.~Treatment)+
    scale_fill_gradient2(low = "red", high = "blue", mid = "white",
                         midpoint = mid.pointss, space = "Lab",
                         name="Germination Rate (1/days til germ)") +
    theme_minimal()+
    theme(legend.position = "bottom")+
    ylab("Night Temperature (C)")+
    xlab("Day Temperature (C)")+
    labs(title = paste0(s))
  
  
  CairoPNG(file.path(paste0("Figures/",s,"_Days_til_Germ_raster.png")), width = 12, height = 7, units = "in",
           res = 300)
  print(germplot)
  dev.off()
  
  germplots[[s]] <- germplot
  
}

# Arrange in 2x2 grid
(germplots[["POSE"]] | germplots[["ACMI"]]) /
  (germplots[["ARTR"]] | germplots[["ELEL"]])



