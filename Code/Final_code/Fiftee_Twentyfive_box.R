##PLOTS

subset_listgerm_night <- list()

for (s in species){
  
  sub_agg_germ <- short.data.list2.Germ[[s]]
  
  high_range_germ_night <- subset.data.frame(sub_agg_germ, Night_Scale_7 > 15)
  
  subset_listgerm_night[[s]] <- high_range_germ_night
}

subset_listgerm_day <- list()

for (s in species){
  
  sub_agg_germ <- short.data.list2.Germ[[s]]
  
  high_range_germ_day <- subset.data.frame(sub_agg_germ, Day_Scale_7 > 15)
  
  subset_listgerm_day[[s]] <- high_range_germ_day
}

subset_listemerg_night <- list()

for (s in species){
  
  sub_agg_emerg <- short.data.list2.Emerg[[s]]
  
  high_range_emerg_night <- subset.data.frame(sub_agg_emerg, Night_Scale_7 > 15)
  
  subset_listemerg_night[[s]] <- high_range_emerg_night
}

subset_listemerg_day <- list()

for (s in species){
  
  sub_agg_emerg <- short.data.list2.Emerg[[s]]
  
  high_range_emerg_day <- subset.data.frame(sub_agg_emerg, Day_Scale_7 > 15)
  
  subset_listemerg_day[[s]] <- high_range_emerg_day
}

#Germination x Night

plotshigh <- list()

for (s in species){
  
  temp_agg3 <- subset_listgerm_night[[s]] 
  
  night.t.germ.boxplots <- ggplot(temp_agg3, aes(factor(round(Night_Scale_7,1)), days_til_germ,fill=Treatment)) +
    geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,linewidth=0.5)+ ##put the line int he background
    geom_boxplot() +
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0(s),
          x = "Night Temperature",
          y = "Days Until Germination")+
    theme_minimal(base_size = 10)+
    scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))#+
  #theme(legend.position="bottom")
  
  plotshigh[[s]] <- night.t.germ.boxplots
}  


Night_temp_germ_4boxhigh <- plotshigh[["POSE"]] + plotshigh[["ACMI"]] + 
  plotshigh[["ARTR"]] + plotshigh[["ELEL"]]+
  plot_layout(ncol = 2, axis_titles = 'collect', guides = 'collect') &
  theme(legend.position="bottom")
#plot_annotation(title = "Days Until Germination x Night Scale")

ggsave("Figures/RootDays_NightTemp_Boxhigh.png", Night_temp_germ_4boxhigh)


#Emergence x Night

plotshigh1 <- list()

for (s in species){
  
  temp_agg3 <- subset_listemerg_night[[s]]
  
  night.t.emerg.boxplots <- ggplot(temp_agg3, aes(factor(round(Night_Scale_7,1)), days_til_emerg,fill=Treatment)) +
    geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,linewidth=0.5)+ ##put the line int he background
    geom_boxplot() +
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0(s),
          x = "Night Temperature",
          y = "Days Until Emergence")+
    theme_minimal(base_size = 10)+
    scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    theme(legend.position="bottom")
  
  plotshigh1[[s]] <- night.t.emerg.boxplots
  
}  

Night_temp_emerg_4boxhigh <- plotshigh1[["POSE"]] + plotshigh1[["ACMI"]] + 
  plotshigh1[["ARTR"]] + plotshigh1[["ELEL"]]+
  plot_layout(ncol= 2, axis_titles = 'collect', guides = 'collect')&
  theme(legend.position="bottom")
#plot_annotation(title = "Days Until Emergence x Night Scale")

ggsave("Figures/ShootDays_NightTemp_Boxhigh.png", Night_temp_emerg_4boxhigh)

##DAY
# Germination x Day

plotshigh2 <- list()

for (s in species){
  
  temp_agg3 <- subset_listgerm_day[[s]]
  
  day.t.germ.boxplots <- ggplot(temp_agg3, aes(factor(round(Day_Scale_7,1)), days_til_germ,fill=Treatment)) +
    geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,linewidth=0.5)+ ##put the line int he background
    geom_boxplot() +
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0(s),
          x = "Day Temperature",
          y = "Days Until Germination")+
    theme_minimal(base_size = 10)+
    scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    theme(legend.position="bottom")
  
  plotshigh2[[s]] <- day.t.germ.boxplots
}  

Day_temp_germ_4boxhigh <- plotshigh2[["POSE"]] + plotshigh2[["ACMI"]] + 
  plotshigh2[["ARTR"]] + plotshigh2[["ELEL"]]+
  plot_layout(ncol = 2, axis_titles = 'collect', guides = 'collect')&
  theme(legend.position="bottom")
#plot_annotation(title = "Days Until Germination x Day Scale")

ggsave("Figures/RootDays_DayTemp_Boxhigh.png", Day_temp_germ_4boxhigh)

# Emergence x Day

plotshigh3 <- list()

for (s in species){
  
  temp_agg3 <- subset_listemerg_day[[s]]
  
  day.t.emerg.boxplots <- ggplot(temp_agg3, aes(factor(round(Day_Scale_7,1)), days_til_emerg,fill=Treatment)) +
    geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,linewidth=0.5)+ ##put the line int he background
    geom_boxplot() +
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0(s),
          x = "Day Temperature",
          y = "Days Until Emergence")+
    theme_minimal(base_size = 10)+
    scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    theme(legend.position="bottom")
  
  plotshigh3[[s]] <- day.t.emerg.boxplots
  
}  

Day_temp_emerg_4boxhigh <- plotshigh3[["POSE"]] + plotshigh3[["ACMI"]] + 
  plotshigh3[["ARTR"]] + plotshigh3[["ELEL"]]+
  plot_layout(ncol= 2, axis_titles = 'collect', guides = 'collect')&
  theme(legend.position="bottom")
#plot_annotation(title = "Days Until Emergence x Day Scale")

ggsave("Figures/ShootDays_DayTemp_Boxhigh.png", Day_temp_emerg_4boxhigh)

#######
##RATES

subset_listgrate_night <- list()

for (s in species){
  
  sub_agg_germ <- agg.list.sub.RGrate[[s]]
  
  high_range_grate_night <- subset.data.frame(sub_agg_germ, Night_Scale_7 > 15)
  
  subset_listgrate_night[[s]] <- high_range_grate_night
}

subset_listgrate_day <- list()

for (s in species){
  
  sub_agg_germ <- agg.list.sub.RGrate[[s]]
  
  high_range_grate_day <- subset.data.frame(sub_agg_germ, Day_Scale_7 > 15)
  
  subset_listgrate_day[[s]] <- high_range_grate_day
}

subset_listerate_night <- list()

for (s in species){
  
  sub_agg_erate <- agg.list.sub.EGrate[[s]]
  
  high_range_erate_night <- subset.data.frame(sub_agg_erate, Night_Scale_7 > 15)
  
  subset_listerate_night[[s]] <- high_range_erate_night
}

subset_listerate_day <- list()

for (s in species){
  
  sub_agg_erate <- agg.list.sub.EGrate[[s]]
  
  high_range_erate_day <- subset.data.frame(sub_agg_erate, Day_Scale_7 > 15)
  
  subset_listerate_day[[s]] <- high_range_erate_day
}
### Post Emergence/Germination Growth Rates
#Night Scale

#Germination Rate x Night Scale

plotshigh4 <- list()

for (s in species){
  
  temp_agg3 <- agg.list.sub.RGrate[[s]] 
  
  night.t.rgrate.boxplots <- ggplot(temp_agg3, aes(factor(round(Night_Scale_7,1)), Avg_Root_rate,fill=Treatment)) +
    geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,linewidth=0.5)+ ##put the line int he background
    geom_boxplot() +
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0(s),
          x = "Night Temperature",
          y = "Root Growth Rate")+
    theme_minimal(base_size = 10)+
    scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    theme(legend.position="bottom")
  
  plotshigh4[[s]] <- night.t.rgrate.boxplots
  
}  

Night_temp_rgrate_4boxhigh <- plotshigh4[["POSE"]] + plotshigh4[["ACMI"]] + 
  plotshigh4[["ARTR"]] + plotshigh4[["ELEL"]]+
  plot_layout(ncol= 2, axis_titles = 'collect', guides = 'collect')&
  theme(legend.position="bottom")
#plot_annotation(title = "Root Growth Rate x Night Scale")

ggsave("Figures/Rootrate_NightTemp_Boxhigh.png", Night_temp_rgrate_4boxhigh)


# Shoot Growth Rate x Night Scale

plotshigh5 <- list()

for (s in species){
  
  temp_agg3 <- agg.list.sub.EGrate[[s]] 
  
  night.t.egrate.boxplots <- ggplot(temp_agg3, aes(factor(round(Night_Scale_7,1)), Avg_Root_rate,fill=Treatment)) +
    geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,linewidth=0.5)+ ##put the line int he background
    geom_boxplot() +
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0(s),
          x = "Night Temperature",
          y = "Shoot Growth Rate")+
    theme_minimal(base_size = 10)+
    scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    theme(legend.position="bottom")
  
  plotshigh5[[s]] <- night.t.egrate.boxplots
  
} 

Night_temp_egrate_4boxhigh <- plotshigh5[["POSE"]] + plotshigh5[["ACMI"]] + 
  plotshigh5[["ARTR"]] + plotshigh5[["ELEL"]]+
  plot_layout(ncol= 2, axis_titles = 'collect', guides = 'collect')&
  theme(legend.position="bottom")
#plot_annotation(title = "Shoot Growth Rate x Night Scale")

ggsave("Figures/Shootrate_NightTemp_Boxhigh.png", Night_temp_egrate_4boxhigh)

#Day Scale

# Root Growth Rate x Day Scale

plotshigh6 <- list()

for (s in species){
  
  temp_agg3 <- agg.list.sub.RGrate[[s]] 
  
  day.t.rgrate.boxplots <- ggplot(temp_agg3, aes(factor(round(Day_Scale_7,1)), Avg_Root_rate,fill=Treatment)) +
    geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,linewidth=0.5)+ ##put the line int he background
    geom_boxplot() +
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0(s),
          x = "Day Temperature",
          y = "Root Growth Rate")+
    theme_minimal(base_size = 10)+
    scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    theme(legend.position="bottom")
  
  plotshigh6[[s]] <- day.t.rgrate.boxplots
  
}  

Day_temp_rgrate_4boxhigh <- plotshigh6[["POSE"]] + plotshigh6[["ACMI"]] + 
  plotshigh6[["ARTR"]] + plotshigh6[["ELEL"]]+
  plot_layout(ncol= 2, axis_titles = 'collect', guides = 'collect')&
  theme(legend.position="bottom")
#plot_annotation(title = "Root Growth Rate x Day Scale")

ggsave("Figures/Rootrate_DayTemp_Boxhigh.png", Day_temp_rgrate_4boxhigh)

# Shoot Growth Rate x Day Scale

plotshigh7 <- list()

for (s in species){
  
  temp_agg3 <- agg.list.sub.EGrate[[s]] 
  
  day.t.egrate.boxplots <- ggplot(temp_agg3, aes(factor(round(Day_Scale_7,1)), Avg_Root_rate,fill=Treatment)) +
    geom_smooth(se=F,aes(group=Treatment,col=Treatment),span = 1.5,linewidth=0.5)+ ##put the line in the background
    geom_boxplot() +
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0(s),
          x = "Day Temperature",
          y = "Shoot Growth Rate")+
    theme_minimal(base_size = 10)+
    scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    theme(legend.position="bottom")
  
  plotshigh7[[s]] <- day.t.egrate.boxplots
  
} 

Day_temp_egrate_4boxhigh <- plotshigh7[["POSE"]] + plotshigh7[["ACMI"]] + 
  plotshigh7[["ARTR"]] + plotshigh7[["ELEL"]]+
  plot_layout(ncol= 2, axis_titles = 'collect', guides = 'collect')&
  theme(legend.position="bottom")
#plot_annotation(title = "Shoot Growth Rate x Day Scale")

ggsave("Figures/Shootrate_DayTemp_Boxhigh.png", Day_temp_egrate_4boxhigh)