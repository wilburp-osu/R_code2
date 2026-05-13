##PLOTS

#Germination x Night

plots <- list()

for (s in species){
  
  temp_agg3 <- short.data.list2.Germ[[s]] 
  
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
  
  plots[[s]] <- night.t.germ.boxplots
}  

Night_temp_germ_4box <- plots[["POSE"]] + plots[["ACMI"]] + plots[["ARTR"]] +
  plots[["ELEL"]]+
  plot_layout(ncol = 2, axis_titles = 'collect', guides = 'collect') &
  theme(legend.position="bottom")
  #plot_annotation(title = "Days Until Germination x Night Scale")

ggsave("Figures/RootDays_NightTemp_Box.png", Night_temp_germ_4box)


#Emergence x Night

plots1 <- list()

for (s in species){
  
  temp_agg3 <- short.data.list2.Emerg[[s]] 
  
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
  
  plots1[[s]] <- night.t.emerg.boxplots
  
}  

Night_temp_emerg_4box <- plots1[["POSE"]] + plots1[["ACMI"]] + 
  plots1[["ARTR"]] + plots1[["ELEL"]]+
  plot_layout(ncol= 2, axis_titles = 'collect', guides = 'collect')&
  theme(legend.position="bottom")
  #plot_annotation(title = "Days Until Emergence x Night Scale")

ggsave("Figures/ShootDays_NightTemp_Box.png", Night_temp_emerg_4box)

##DAY
# Germination x Day

plots2 <- list()

for (s in species){
  
  temp_agg3 <- short.data.list2.Germ[[s]] 
  
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
  
  plots2[[s]] <- day.t.germ.boxplots
}  

Day_temp_germ_4box <- plots2[["POSE"]] + plots2[["ACMI"]] + plots2[["ARTR"]] +
  plots2[["ELEL"]]+
  plot_layout(ncol = 2, axis_titles = 'collect', guides = 'collect')&
  theme(legend.position="bottom")
  #plot_annotation(title = "Days Until Germination x Day Scale")

ggsave("Figures/RootDays_DayTemp_Box.png", Day_temp_germ_4box)

# Emergence x Day

plots3 <- list()

for (s in species){
  
  temp_agg3 <- short.data.list2.Emerg[[s]] 
  
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
  
  plots3[[s]] <- day.t.emerg.boxplots
  
}  

Day_temp_emerg_4box <- plots3[["POSE"]] + plots3[["ACMI"]] + 
  plots3[["ARTR"]] + plots3[["ELEL"]]+
  plot_layout(ncol= 2, axis_titles = 'collect', guides = 'collect')&
  theme(legend.position="bottom")
  #plot_annotation(title = "Days Until Emergence x Day Scale")

ggsave("Figures/ShootDays_DayTemp_Box.png", Day_temp_emerg_4box)

### Post Emergence/Germination Growth Rates
#Night Scale

#Germination Rate x Night Scale

plots4 <- list()

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
  
  plots4[[s]] <- night.t.rgrate.boxplots
  
}  

Night_temp_rgrate_4box <- plots4[["POSE"]] + plots4[["ACMI"]] + 
  plots4[["ARTR"]] + plots4[["ELEL"]]+
  plot_layout(ncol= 2, axis_titles = 'collect', guides = 'collect')&
  theme(legend.position="bottom")
  #plot_annotation(title = "Root Growth Rate x Night Scale")

ggsave("Figures/Rootrate_NightTemp_Box.png", Night_temp_rgrate_4box)


# Shoot Growth Rate x Night Scale

plots5 <- list()

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
  
  plots5[[s]] <- night.t.egrate.boxplots
  
} 
 
Night_temp_egrate_4box <- plots5[["POSE"]] + plots5[["ACMI"]] + 
  plots5[["ARTR"]] + plots5[["ELEL"]]+
  plot_layout(ncol= 2, axis_titles = 'collect', guides = 'collect')&
  theme(legend.position="bottom")
  #plot_annotation(title = "Shoot Growth Rate x Night Scale")

ggsave("Figures/Shootrate_NightTemp_Box.png", Night_temp_egrate_4box)

#Day Scale

# Root Growth Rate x Day Scale

plots6 <- list()

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
  
  plots6[[s]] <- day.t.rgrate.boxplots
  
}  

Day_temp_rgrate_4box <- plots6[["POSE"]] + plots6[["ACMI"]] + 
  plots6[["ARTR"]] + plots6[["ELEL"]]+
  plot_layout(ncol= 2, axis_titles = 'collect', guides = 'collect')&
  theme(legend.position="bottom")
  #plot_annotation(title = "Root Growth Rate x Day Scale")

ggsave("Figures/Rootrate_DayTemp_Box.png", Day_temp_rgrate_4box)

# Shoot Growth Rate x Day Scale

plots7 <- list()

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
  
  plots7[[s]] <- day.t.egrate.boxplots
  
} 

Day_temp_egrate_4box <- plots7[["POSE"]] + plots7[["ACMI"]] + 
  plots7[["ARTR"]] + plots7[["ELEL"]]+
  plot_layout(ncol= 2, axis_titles = 'collect', guides = 'collect')&
  theme(legend.position="bottom")
  #plot_annotation(title = "Shoot Growth Rate x Day Scale")

ggsave("Figures/Shootrate_DayTemp_Box.png", Day_temp_egrate_4box)

# Big_Plot <- Day_temp_germ_4box + Day_temp_emerg_4box + Night_temp_emerg_4box +
#   Night_temp_germ_4box + plot_layout(ncol=2, nrow=2)
# 
# Big_Plot
