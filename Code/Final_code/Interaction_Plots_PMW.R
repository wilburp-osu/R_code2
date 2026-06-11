################################3
####Days Until Germ
#Only Significant for ACMI###

intplots <-list()

for (s in species){
  temp.interaction<-short.data.list2.Germ[[s]]
  temp.interaction$interaction_axis<- log(temp.interaction$Day_Scale_7/temp.interaction$Night_Scale_7)
  
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.9737892] <- "2.648~3.909"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.5842735 & temp.interaction$interaction_axis < 0.9737892] <- "1.794~2.648"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.1947578 & temp.interaction$interaction_axis < 0.5842735] <- "1.215~1.794"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.1947578 & temp.interaction$interaction_axis < 0.1947578] <- "0.823~1.215"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.5842735 & temp.interaction$interaction_axis < -0.1947578] <- "0.558~0.823"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.9737892 & temp.interaction$interaction_axis < -0.5842735] <- "0.378~0.558"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis < -0.9737892] <- "0.256~0.378"
  
  Temp.int.germ.box <- ggplot(temp.interaction, aes(x = interaction_groups, y = days_til_germ, fill = Treatment))+
    geom_smooth(#method="lm",
      se=F,aes(group=Treatment,col=Treatment),linewidth=0.5,span=1.5)+ ##put the line int he background
    geom_boxplot()+
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 4, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"),
                 show.legend = F)+
    labs( title = paste0(s),
          x = "Day Temperature:Night Temperature Ratio",
          y = "Days Until Germination")+
    theme_minimal(base_size = 10)+
    scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    theme(legend.position="bottom", axis.text.x = element_text(angle = 30, vjust = 1, hjust = 0.9, color = "black"))
 
  intplots[[s]] <- Temp.int.germ.box
  
}

Temp_germ_int_4box <- intplots[["POSE"]] + intplots[["ACMI"]] + 
  intplots[["ARTR"]] + intplots[["ELEL"]]+
  plot_layout(ncol = 2, axis_titles = 'collect', guides = 'collect') &
  theme(legend.position="bottom")
  #plot_annotation(title = "Days Until Germ x Day Scale:Night Scale")

ggsave("Figures/RootDays_TempInt_Box.png", Temp_germ_int_4box)

intplotspp <-list()

for (s in species){
  temp.interactionpp<-short.data.list2.Germ[[s]]
  temp.interactionpp$interaction_axis<- log(temp.interaction$Day_Scale_7/temp.interaction$Night_Scale_7)
  
  temp.interactionpp$interaction_groups[temp.interactionpp$interaction_axis >= 0.9737892] <- "Hot Day, Cold Night"
  temp.interactionpp$interaction_groups[temp.interactionpp$interaction_axis >= 0.5842735 & temp.interaction$interaction_axis < 0.9737892] <- "Warm Day, Cool Night"
  temp.interactionpp$interaction_groups[temp.interactionpp$interaction_axis >= 0.1947578 & temp.interaction$interaction_axis < 0.5842735] <- "Tepid Day, Chilly Night"
  temp.interactionpp$interaction_groups[temp.interactionpp$interaction_axis >= -0.1947578 & temp.interaction$interaction_axis < 0.1947578] <- "Similar Day & Night"
  temp.interactionpp$interaction_groups[temp.interactionpp$interaction_axis >= -0.5842735 & temp.interaction$interaction_axis < -0.1947578] <- "Chilly Day, Tepid Night"
  temp.interactionpp$interaction_groups[temp.interactionpp$interaction_axis >= -0.9737892 & temp.interaction$interaction_axis < -0.5842735] <- "Cool Day, Warm Night"
  temp.interactionpp$interaction_groups[temp.interactionpp$interaction_axis < -0.9737892] <- "Cold Day, Hot Night"
  
  Temp.int.germ.boxpp <- ggplot(temp.interactionpp, aes(x = interaction_groups, y = days_til_germ, fill = Treatment))+
    geom_smooth(#method="lm",
      se=F,aes(group=Treatment,col=Treatment),linewidth=1.5,span=1.5)+ ##put the line int he background
    geom_boxplot()+
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 4, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"),
                 show.legend = F)+
    labs(title = paste0(s),
          x = "Day Temperature:Night Temperature Ratio",
          y = "Days Until Germination")+
    theme_minimal(base_size = 10)+
    scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    theme(legend.position="bottom", axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, color = "black"),
          axis.title.x= element_text(vjust = -2))
  
  intplotspp[[s]] <- Temp.int.germ.boxpp
  
}

Temp_germ_int_4box_powerpoint <- intplotspp[["POSE"]] + intplotspp[["ACMI"]] + 
  intplotspp[["ARTR"]] + intplotspp[["ELEL"]]+
  plot_layout(ncol = 2, axis_titles = 'collect', guides = 'collect') &
  theme(legend.position="bottom")
#plot_annotation(title = "Days Until Germ x Day Scale:Night Scale")

ggsave("Figures/RootDays_TempInt_Box_powerpoint.png", Temp_germ_int_4box_powerpoint)

###################
###Days til Emerg
##Significant for NONE

intplots1 <- list()

for (s in species){
  temp.interaction<-short.data.list2.Emerg[[s]]
  temp.interaction$interaction_axis<- log(temp.interaction$Day_Scale_7/temp.interaction$Night_Scale_7)
  
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.9737892] <- "2.648~3.909"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.5842735 & temp.interaction$interaction_axis < 0.9737892] <- "1.794~2.648"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.1947578 & temp.interaction$interaction_axis < 0.5842735] <- "1.215~1.794"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.1947578 & temp.interaction$interaction_axis < 0.1947578] <- "0.823~1.215"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.5842735 & temp.interaction$interaction_axis < -0.1947578] <- "0.558~0.823"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.9737892 & temp.interaction$interaction_axis < -0.5842735] <- "0.378~0.558"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis < -0.9737892] <- "0.256~0.378"
  
  Temp.int.emerg.box <- ggplot(temp.interaction, aes(x = interaction_groups, y = days_til_emerg, fill = Treatment))+
    geom_smooth(#method="lm",
      se=F,aes(group=Treatment,col=Treatment),linewidth=0.5, span=1.5)+ ##put the line int he background
    geom_boxplot()+
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"),
                 show.legend = F)+
    labs( title = paste0(s),
          x = "Day Temperature:Night Temperature Ratio",
          y = "Days Until Emergence")+
    theme_minimal(base_size = 10)+
    scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4" ))+
    theme(legend.position="bottom", axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, color = "black"))
  
  intplots1[[s]] <- Temp.int.emerg.box
  
}

Temp_emerg_int_4box <- intplots1[["POSE"]] + intplots1[["ACMI"]] + 
  intplots1[["ARTR"]] + intplots1[["ELEL"]]+
  plot_layout(ncol = 2, axis_titles = 'collect', guides = 'collect') &
  theme(legend.position="bottom")
  #plot_annotation(title = "Days Until Emerg x Day Scale:Night Scale")

ggsave("Figures/ShootDays_TempInt_Box.png", Temp_emerg_int_4box)

###################
###Germ rate
##Significant for ACMI and ELEL

intplots2 <- list()

for (s in species){
  temp.interaction<-agg.list.sub.RGrate[[s]]
  temp.interaction$interaction_axis<- log(temp.interaction$Day_Scale_7/temp.interaction$Night_Scale_7)
  
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.9737892] <- "2.648~3.909"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.5842735 & temp.interaction$interaction_axis < 0.9737892] <- "1.794~2.648"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.1947578 & temp.interaction$interaction_axis < 0.5842735] <- "1.215~1.794"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.1947578 & temp.interaction$interaction_axis < 0.1947578] <- "0.823~1.215"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.5842735 & temp.interaction$interaction_axis < -0.1947578] <- "0.558~0.823"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.9737892 & temp.interaction$interaction_axis < -0.5842735] <- "0.378~0.558"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis < -0.9737892] <- "0.256~0.378"
  
  Temp.int.rgrate.box <- ggplot(temp.interaction, aes(x = interaction_groups, y = Avg_Root_rate, fill = Treatment))+
    geom_smooth(#method="lm",
      se=F,aes(group=Treatment,col=Treatment),linewidth=0.5, span = 1.5)+ ##put the line int he background
    geom_boxplot()+
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"),
                 show.legend = F)+
    labs( title = paste0(s),
          x = "Day Temperature:Night Temperature Ratio",
          y = "Root Growth Rate \n")+
    theme_minimal(base_size = 10)+
    scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    theme(legend.position="bottom", axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, color = "black"))
  
  intplots2[[s]] <- Temp.int.rgrate.box
  
}

Temp_rgrate_int_4box <- intplots2[["POSE"]] + intplots2[["ACMI"]] + 
  intplots2[["ARTR"]] + intplots2[["ELEL"]]+
  plot_layout(ncol = 2, axis_titles = 'collect', guides = 'collect') &
  theme(legend.position="bottom")
  #plot_annotation(title = "Root Growth Rate x Day Scale:Night Scale")

ggsave("Figures/RootRate_TempInt_Box.png", Temp_rgrate_int_4box)

##################
###Emerg rate
##Significant for ELEL                     

intplots3 <- list()

for (s in species){
  temp.interaction<-agg.list.sub.EGrate[[s]]
  temp.interaction$interaction_axis<- log(temp.interaction$Day_Scale_7/temp.interaction$Night_Scale_7)
  
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.9737892] <- "2.648~3.909"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.5842735 & temp.interaction$interaction_axis < 0.9737892] <- "1.794~2.648"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.1947578 & temp.interaction$interaction_axis < 0.5842735] <- "1.215~1.794"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.1947578 & temp.interaction$interaction_axis < 0.1947578] <- "0.823~1.215"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.5842735 & temp.interaction$interaction_axis < -0.1947578] <- "0.558~0.823"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.9737892 & temp.interaction$interaction_axis < -0.5842735] <- "0.378~0.558"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis < -0.9737892] <- "0.256~0.378"
  
  Temp.int.egrate.box <- ggplot(temp.interaction, aes(x = interaction_groups, y = Avg_Root_rate, fill = Treatment))+
    geom_smooth(#method="lm",
      se=F,aes(group=Treatment,col=Treatment),linewidth=0.5, span=1.5)+ ##put the line int he background
    geom_boxplot()+
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0(s),
          x = "Day Temperature:Night Temperature Ratio",
          y = "Shoot Growth Rate \n")+
    theme_minimal(base_size = 10)+
    scale_color_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    scale_fill_manual(values=c("#CCCCCC", "#33a02c","#b2df8a","#1f78b4"))+
    theme(legend.position="bottom", axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, color = "black"))
  
  intplots3[[s]] <- Temp.int.egrate.box
  
}

Temp_egrate_int_4box <- intplots3[["POSE"]] + intplots3[["ACMI"]] + 
  intplots3[["ARTR"]] + intplots3[["ELEL"]]+
  plot_layout(ncol = 2, axis_titles = 'collect', guides = 'collect') &
  theme(legend.position="bottom")
  #plot_annotation(title = "Shoot Growth Rate x Day Scale:Night Scale")

ggsave("Figures/ShootRate_TempInt_Box.png", Temp_egrate_int_4box)