##testing things

test<-short.data.list2[[1]]
test$interaction_axis<- log(test$Day_Scale_7/test$Night_Scale_7)

#binning
test$interaction_groups <- c("Low_High","Med_High", "High_Med", "High_Low" )

test$interaction_groups[test$interaction_axis >= 0.75] <- "High_Low"
test$interaction_groups[test$interaction_axis >= 0 & test$interaction_axis < 0.75] <- "High_Med"
test$interaction_groups[test$interaction_axis >= -0.75 & test$interaction_axis < 0] <- "Med_High"
test$interaction_groups[test$interaction_axis < -0.75] <- "Low_High"

#finding intervals for bins
minint <- min(test$interaction_axis)
maxint <- max(test$interaction_axis)

rngint <- maxint - minint

binsint <-rngint/7

groups <- c(maxint, maxint-binsint, maxint-(2*binsint), maxint-(3*binsint), maxint-(4*binsint), 
  maxint-(5*binsint), maxint-(6*binsint), maxint-(7*binsint))

exp(groups)

###Days til Germ ##TEM just edited this one.

for (s in species){
  temp.interaction<-short.data.list2[[s]]
  temp.interaction$interaction_axis<- log(temp.interaction$Day_Scale_7/temp.interaction$Night_Scale_7)

  
  #make labels for temp ratios
  temp.interaction$interaction_groups <- c("0.256~0.378 \n Cold Day, Hot Night","0.378~0.558 ", "0.558~0.823", "0.823~1.215",
                               "1.215~1.794", "1.794~2.648", "2.648~3.909")
  #fill with intervals 
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.9737892] <- "2.648~3.909 \n Hot Day, Cold Night"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.5842735 & temp.interaction$interaction_axis < 0.9737892] <- "1.794~2.648 \n Warm Day, Cool Night"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.1947578 & temp.interaction$interaction_axis < 0.5842735] <- "1.215~1.794 \n Tepid Day, Chilly Night"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.1947578 & temp.interaction$interaction_axis < 0.1947578] <- "0.823~1.215 \n Similar Day & Night"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.5842735 & temp.interaction$interaction_axis < -0.1947578] <- "0.558~0.823 \n Chilly Day, Tepid Night"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.9737892 & temp.interaction$interaction_axis < -0.5842735] <- "0.378~0.558 \n Cool Day, Warm Night"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis < -0.9737892] <- "0.256~0.378 \n Cold Day, Hot Night"

  Temp.int.germ.box <- ggplot(temp.interaction, aes(x = interaction_groups, y = days_til_germ, fill = Treatment))+
    geom_smooth(#method="lm",
      se=F,aes(group=Treatment,col=Treatment),size=1.5,span=1.5)+ ##put the line int he background
    geom_boxplot()+
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",stroke=2,
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0("Days Until Germ x Day Scale:Night Scale (", s, ")"),
          x = "Day Temperature:Night Temperature Ratio",
          y = "Days Until Germination")+
    theme_minimal(base_size = 15)+
    scale_color_manual(values=c("grey60","green3","orchid3","dodgerblue3"))+
    scale_fill_manual(values=c("grey80","green","orchid1","skyblue"))
    
    print(Temp.int.germ.box)

}

###Days til Emerg

for (s in species){
  temp.interaction<-short.data.list2[[s]]
  temp.interaction$interaction_axis<- log(temp.interaction$Day_Scale_7/temp.interaction$Night_Scale_7)
  
  temp.interaction$interaction_groups <- c("0.256~0.378","0.378~0.558", "0.558~0.823", "0.823~1.215",
                                           "1.215~1.794", "1.794~2.648", "2.648~3.909")
  
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.9737892] <- "2.648~3.909"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.5842735 & temp.interaction$interaction_axis < 0.9737892] <- "1.794~2.648"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.1947578 & temp.interaction$interaction_axis < 0.5842735] <- "1.215~1.794"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.1947578 & temp.interaction$interaction_axis < 0.1947578] <- "0.823~1.215"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.5842735 & temp.interaction$interaction_axis < -0.1947578] <- "0.558~0.823"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.9737892 & temp.interaction$interaction_axis < -0.5842735] <- "0.378~0.558"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis < -0.9737892] <- "0.256~0.378"
  
  Temp.int.emerg.box <- ggplot(temp.interaction, aes(x = interaction_groups, y = days_til_emerg, fill = Treatment))+
    geom_smooth(#method="lm",
      se=F,aes(group=Treatment,col=Treatment),size=1)+ ##put the line int he background
    geom_boxplot()+
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0("Days Until Emerg x Day Scale/Night Scale (", s, ")"),
          x = "Day Temperature / Night Temperature",
          y = "Days Until Emergence")+
    theme_minimal(base_size = 15)+
    scale_color_manual(values=c("grey60","green3","orchid3","dodgerblue3"))+
    scale_fill_manual(values=c("grey80","green","orchid1","skyblue"))
  
  print(Temp.int.emerg.box)
  
}

###Germ rate

for (s in species){
  temp.interaction<-short.data.list2[[s]]
  temp.interaction$interaction_axis<- log(temp.interaction$Day_Scale_7/temp.interaction$Night_Scale_7)
  
  temp.interaction$interaction_groups <- c("0.256~0.378","0.378~0.558", "0.558~0.823", "0.823~1.215",
                                           "1.215~1.794", "1.794~2.648", "2.648~3.909")
  
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.9737892] <- "2.648~3.909"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.5842735 & temp.interaction$interaction_axis < 0.9737892] <- "1.794~2.648"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.1947578 & temp.interaction$interaction_axis < 0.5842735] <- "1.215~1.794"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.1947578 & temp.interaction$interaction_axis < 0.1947578] <- "0.823~1.215"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.5842735 & temp.interaction$interaction_axis < -0.1947578] <- "0.558~0.823"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.9737892 & temp.interaction$interaction_axis < -0.5842735] <- "0.378~0.558"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis < -0.9737892] <- "0.256~0.378"
  
  Temp.int.grate.box <- ggplot(temp.interaction, aes(x = interaction_groups, y = germ_rate, fill = Treatment))+
    geom_smooth(#method="lm",
      se=F,aes(group=Treatment,col=Treatment),size=1)+ ##put the line int he background
    geom_boxplot()+
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0("Germ Rate x Day Scale/Night Scale (", s, ")"),
          x = "Day Temperature / Night Temperature",
          y = "Germination Rate")+
    theme_minimal(base_size = 15)+
    scale_color_manual(values=c("grey60","green3","orchid3","dodgerblue3"))+
    scale_fill_manual(values=c("grey80","green","orchid1","skyblue"))
  
  print(Temp.int.grate.box)
  
}

###Emerg rate

for (s in species){
  temp.interaction<-short.data.list2[[s]]
  temp.interaction$interaction_axis<- log(temp.interaction$Day_Scale_7/temp.interaction$Night_Scale_7)
  
  temp.interaction$interaction_groups <- c("0.256~0.378","0.378~0.558", "0.558~0.823", "0.823~1.215",
                                           "1.215~1.794", "1.794~2.648", "2.648~3.909")
  
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.9737892] <- "2.648~3.909"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.5842735 & temp.interaction$interaction_axis < 0.9737892] <- "1.794~2.648"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.1947578 & temp.interaction$interaction_axis < 0.5842735] <- "1.215~1.794"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.1947578 & temp.interaction$interaction_axis < 0.1947578] <- "0.823~1.215"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.5842735 & temp.interaction$interaction_axis < -0.1947578] <- "0.558~0.823"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.9737892 & temp.interaction$interaction_axis < -0.5842735] <- "0.378~0.558"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis < -0.9737892] <- "0.256~0.378"
  
  Temp.int.erate.box <- ggplot(temp.interaction, aes(x = interaction_groups, y = emerg_rate, fill = Treatment))+
    geom_smooth(#method="lm",
      se=F,aes(group=Treatment,col=Treatment),size=1)+ ##put the line int he background
    geom_boxplot()+
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 24,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"))+
    labs( title = paste0("Emerg Rate x Day Scale/Night Scale (", s, ")"),
          x = "Day Temperature / Night Temperature",
          y = "Emergence Rate")+
    theme_minimal(base_size = 15)+
    scale_color_manual(values=c("grey60","green3","orchid3","dodgerblue3"))+
    scale_fill_manual(values=c("grey80","green","orchid1","skyblue"))
  
  print(Temp.int.erate.box)
  
}



