##testing things

##Checking Palettes
#Original

install.packages("colorblindcheck")
library(colorblindcheck)

pal1 = c("#00CD00","#CD69C9","#1874cd") #was able to find the hexes for the original color scheme

palette_check(pal1, plot = TRUE)

#New

Newpal = c("#33a02c","#b2df8a","#1f78b4") ##this would work best with grey
palette_check(Newpal, plot = T)

Newpal1 = c("#1b9e77","#d95f02","#7570b3")
palette_check(Newpal1, plot = T)

Newpal2 = c("#66c2a5","#fc8d62","#8da0cb")
palette_check(Newpal2, plot = T)

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

## Also using these to test color palettes
# retrieved palettes from colorbrewer

###Days til Germ
## I think this color palette could work best
# It does not pop a whole lot but avoids color blindness issues and helps to show that
# the "H" and "L" are related by being variations of the same hue
## also edited the mean points, This could work

for (s in species){
  temp.interaction<-short.data.list2[[s]]
  temp.interaction$interaction_axis<- log(temp.interaction$Day_Scale_7/temp.interaction$Night_Scale_7)

  
  #make labels for temp ratios
  temp.interaction$interaction_groups <- c("0.256~0.378","0.378~0.558", "0.558~0.823", "0.823~1.215",
                               "1.215~1.794", "1.794~2.648", "2.648~3.909")
  #fill with intervals 
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.9737892] <- "2.648~3.909"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.5842735 & temp.interaction$interaction_axis < 0.9737892] <- "1.794~2.648"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= 0.1947578 & temp.interaction$interaction_axis < 0.5842735] <- "1.215~1.794"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.1947578 & temp.interaction$interaction_axis < 0.1947578] <- "0.823~1.215"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.5842735 & temp.interaction$interaction_axis < -0.1947578] <- "0.558~0.823"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis >= -0.9737892 & temp.interaction$interaction_axis < -0.5842735] <- "0.378~0.558"
  temp.interaction$interaction_groups[temp.interaction$interaction_axis < -0.9737892] <- "0.256~0.378"

  Temp.int.germ.box <- ggplot(temp.interaction, aes(x = interaction_groups, y = days_til_germ, fill = Treatment))+
    geom_smooth(#method="lm",
      se=F,aes(group=Treatment,col=Treatment),size=1)+ ##put the line int he background
    geom_boxplot()+
    stat_boxplot(geom ='errorbar',col="black") +
    stat_summary(fun = "mean", geom = "point", size = 4, shape = 17,col="black",
                 aes(fill=Treatment),position = position_dodge2 (width = 0.75, preserve = "single"),
                 show.legend = F)+
    labs( title = paste0("Days Until Germ x Day Scale/Night Scale (", s, ")"),
          x = "Day Temperature / Night Temperature",
          y = "Days Until Germination")+
    theme_minimal(base_size = 15)+
    scale_color_manual(values=c("grey60","#33a02c","#b2df8a","#1f78b4"))+
    scale_fill_manual(values=c("grey80","#33a02c","#b2df8a","#1f78b4"))
    
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
    scale_color_manual(values=c("grey60","#1b9e77","#d95f02","#7570b3"))+
    scale_fill_manual(values=c("grey80","#1b9e77","#d95f02","#7570b3"))
  
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
    scale_color_manual(values=c("grey60","#66c2a5","#fc8d62","#8da0cb"))+
    scale_fill_manual(values=c("grey80","#66c2a5","#fc8d62","#8da0cb"))
  
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



