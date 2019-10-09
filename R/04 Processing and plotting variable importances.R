## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Processing and plotting measures of variable importance

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

rm(list=ls())
.libPaths("X:\\Old N drive\\R\\R-2.14.0\\library")

##Set working directory
setwd("X:\\Species rarity\\Rarity- Global analysis\\")

## Libraries
library(ggplot2)
library(viridis)
library(gridExtra)
library(grid)
library(png)
library(plyr)

## Results path
results.dir.path<-"X:\\Species rarity\\Rarity- Global analysis\\Model Output"

## Model iteration
mod<-"RF mods - main"

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##      Model Fit
## First establish which models have worked well
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

GlobMod<-read.table(paste0(results.dir.path,"\\",mod,"\\Tetrapod\\Global\\Model Performance.txt"))
files<-list.files(paste0(results.dir.path,"\\",mod), recursive=T)
files<-grep("Tetrapod",files,  value=T)
files<-grep("Model Performance",files,  value=T)
ModPer<-do.call(rbind, lapply(files, function(x){
  
  dat<-read.table(paste(results.dir.path,mod, x, sep="\\"))
  dat$region<-strsplit(x,"/")[[1]][2]
  return(dat)
}))

ModPer$region<-gsub("[.]", "- \n", ModPer$region)
ModPer$region<-gsub("Sibserian","Siberian", ModPer$region)
dim(ModPer);ModPer<-ModPer[ModPer$Rsq>0.25,]; dim(ModPer)
ModPer$region<-as.factor(ModPer$region)
ModPer$region<-factor(ModPer$region, levels=levels(ModPer$region)[c(7,1:6,8:20)])

ggplot()+
  geom_boxplot(data=ModPer, aes(x=region, y=Rsq))+
  scale_y_continuous(limits=c(0,1))+
  labs(y=  expression (~R^2))+
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title.x = element_blank())

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Figure 2- Variable importance measures from the global models of
## tetrapod, amphibian, reptile, bird and mammal species richness
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## 1. Collate the variable importance scores for the global models
files<-list.files(paste0(results.dir.path,"\\",mod), recursive=T)
files<-grep("bootstrap", files, value = T)
files<-grep("Global", files, value = T)
f<-grep("Tetrapod", files, value = T)

VarImp<-read.table(paste(results.dir.path,mod,f, sep="\\" ))
VarImp$VI<-sqrt((VarImp$MSErand-VarImp$MSEasis)/VarImp$MSEasis)
VariableNames<-as.data.frame(levels(VarImp$var))
colnames(VariableNames)[1]<-"var"
VariableNames$var<-sort(VariableNames$var)
VariableNames$names<-c("Annual Precipitation","Area Anthropogenic Land Use","Long Term Land Cover Change",
                       "Habitat Diversity", "Human Influence Index","Invasive Alien Species","Insularity",
                       "Area Protected Land","Mean Temperature","Long Term Climate Stability",
                       "Minimum Elevation","Precipation Seasonality","Std. Dev. Elevation",
                       "Short Term Land Cover Change",
                       "Temperature Seasonality","Total Species Richness")
VariableNames$class<-c("Environmental","Anthropogenic","Anthropogenic",
                       "Anthropogenic","Anthropogenic","Anthropogenic","Environmental",
                       "Anthropogenic","Environmental","Environmental",
                       "Environmental","Environmental","Environmental","Anthropogenic",
                       "Environmental","Other")
VariableNames$class2<-c("E","H","H","H","H","H","E","H","E","E","E","E","E","H","E","O")
VariableNames
VarImp<-merge(VarImp, VariableNames, by="var")
VarImp$names<-factor(VarImp$names, names(sort(tapply(VarImp$VI,VarImp$names, median, na.rm=T))))

mypal<-c(viridis_pal()(3), "grey")

VarImp$class<-as.factor(VarImp$class)
VarImp$class<-factor(VarImp$class, levels = c("Anthropogenic","Environmental", "Other"))
VarImp$names<-paste(VarImp$names,VarImp$class2, sep=" -")
GlobalVI<-VarImp
VarImp<-droplevels(VarImp)
VarImp$names<-base::factor(VarImp$names, names(sort(tapply(GlobalVI$VI,GlobalVI$names, median, na.rm=T))))

Ind<-ggplot()+
  geom_boxplot(data=VarImp, aes(x=names, y=VI), fill="grey")+ #, fill=class
  coord_flip()+
  labs(y="Variable Importance", x="", fill="Variable Class")+
  scale_y_continuous(limits=c(0, max(VarImp$VI)))+
  scale_fill_manual(values=mypal)+
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(vjust=-0.5, size=12),
        axis.title.y = element_text(vjust=1.5, size=12),
        axis.text = element_text(size=12),
        legend.key = element_blank(),
        legend.position ="none") 

taxa<-c("Amphibian","Reptile","Bird","Mammal")
taxaPlots<-lapply(taxa, function(t){
  
  VarImp<-read.table(paste0(results.dir.path,"\\",mod,"\\",t,"\\Global\\Variable Importance bootstrap.txt"))
  VarImp$VI<-sqrt((VarImp$MSErand-VarImp$MSEasis)/VarImp$MSEasis)
  VariableNames
  VarImp<-merge(VarImp, VariableNames, by="var")
  VarImp$names<-factor(VarImp$names, names(sort(tapply(VarImp$VI,VarImp$names, median, na.rm=T))))
  blank <- grid.rect(gp=gpar(col="white"))
  VarImp$class<-as.factor(VarImp$class)
  VarImp$class<-factor(VarImp$class, levels = c("Anthropogenic","Environmental", "Other"))
  VarImp<-droplevels(VarImp)
  VarImp$names<-paste(VarImp$names,VarImp$class2, sep=" -")
  VarImp$names<-base::factor(VarImp$names, names(sort(tapply(GlobalVI$VI,GlobalVI$names, median, na.rm=T))))
  
  img <- readPNG(paste0("X:\\Species rarity\\Silhouettes\\",t,"s.png"))
  g <- rasterGrob(img, interpolate=TRUE)
  
  Ind<-ggplot()+
    geom_boxplot(data=VarImp, aes(x=names, y=VI), fill="grey")+
    coord_flip()+
    labs(y="Variable Importance", x="", fill="Variable Class")+
    scale_y_continuous(limits=c(0, 8))+
    annotation_custom(g, ymin=5, ymax =8,
                      xmin=1, xmax=3)+
    theme(axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          # text = element_text(size=10),
          axis.title.x = element_text(vjust=-0.5, size=10),
          axis.title.y = element_text(vjust=1.5),
          axis.text = element_text(size=12),
          legend.key = element_blank(),
          legend.title.align=0.5,
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.position = "bottom")
  return(Ind)
})

# jpeg("X:\\Species rarity\\Rarity- Global analysis\\Plots\\Variable Importance Figure 2.jpg",
#      height=29, width=25.5, res=300, units="cm")
# grid.arrange(arrangeGrob(Ind+labs(title="a.")+theme(axis.title.x = element_blank()), 
#                          arrangeGrob(taxaPlots[[1]]+labs(title="b.")+theme(legend.position = "none", axis.title.x = element_blank()),
#                                      taxaPlots[[2]]+labs(title="c.")+theme(legend.position = "none", axis.text.y = element_blank(), 
#                                                                            axis.ticks.y = element_blank(), axis.title.x = element_blank()),
#                                      ncol=2, widths=c(3,1.8)),
#                          arrangeGrob(taxaPlots[[3]]+labs(title="d.")+theme(legend.position = "none", axis.title.x = element_blank()),
#                                      taxaPlots[[4]]+labs(title="e.")+theme(legend.position = "none", axis.text.y = element_blank(), 
#                                                                            axis.ticks.y = element_blank(), axis.title.x = element_blank()),
#                                      ncol=2, widths=c(3,1.8)),
#                          bottom=textGrob("       Variable Importance", gp=gpar(fontsize=12))))
# dev.off()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Figure 3: Heat map of individual variaable importance for the global and regional models
## of threatened tetrapod species richness
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

files<-list.files(paste0(results.dir.path,"\\",mod), recursive=T)
files<-grep("bootstrap", files, value = T)
files<-grep("Tetrapod", files, value = T)
files<-files[-grep("Global",files)]

VarImp<-do.call(rbind, lapply(files, function(f){
  region<-strsplit(f,"/")[[1]][[2]]
  dat<-read.table(paste(results.dir.path,mod,f, sep="\\"))
  
  ## want to get rid of the blocks with poor models
  perf<-read.table(paste(results.dir.path,mod,"Tetrapod",region,"Model Performance.txt",sep="\\"))
  perf<-perf[perf$Rsq>0.25,]
  dat<-dat[(dat$block %in% perf$block),]
  dat$VI<-sqrt((dat$MSErand-dat$MSEasis)/dat$MSEasis)
  dat$VI<-ifelse(is.na(dat$VI),0, dat$VI)
  dat<-merge(dat, VariableNames, by="var", all=T)
  VI<-as.data.frame(tapply(dat$VI, dat$var, median))
  VI$class<-rownames(VI); rownames(VI)<-NULL
  VI$region<-region
  colnames(VI)[1]<-"VI"
  return(VI)
}))

VarImp<-as.data.frame(VarImp)
head(VarImp)

## Summarise the data
Test<-ddply(VarImp, .(class,region),summarise, mean=mean(VI))
Test$rescale<-(Test$mean)
Test$rescale<-as.numeric(Test$rescale)
Test$region<-gsub("[.]","- \n", Test$region)
Test$region<-gsub("Sibserian","Siberian", Test$region)
## order the variables by the mean importance in the global model
head(GlobalVI)
Test$var<-Test$class
Test<-merge(Test[,c(2:5)], VariableNames, by.x="var",by.y="var")
Test$names<-paste(Test$names,Test$class2, sep=" -")
Test$names<-base::factor(Test$names, names(sort(tapply(GlobalVI$VI,GlobalVI$names, median, na.rm=T))))
Test$lab<-Test$names#paste(Test$names,Test$class2, sep=" - ")
odr<-as.data.frame(levels(as.factor(Test$names)))
odr$order<-seq(1,dim(odr)[1],1)
colnames(odr)[1]<-"names"
Test<-merge(Test, odr, by="names")
Test$lab<-factor(Test$lab, levels=c(levels(Test$lab),"Rsq"))

## want to put the model performance along the top
files<-list.files(paste0(results.dir.path, "\\",mod), recursive=T)
files<-grep("Model Performance",files, value=T)
files<-grep("Tetrapod",files, value=T)
files<-files[-grep("Global",files)]

Perf<-as.data.frame(do.call(rbind,lapply(files, function(f){
  region<-strsplit(f,"/")[[1]][2]
  dat<-read.table(paste(results.dir.path,mod,"Tetrapod",region,"Model Performance.txt",sep="\\"))
  dat<-dat[dat$Rsq>0.25,]
  return(c(region,mean(dat$Rsq)))})),stringsAsFactors = F)

colnames(Perf)<-c("region","rsq")
Perf$region<-gsub("[.]","- \n", Perf$region)
Perf$region<-gsub("Sibserian","Siberian", Perf$region)
Perf$y<-dim(odr)[1]+1
Perf$rsq<-round(as.numeric(Perf$rsq),2)
Perf<-Perf[order(Perf$rsq),]
Perf$order<-seq(1,nrow(Perf),1)
Perf$region<-reorder(Perf$region, Perf$order)

Test$region<-factor(Test$region,ordered=T, levels=rev(levels(Perf$region)))

##Plotting
Region<-ggplot(Test, aes(x=region, y=lab)) +
  geom_tile(aes(fill = rescale), colour = "white", width=0.9) +
  annotate("text",label=Perf$rsq, x=Perf$region, y=Perf$y)+
  scale_fill_viridis(option="B", direction=-1)+
  scale_y_discrete(drop=F)+  
  labs(fill="Variable \nImportance")+
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size=12),
        legend.key = element_blank(),
        legend.title.align=0.5,
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1, size=12))

## GLOBAL
files<-list.files(paste0(results.dir.path,"\\", mod), recursive=T)
files<-grep("Variable Importance",files, value=T)
files<-grep("Global",files, value=T)
files<-grep("Tetrapod",files, value=T)

VarImp<-do.call(rbind, lapply(files, function(f){
  region<-strsplit(f,"/")[[1]][2]
  dat<-read.table(paste(results.dir.path,mod,"Tetrapod",region,"Variable Importance bootstrap.txt",sep="\\"))
  
  dat$VI<-sqrt((dat$MSErand-dat$MSEasis)/dat$MSEasis)
  dat<-merge(dat, VariableNames, by="var", all=T)
  VI<-as.data.frame(tapply(dat$VI, dat$var, median))
  VI$class<-rownames(VI); rownames(VI)<-NULL
  VI$region<-region
  colnames(VI)[1]<-"VI"
  return(VI)
}))

VarImp<-as.data.frame(VarImp)
head(VarImp)

## Summarise the data
Test<-ddply(VarImp, .(class,region),summarise, mean=mean(VI))
Test$rescale<-(Test$mean)
Test$rescale<-as.numeric(Test$rescale)
Test$region<-gsub("[.]","- \n", Test$region)
Test$region<-gsub("Sibserian","Siberian", Test$region)
## order the variables by the mean importance in the global model
head(GlobalVI)
Test$var<-Test$class
Test<-merge(Test[,c(2:5)], VariableNames, by.x="var",by.y="var")
Test$names<-paste(Test$names,Test$class2, sep=" -")
Test$names<-base::factor(Test$names, names(sort(tapply(GlobalVI$VI,GlobalVI$names, median, na.rm=T))))
Test$lab<-Test$names#paste(Test$names,Test$class2, sep=" - ")
odr<-as.data.frame(levels(as.factor(Test$names)))
odr$order<-seq(1,dim(odr)[1],1)
colnames(odr)[1]<-"names"
Test<-merge(Test, odr, by="names")
Test$region<-paste0("           ",Test$region)
Test$lab<-factor(Test$lab, levels=c(levels(Test$lab),"R-squared"))

dat<-read.table(paste0(results.dir.path,"\\",mod,"\\Tetrapod\\Global\\Model Performance.txt"))
R<-round(mean(dat$Rsq),2)
n<-67340

##Plotting
Global<-ggplot(Test, aes(x=region, y=lab)) +
  geom_tile(aes(fill = rescale), colour = "white", width=0.9) +
  annotate("text",label=R, x=1, y=17)+
  scale_fill_viridis(option="B", direction=-1, limits=)+
  scale_y_discrete(drop=F)+
  labs(fill="Variable \nImportance")+
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size=12),
        legend.key = element_blank(),
        legend.title.align=0.5,
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1, size=12))


# jpeg(paste0(getwd(),"\\Plots\\Heat map Figure 3.jpg"),
#      height=20, width=40, res=300, units="cm")
# 
# grid.arrange(Global+theme(legend.position = "none"),
#              Region+theme(axis.text.y = element_blank(),
#                           axis.ticks.y =element_blank())
#              , nrow=1, widths=c(1.3,5))
# dev.off()

