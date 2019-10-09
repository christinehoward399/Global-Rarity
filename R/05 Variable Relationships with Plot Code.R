##~~~~~~~~~~~~~~~~~~~~~ ##

## Variable Relationships

##~~~~~~~~~~~~~~~~~~~~~ ##

rm(list=ls())
.libPaths("X:\\Old N drive\\R\\R-2.14.0\\library")
setwd("X:\\Species rarity\\Rarity- Global analysis\\")

library(ggplot2)
library(gridExtra)
library(grid)
library(randomForest)
library(sp)
library(raster)
library(viridis)
library(png)

## Taxonomic group
taxa<-c("Tetrapod","Amphibian","Reptile","Bird","Mammal")

## Results path
results.dir.path<-paste0(getwd(),"\\Model Output\\Total Species Richness")

## Read in data
my.data<-read.csv(paste0(getwd(),"\\Data\\Combined Gridded Data.csv"))

## Environmental and Anthropogenic Variables
clim.var<-c("Annual.Precip","Temp.seasonality","Precip.seasonality","Mean.Temp","MinElev",
            "meanED","HabDiv","Islandness","STDElev", "HII", "Anthro","GCCchange", "IAS","IUCN","STLCchange")

## All Explanatory variables
CoVars<-c(clim.var)

## Response variable- threatened species richness, total species richness
Response<-"AtRisk"

## Does total species include DD species
TOT<-c("AtRisk","Secure")

## Species included in At-risk
RISK<-c("AtRisk")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## LOOP 
## Bit awkward but i've only magaed to get the partialPlot function to run in sequence
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

Files<-list.files(results.dir.path, recursive=T)
Files<-grep("model.output.RF",Files,  value=T)

for ( x in  1:length(Files)){
  f<-Files[x]
  
  t<-strsplit(f,"/")[[1]][1]
  r<-strsplit(f,"/")[[1]][2]
  
  ## Format species data
  MyDat<-my.data[,c("x","y",clim.var,r,grep(t, colnames(my.data), value=T))]
  colnames(MyDat)[grep("AtRisk",colnames(MyDat))]<-c("AtRisk")
  colnames(MyDat)[grep("Secure",colnames(MyDat))]<-c("Secure")
  colnames(MyDat)[grep("DD",colnames(MyDat))]<-c("DD")
  MyDat$Total<-apply(MyDat[,colnames(MyDat) %in% TOT],1,sum)
  
  ## Subset down to the realm of interest
  MyDat<-MyDat[complete.cases(MyDat[,r]),]
  
  ## output directory
  dir.create(paste(results.dir.path,t,r,"Partial residuals", sep="\\"))
  
  doneVars<-list.files(paste(results.dir.path,t,r,"Partial residuals", sep="\\"))
  
  if(length(doneVars)<length(CoVars)){
    
    ## Load up models
    load(paste(results.dir.path,f, sep="/"))
    
    ## Variables
    clim.var1<-CoVars

    for (v in 1:length(clim.var1)){
      
      var<-clim.var1[v]
      
      ## calculate the partial residuals
      modresults<-lapply(1:length(mod), function(i){
        dat<-partialPlot(mod[[i]]$mod, x.var=as.character(clim.var1[v]),
                         pred.data=MyDat, add=F, plot=F)
        dat<-as.data.frame(cbind(dat$x,dat$y))
        colnames(dat)<-c("x", paste0("mod",i))
        return(dat)
      })
      
      ## Format data for plotting
      ReducedDAT<-Reduce(function(...) merge(..., all=T), modresults)
      ReducedDAT$m<-apply(ReducedDAT[,grep("mod", colnames(ReducedDAT), value=T)],1,mean)
      ReducedDAT$sd<-apply(ReducedDAT[,grep("mod", colnames(ReducedDAT), value=T)],1,sd) 
      write.table(ReducedDAT, paste0(paste(results.dir.path,t,r,"Partial residuals", sep="\\"),"\\",var,"_partials.txt"))
      rm(ReducedDAT); rm(modresults)
    }
  }
  print(f)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## PLOTS
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
vars<-as.data.frame(clim.var)
vars$Name<-c("Annual precipitation (mm)", 
             "Temperature seasonality",
             "Precipitation seasonality",
             "Mean temperature",
             "Minimum elevation (m)",
             "Long term climate stability",
             "Habitat diversity",
             "Insularity",
             "Elevation std (m)",
             "Human Influence Index",  
             "Area of anthropogenic land use",
             "Long term land use change",
             "Invasive alien species",
             "Area of protected land",
             "Short term land use change")
colnames(vars)[1]<-"clim.var"

mod<-"Threatened Species Richness"

## Results path
results.dir.path<-paste0(getwd(),"\\Model Output\\")
test<-lapply(clim.var, function(v){  
  
  plotData<-lapply(taxa, function(t){
    
    # Read in data
    ReducedDAT<-read.table(paste0(results.dir.path, mod,"\\",t, "\\Global\\Partial residuals\\",v,"_partials.txt"))
    for ( i in grep("mod", colnames(ReducedDAT))){ReducedDAT[,i]<-scale(ReducedDAT[,i])}
    ReducedDAT$m<-apply(ReducedDAT[,grep("mod", colnames(ReducedDAT))], 1,mean)
    ReducedDAT$sd<-apply(ReducedDAT[,grep("mod", colnames(ReducedDAT))], 1,sd)
    ReducedDAT<-ReducedDAT[ReducedDAT$x<up,]
    ReducedDAT<-ReducedDAT[ReducedDAT$x>lower,]
    ReducedDAT$upperSD<-ReducedDAT$m+ReducedDAT$sd
    ReducedDAT$lowerSD<-ReducedDAT$m-ReducedDAT$sd
    
    nbcst<-append(ReducedDAT[,"x"], rev(ReducedDAT[,"x"]))
    nbstanders<-append(ReducedDAT$upperSD,rev(ReducedDAT$lowerSD))
    nb<-as.data.frame(cbind(nbcst, nbstanders), strings.as.factors=F)
    nb$taxa<-t
    
    ReducedDAT<-ReducedDAT[,c("x","m")] 
    ReducedDAT$taxa<-t
    return(list(ReducedDAT, nb))
  })
  
  ReducedDAT<-sapply(plotData, "[",1)
  ReducedDAT<-do.call("rbind", ReducedDAT)
  ReducedDAT$taxa<-as.factor(ReducedDAT$taxa)
  
  nb<-sapply(plotData, "[",2)
  nb <- do.call("rbind", nb) 
  nb$taxa<-as.factor(nb$taxa)
  
  number_ticks <- function(n) {function(limits) pretty(limits, n)}
  
  xlab<-vars[vars$clim.var==v,"Name"]
  
  p<-ggplot()+
    geom_polygon(data=nb,aes(x=nbcst, y=nbstanders, fill=taxa), alpha=0.5)+
    geom_line(data=ReducedDAT, aes_string(x=ReducedDAT[,"x"], y=ReducedDAT[,"m"], col=ReducedDAT[,"taxa"]))+
    labs(x=xlab, y="", fill="Taxa")+
    scale_x_continuous(breaks=number_ticks(4))+
    scale_y_continuous(breaks=number_ticks(4))+
    scale_color_discrete(guide=F)+
    guides(fill=guide_legend(ncol=2))+
    theme(
      panel.background=element_blank(),
      panel.grid.minor=element_blank(),
      axis.line.x=element_line(colour="black"),
      axis.line.y=element_line(colour="black"),
      axis.text=element_text(colour="black",size=12),
      axis.title=element_text(colour="black",size=12),
      #axis.ticks.x=element_line(colour="black"),
      legend.text=element_text(size=12),
      legend.title=element_text(size=12)
    )
  
  return(p)
})

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

legend <- g_legend(test[[15]])

xlab <- "Mean Temperature (°C)"

p3<-arrangeGrob(grid.arrange(test[[1]]+theme(legend.position="none")+labs(title="a."), 
                             test[[2]]+theme(legend.position="none")+labs(title="b."),
                             test[[3]]+theme(legend.position="none")+labs(title="c."),
                             test[[4]]+xlab(xlab)+theme(legend.position="none")+labs(title="d."),
                             test[[5]]+theme(legend.position="none")+labs(title="e."),
                             test[[6]]+theme(legend.position="none")+labs(title="f."),
                             test[[7]]+theme(legend.position="none")+labs(title="g."),
                             test[[8]]+xlab(expression("Insularity "~(km^2)))+theme(legend.position="none")+labs(title="h."),
                             test[[9]]+theme(legend.position="none")+labs(title="i."),
                             test[[10]]+theme(legend.position="none")+labs(title="j."),
                             test[[11]]+xlab(expression("Area of anthropogneic land use "~(km^2)))+theme(legend.position="none")+labs(title="k."),
                             test[[12]]+xlab(expression("Long term land use change (% change)"))+theme(legend.position="none")+labs(title="l."),ncol=4),
                grid.arrange(test[[13]]+theme(legend.position="none")+labs(title="m."),
                             test[[14]]+xlab(expression("Area of protected land "~(km^2)))+theme(legend.position="none")+labs(title="n."),
                             test[[15]]+xlab(expression("Short term land use change (% change)"))+theme(legend.position="none")+labs(title="o."),
                             legend, ncol=4),
                nrow=2, heights=c(3,1), 
                left=textGrob("Standardised total species richness", 
                              gp=gpar(fontsize=13), rot=90))

grid.arrange(p3)

# jpeg(file=paste0("X:\\Species rarity\\Rarity- Global analysis\\Plots\\ Supp figure partials ",mod,".jpg"),
#      height=20, width=38, res=300, units="cm")
# grid.arrange(p3)
# dev.off()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Figure of partials for poster and ms
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
mypal<-viridis_pal(option = "inferno")(5)
mypal<-mypal[c(2,4)]

results.dir.path<-paste0(getwd(),"\\Model Output\\")
taxa<-c("Amphibian","Reptile","Bird","Mammal")

fig4<-lapply(taxa, function(t){
  
  clim.var<-ifelse(c(t=="Amphibian",t=="Amphibian"),c("STDElev","GCCchange"),
                   ifelse(c(t=="Reptile",t=="Reptile"),c("Mean.Temp","GCCchange"),
                          ifelse(c(t=="Bird",t=="Bird"),c("Temp.seasonality","Anthro"),
                                 ifelse(c(t=="Mammal",t=="Mammal"),c("Islandness","Anthro"),clim.var))))
  
  test<-lapply(clim.var, function(v){  
    print(v)
    
    ## Threatened species richness
    # Read in data
    ReducedDAT<-read.table(paste0(results.dir.path, "Threatened Species Richness\\",t, "\\Global\\Partial residuals\\",v,"_partials.txt"))
    for ( i in grep("mod", colnames(ReducedDAT))){ReducedDAT[,i]<-scale(ReducedDAT[,i])}
    ReducedDAT$m<-apply(ReducedDAT[,grep("mod", colnames(ReducedDAT))], 1,mean)
    ReducedDAT$sd<-apply(ReducedDAT[,grep("mod", colnames(ReducedDAT))], 1,sd)
    ReducedDAT<-ReducedDAT[ReducedDAT$x<up,]
    ReducedDAT<-ReducedDAT[ReducedDAT$x>lower,]
    ReducedDAT$upperSD<-ReducedDAT$m+ReducedDAT$sd
    ReducedDAT$lowerSD<-ReducedDAT$m-ReducedDAT$sd
    
    nbcst<-append(ReducedDAT[,"x"], rev(ReducedDAT[,"x"]))
    nbstanders<-append(ReducedDAT$upperSD,rev(ReducedDAT$lowerSD))
    nb<-as.data.frame(cbind(nbcst, nbstanders), strings.as.factors=F)
    nb$RSP<-"Threatened Species Richness"
    
    number_ticks <- function(n) {function(limits) pretty(limits, n)}
    xlab<-vars[vars$clim.var==v,"Name"]
    tempLab<-"Mean Temperature (?C)"
    xlab<-ifelse(v=="Anthro",(expression("Area of anthropogneic land use "~(km^2))),
                 ifelse (v=="GCCchange",(expression("Long term land use change (% change)")),
                         ifelse(v=="Islandness",(expression("Insularity - land mass area"~(km^2))),
                                ifelse(v=="Mean.Temp",tempLab,xlab))))
    
    p<-ggplot()+
      geom_polygon(data=nb,aes(x=nbcst, y=nbstanders, fill=RSP), alpha=0.5)+
      geom_line(data=ReducedDAT, aes_string(x=ReducedDAT[,"x"], y=ReducedDAT[,"m"]), col=mypal[1])+
      scale_x_continuous(breaks=number_ticks(4))+
      scale_color_manual(values=mypal)+
      scale_fill_manual(labels = c("Threatened Species Richness","Total Species Richness"), 
                        values=alpha(mypal, 0.5))+
      guides(fill=guide_legend(nrow=2),color = FALSE)+
      labs(fill="Response", y="Threatened")+
      xlab(xlab)+
      theme(
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"),
        axis.text=element_text(colour="black",size=12),
        axis.title=element_text(colour="black",size=12),
        #axis.ticks.x=element_line(colour="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.position="bottom"
      )
    
    
    ## Total species richness
    # Read in data
    ReducedDAT1<-read.table(paste0(results.dir.path, "Total Species Richness\\",t, "\\Global\\Partial residuals\\",v,"_partials.txt"))
    ReducedDAT1$m<-apply(ReducedDAT1[,grep("mod", colnames(ReducedDAT1))], 1,mean)
    xs<-scale(ReducedDAT1$m)
    for ( i in grep("mod", colnames(ReducedDAT1))){ReducedDAT1[,i]<-scale(ReducedDAT1[,i])}
    ReducedDAT1$m<-apply(ReducedDAT1[,grep("mod", colnames(ReducedDAT1))], 1,mean)
    ReducedDAT1$sd<-apply(ReducedDAT1[,grep("mod", colnames(ReducedDAT1))], 1,sd)
    
    # Subset to the middle 90% of data for each varible
    lower<-unname(quantile(ReducedDAT1$x,  probs = c(5)/100))
    up<-unname(quantile(ReducedDAT1$x,  probs = c(95)/100))
    ReducedDAT1<-ReducedDAT1[ReducedDAT1$x<up,]
    ReducedDAT1<-ReducedDAT1[ReducedDAT1$x>lower,]
    ReducedDAT1$upperSD<-ReducedDAT1$m+ReducedDAT1$sd
    ReducedDAT1$lowerSD<-ReducedDAT1$m-ReducedDAT1$sd
    
    nbcst1<-append(ReducedDAT1[,"x"], rev(ReducedDAT1[,"x"]))
    nbstanders1<-append(ReducedDAT1$upperSD,rev(ReducedDAT1$lowerSD))
    nb1<-as.data.frame(cbind(nbcst1, nbstanders1), strings.as.factors=F)
    nb1$RSP<-"Total Species Richness"
    
    number_ticks <- function(n) {function(limits) pretty(limits, n)}
    
    unscaled_vals <- xs + attr(xs, 'scaled:scale') + attr(xs, 'scaled:center')
    
    p<-p+geom_polygon(data=nb1,aes(x=nbcst1, y=nbstanders1, fill=RSP), alpha=0.5)+
      geom_line(data=ReducedDAT1, aes_string(x=ReducedDAT1[,"x"], y=ReducedDAT1[,"m"]), col=mypal[2])+
      scale_y_continuous(sec.axis = sec_axis(~. + attr(xs, 'scaled:scale') + attr(xs, 'scaled:center'), name = "Total"))
    
    return(p)
    print(v)
  })
  
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
  }
  
  legend <- g_legend(test[[1]])
  blank <- grid.rect(gp=gpar(col="white"))
  
  img <- readPNG(paste0("X:\\Species rarity\\Silhouettes\\",t,"s.png"))
  g <- rasterGrob(img, interpolate=TRUE)
  
  TitleLab<-ifelse(c(t=="Amphibian",t=="Amphibian"),c("a.","b."),
                   ifelse(c(t=="Reptile",t=="Reptile"),c("c.","d."),
                          ifelse(c(t=="Bird",t=="Bird"),c("e.","f."),
                                 ifelse(c(t=="Mammal",t=="Mammal"),c("g.","h."),clim.var))))
  
  p3<-arrangeGrob(grid.arrange(test[[1]]+labs(title= TitleLab[1])+theme(legend.position="none"), 
                               test[[2]]+labs(title= TitleLab[2])+theme(legend.position="none"),ncol=2),
                  grid.arrange(blank,g,ncol=1, heights=c(0.8,1)),ncol=2, widths=c(2,0.2))

  return(p3)
  
})

# jpeg(file="X:\\Species rarity\\Rarity- Global analysis\\Plots\\Figure 4 partials.jpg",
#      width=20,height=22, units="cm", res=300)
# grid.arrange(arrangeGrob(fig4[[1]],fig4[[2]],fig4[[3]],fig4[[4]],legend, ncol=1, heights=c(1,1,1,1,0.3)),
#              left=textGrob("Standardised species richness", gp=gpar(fontsize=13), rot=90))
# dev.off()
