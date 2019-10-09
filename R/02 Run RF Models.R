## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Run Random Forest Models

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

rm(list=ls())

## ~~~~~~~~~~~~~~~ Library etc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

.libPaths("N:\\R\\R-2.14.0\\library")

library(randomForest)
library(snowfall)

## ~~~~~~~~~~~~~~~ Set parameters, load data, etc~~~~~~~~~~~~ ##

##Set working directory
setwd("X:\\Species rarity\\Rarity- Global analysis\\")

## Taxonomic group
taxaList<-c("Tetrapod")##,

## Numer of CPUS to run on
CP<-4

## Excluded Species
excluded.species<-"x"

## Results path
results.dir.path<-paste0(getwd(),"\\Model Output\\Total Species Richness")#With DD as covariate")#)
#All DD as secure")#Recent land cover change")

## Create overall model output director
if(file.exists(results.dir.path) == FALSE){dir.create(results.dir.path)}

## Read in data
my.data<-read.csv(paste0(getwd(),"\\Data\\Combined Gridded Data.csv"))

## Environmental and Anthropogenic Variables
clim.var<-c("Annual.Precip","Temp.seasonality","Precip.seasonality","Mean.Temp","MinElev",
            "meanED","HabDiv","Islandness","STDElev", "HII", "Anthro","GCCchange", "IAS","IUCN","STLCchange")

## All Explanatory variables
CoVars<-c(clim.var)#, "Total")#,"DD")

## Zoogeographic Regions
regionList<-c("African","Amazonian","Arctico.Sibserian",
           "Australian","Chinese","Eurasian","Guinea.Congolian","Indo.Malayan","Japanese","Madagascan","Mexican","North.American",
           "Novozelandic","Oriental","Panamanian","Papua.Melanesian","Saharo.Arabian","South.American","Tibetan","Global")

## Response variable- threatened species richness, residual threatened species richness, DD species?? 
Response<-"AtRisk"

## Does total species include DD species
TOT<-c("AtRisk","Secure")

## Species included in At-risk
RISK<-c("AtRisk")

## Source functions
source("R/RF_function.R")
source("R/Output_function.R")

## ~~~~~~~~~~~~~~~ MODEL~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Working through each taxonomic group
for (t in taxaList){
  print(t)
  
  for (f in regionList){
    print(f)
    
    ## Format species data
    MyDat<-my.data[,c("x","y",clim.var,f,grep(t, colnames(my.data), value=T))]
    colnames(MyDat)[grep("AtRisk",colnames(MyDat))]<-c("AtRisk")
    colnames(MyDat)[grep("Secure",colnames(MyDat))]<-c("Secure")
    colnames(MyDat)[grep("DD",colnames(MyDat))]<-c("DD")
    MyDat$Total<-apply(MyDat[,colnames(MyDat) %in% TOT],1,sum)
    MyDat$AtRisk<-apply(MyDat[,colnames(MyDat) %in% RISK],1,sum)
    
    ## Subset down to the realm of interest
    MyDat<-MyDat[complete.cases(MyDat[,f]),]
    
    ##Set up output directory
    Output_files(outDir=results.dir.path,realm=f, taxa=t)
    
    #Summary data
    ncells.pres <- c("T.cells",length(subset(MyDat,AtRisk>=1)[,"AtRisk"])) #Number of cells with species present
    block.sum <- aggregate(AtRisk~get(f),data=MyDat, function(c)sum(c!=0)); 
    no.block.zero <- length(subset(block.sum,AtRisk!=0)[,"AtRisk"])#Number of blocks contain species presences
    write.table(rbind(block.sum,ncells.pres),paste(results.dir.path,"\\",t,"\\", f,"\\",t,".",f,".summary.stats.txt",sep=""),sep="",quote=FALSE,row.names=FALSE)
    
    block.include <- subset(block.sum,AtRisk > 0)[,1]
    
    if(as.numeric(ncells.pres[2]) >= 5 & no.block.zero > 1){  
      
      block.include <- subset(block.sum,AtRisk > 0)[,1]
      
      if(file.exists(paste(results.dir.path,t,f,paste(t,f,"model.output.RF.rda",sep="."),sep="/"))==FALSE){
        
        RF_funct(source.data=MyDat,
                 outDir=results.dir.path,
                 realm=f,
                 taxa=t,
                 blocks=block.include,
                 Response=Response)
        
      }else{excluded.species <- c(excluded.species,t)}
      
    }
  }
}