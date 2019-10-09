## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Extract Model performance and variable importance scores

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

rm(list=ls())

##Libraries
.libPaths("N:\\R\\R-2.14.0\\library")
library(snowfall)
library(randomForest)

setwd("X:\\Species rarity\\Rarity- Global analysis\\")

#Taxonomic group
taxa<-c("Amphibian","Reptile","Bird","Mammal","Tetrapod")

## Numer of CPUS to run on
CP<-4

## Excluded Species
excluded.species<-"x"

## Results path
results.dir.path<-"X:\\Species rarity\\Rarity- Global analysis\\Model Output"

## Read in data
my.data<-read.csv(paste0(getwd(),"\\Data\\Combined Gridded Data.csv"))

## Environmental and Anthropogenic Variables
clim.var<-c("Annual.Precip","Temp.seasonality","Precip.seasonality","Mean.Temp","MinElev",
            "meanED","HabDiv","Islandness","STDElev", "HII", "Anthro","GCCchange", "IAS","IUCN","STLCchange")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##        MODEL PROCESSING
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
outpath<-results.dir.path
files<-list.files(outpath, recursive=T)
files<-grep(".model.output",files, value=T)
files<-files[c(261:267)]#20,202,239)]

sfInit(parallel=TRUE, cpus=3)
sfLibrary(randomForest)
sfExport(list=c("files","my.data","outpath","clim.var"))
sfLapply(files, function(f){
  
  ## Format species data
  t<-strsplit(f, "/")[[1]][2]
  z<-strsplit(f, "/")[[1]][3]
  modType<-strsplit(f, "/")[[1]][1]
  
  if(file.exists(paste(outpath,paste(outpath, paste(replace(unlist(strsplit(f,"/")), length(unlist(strsplit(f,"/"))), "Model Performance.txt"), collapse="/"), sep="\\")))==FALSE){
    
    ##~~~~~~~~  Load Model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    load(paste(outpath, f, sep="\\"))
    
    ##~~~~~~~~   Model performance and summary statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    Rsq<-sapply(1:length(mod), function(i){Rsq<-mod[[i]]$Rsq})
    block<-seq(1,length(mod),1)
    perVarexp<-sapply(1:length(mod), function(i){perVarexp<-(mod[[i]]$mod$rsq[length(mod[[i]]$mod$rsq)])*100})
    mse<-sapply(1:length(mod), function(i){mse<-(mod[[i]]$mod$mse[length(mod[[i]]$mod$mse)])})
    nt<-sapply(1:length(mod), function(i){nt<-(mod[[i]]$mod$ntree)})
    mtr<-sapply(1:length(mod), function(i){nt<-(mod[[i]]$mod$mtry )})
    IntRsq<-sapply(1:length(mod), function(i){IntRsq<-mean(mod[[i]]$mod$rsq)})
    dat<-cbind(block, Rsq, perVarexp, mse, nt, mtr, IntRsq)
    write.table(dat,paste(outpath, paste(replace(unlist(strsplit(f,"/")), length(unlist(strsplit(f,"/"))), "Model Performance.txt"), collapse="/"), sep="\\")) 
    
    ##~~~~~~~~~~~~~~~Extract internal Variable Importance metrics for comparison ~~~~~~~~~~~~ ##
    VarImp<-lapply(1:length(mod), function(i){
      VI<-importance(mod[[i]]$mod, scale=F)
      VI<-as.data.frame(VI)
      VI$var<-row.names(VI)
      VI$block<-i
      return(VI)})
    VarImp<-Reduce(function(...) merge(..., all=T), VarImp)
    
    ##~~~~~~~~  Calculate variable importance using equation with bootstrapping ~~~~~~~~~~~~ ##
    
    CoVars<-unique(VarImp$var)
    
    rm(TOT); rm(RISK)
    ## Does total species include DD species
    if(modType=="All DD as secure"){TOT<- c("AtRisk","Secure","DD")}else {TOT<-c("AtRisk","Secure")}
    if(modType=="All DD as threatened"){TOT<- c("AtRisk","Secure","DD")}else {TOT<-TOT}
    
    ## Species included in At-risk - i.e.DD species?
    if(modType=="All DD as threatened"){RISK<- c("AtRisk","DD")}else {RISK<-c("AtRisk")}
    
    ## Format species data
    MyDat<-my.data[,c("x","y",clim.var,z,grep(t, colnames(my.data), value=T))]
    colnames(MyDat)[grep("AtRisk",colnames(MyDat))]<-c("AtRisk")
    colnames(MyDat)[grep("Secure",colnames(MyDat))]<-c("Secure")
    colnames(MyDat)[grep("DD",colnames(MyDat))]<-c("DD")
    MyDat$Total<-apply(MyDat[,colnames(MyDat) %in% TOT],1,sum)
    if(modType=="All DD as threatened"){MyDat$AtRisk<-apply(MyDat[,colnames(MyDat) %in% RISK],1,sum)}
    if(modType=="50% DD as threatened"){MyDat$AtRisk<-round(MyDat$AtRisk+(MyDat$DD*0.5))}
    
    ## Subset down to the realm of interest
    MyDat<-MyDat[complete.cases(MyDat[,z]),]
    
    ## Manual Variable Importance
    VarImpMan<-lapply(1:length(mod), function(i){
      VIman<-do.call(rbind,lapply(CoVars, function(x){
        VImanrep<-do.call(rbind, lapply(1:100, function(v){
          dat<-MyDat#[MyDat[,z]==1,]
          dat$PREDasis<-predict(mod[[i]]$mod, newdata=dat)
          dat$rand<-sample(dat[,x])
          dat[,x]<-dat$rand
          dat$PREDrand<-predict(mod[[i]]$mod, newdata=dat)
          MSEasis<-mean( (dat$PREDasis - dat$AtRisk)^2, na.rm = TRUE)
          MSErand<-mean( (dat$PREDrand - dat$AtRisk)^2, na.rm = TRUE)
          manVI<-sqrt((MSErand-MSEasis)/MSEasis)
          if(is.na(manVI)){manVI<-0} ## i.e. where mse using a randomised variable is < mse from variable as is. Variable Importance = 0
          return(c(x,MSEasis,MSErand, manVI))
        }))
        return(VImanrep)
      }))
      VIman<-as.data.frame(VIman, stringsAsFactors=F)
      colnames(VIman)<-c("var","MSEasis","MSErand","ManVI")
      VIman$block<-i
      return(VIman)
    })
    VarImpMan<-Reduce(function(...) merge(..., all=T), VarImpMan)
    colnames(VarImpMan)<-c("var","MSEasis","MSErand","ManVI","block")
    VI<-merge(VarImp, VarImpMan, by=c("var","block"))
    write.table(VI,paste(outpath, paste(replace(unlist(strsplit(f,"/")), length(unlist(strsplit(f,"/"))), "Variable Importance bootstrap.txt"), collapse="/"), sep="\\")) 
    
  }
})