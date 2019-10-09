## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Blocking Code

## Block the globe and each region

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

rm(list=ls())
setwd("X:\\Species rarity\\Rarity- Global analysis")
.libPaths("N:\\R\\R-2.14.0\\library")

#Libraries
library(raster)
library(rgdal)
library(maptools)
library(plyr)
library(blockTools)

## Import in the climate data
load(paste0(getwd(),"\\Data\\Explanatory layers\\Explanatory variable rasters.Rdata"))

my.data<-lapply(1:length(my.var), function(x){
  dat<-cbind(coordinates(my.var[[x]]), getValues(my.var[[x]]))
  colnames(dat)<-c("x","y",names(my.var[[x]]))
  return(dat)
})

my.data<- Reduce(function(...) merge(..., all=T), my.data)
colnames(my.data)[3:8]<-c("Realms","Mean.Temp","Temp.seasonality","Annual.Precip","Precip.seasonality","HII")
colnames(my.data)[28]<-"Regions"
my.data$Mean.Temp<-my.data$Mean.Temp/10 ## temperature variables come multiplied by 10
realms<-data.frame(id=seq(1:11), realm=c("Afrotropical","Australia","Madagasca","Nearctic","Neotropical",
                                         "Oceania","Oriental","Palearctic","Panamanian","Saharo-Arabian","Sino-Japanese"))

regions<-data.frame(id=seq(1:20),region=c("African","Amazonian","Arctico-Siberian",
                                          "Australian","Chinese","Eurasian",
                                          "Guinea-Congolian","Indo-Malayan","Japanese",
                                          "Madagascan","Mexican","North American",
                                          "Novozelandic","Oriental","Panamanian",
                                          "Papua-Melanesian","Polynesian","Saharo-Arabian",
                                          "South American","Tibetan"))

my.data<-merge(my.data, realms, by.y="id", by.x="Realms")
my.data<-merge(my.data, regions, by.y="id", by.x="Regions")

## Going to block by the climatic variables
vars<-c("Mean.Temp","Temp.seasonality","Annual.Precip","Precip.seasonality")

## Grid 
all.range<-my.var[[1]]

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## grid up ecoregion data
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

#Import ecoregion data
i<-"F:\\Data\\Terrestrial Ecoregions\\official\\wwf_terr_ecos.shp"
ecoReg <- readShapePoly(i, proj4string=CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")) #Set to correct projection (i.e. the polygon projection it with output with)
ecoReg<- spTransform(ecoReg, proj4string(all.range))
ecoReg.num <- ecoReg[9] #ecoregion id
ecoReg.area <- ecoReg[2] #ecoregion area
# plot(ecoReg)
r.grid.l<-raster(extent(all.range),nrows=360, ncols=720, crs=proj4string(all.range))
vals <- 1:ncell(r.grid.l) #create vector of numbers 
r.grid.l <- setValues(r.grid.l, vals) #fill grid squares with numerical value to create label
r.grid.s <- disaggregate(r.grid.l, fact=c(10,10),method="") #disaggregate so that each 2.5' cell has number associating it to a 1 degree cell

##Create grid for dividing large ecoregions
r.grid.20 <- raster(extent(all.range),nrows=36, ncols=72, crs=proj4string(all.range))
vals <- 1:ncell(r.grid.20) #create vector of numbers 
r.grid.20.l <- setValues(r.grid.20, vals) #fill grid squares with numerical value to create label
r.grid.20.s <- disaggregate(r.grid.20.l, fact=c(10,10),method="") #disaggregate so that each 2.5' cell has number associating it to a 1 degree cell
r.grid.20.s <-crop(r.grid.20.s,extent(all.range))
r.grid.20.s

#Rasterize and convert into sampling units with max size of 10*res
raster.feature <- rasterize(ecoReg.num, r.grid.s)
raster.feature.l <- aggregate(raster.feature, 10, modal, progress='text') #if it overlaps with ecor 1, make it ecoregion 1
raster.ecoReg.num <- rasterize(ecoReg.num, r.grid.s,field="ECO_SYM")
raster.ecoReg.num.l <- aggregate(raster.ecoReg.num, 10, modal, progress='text')
raster.ecoReg.area <- rasterize(ecoReg.area, r.grid.s,field="AREA")
raster.ecoReg.area.l <- aggregate(raster.ecoReg.area, 10, modal, progress='text')
r.feat.reg <- stack(raster.feature.l,raster.ecoReg.num.l,r.grid.20.s) #each 0.5 accociated with ecoregion and larger grid cell
reg.id <- function(x,na.rm){as.numeric(paste(x[2],x[3],sep="."))}
sample.unit.id <- stackApply(r.feat.reg,c(1,1,1),fun=reg.id) #dif block dif ids (1 layer of id)

#Convert into data frame and merge with climate data
coord <- coordinates(sample.unit.id)
sample.id <- getValues(sample.unit.id)
sample.area <- getValues(raster.ecoReg.area.l)
sample.unit.df <- as.data.frame(cbind(coord,sample.id,sample.area))
sample.unit.df <- na.omit(sample.unit.df)
names(sample.unit.df)
plot(y~x, data=sample.unit.df,cex=0.1)

## merge in the climate data
sample.unit.df <- merge(sample.unit.df,my.data,by=c("x","y"))#,all.y=TRUE) 

regions<-c("African","Amazonian","Arctico-Siberian",
           "Australian","Chinese","Eurasian","Guinea-Congolian","Indo-Malayan","Japanese",
           "Madagascan","Mexican","North American","Novozelandic","Oriental","Panamanian",
           "Papua-Melanesian","Saharo-Arabian","South American","Tibetan", "Global")

## Now work through each realm and block and create the blocks
for (r in regions){
  
  # Subset down to realm
  if(r !="Global"){ sample.unit.realm<-sample.unit.df[sample.unit.df$region==r,]}else{sample.unit.realm<-sample.unit.df}
  ##Create sample units of only a certain size
  max.size <- 100000
  sample.id <- function(x){if(x["sample.area"] < max.size){strsplit(as.character(x["sample.id"]),".",fixed=TRUE)[[1]][1]}else{x["sample.id"]}} #if ecoregion larger-split it up
  sample.units.id <- apply(sample.unit.realm,1,sample.id)
  sample.units.id <- cbind(sample.unit.realm,sample.units.id)
  # sample.units.id <- sample.units.id[,c("x","y","sample.units.id")]
  colnames(sample.units.id)[1:3] <- c("lon","lat","id.sample")
  
  ##Aggregate to sample regions 
  sample.unit.mean <- aggregate(cbind(Mean.Temp,Temp.seasonality,Annual.Precip,Precip.seasonality)~id.sample,data=sample.units.id,mean)
  colnames(sample.unit.mean) <- c("id.sample",paste0(vars, ".m"))
  sample.unit.var <- aggregate(cbind(Mean.Temp,Temp.seasonality,Annual.Precip,Precip.seasonality)~id.sample,data=sample.units.id,var)
  colnames(sample.unit.var) <- c("id.sample",paste0(vars, ".v"))
  sample.unit.var[is.na(sample.unit.var)] <- 0
  sample.unit.all <- merge(sample.unit.mean,sample.unit.var,by="id.sample",all.x=TRUE)
  
  #Create blocks, dividing polygons into orthogonal blocks on the basis of the climate data.
  blocks <- block(sample.unit.all, n.tr=10, id.vars='id.sample') ## decide the number of blocks to use i.e. n.tr=5 or 10
  blocks <- assignment(blocks, namesCol=as.character(1:10))$assg[[1]][1:10] ## assign subpols to one of 10 blocks 
  blocks <- Reduce(rbind, mapply(function(id.sample, block) data.frame(id.sample, block), id.sample=blocks, block=as.list(1:10),SIMPLIFY=F) )## turn this into a dataframe 
  blocks <- blocks[!is.na(blocks$id.sample),] ## remove any polygons that haven't been assigned to a block (we deal with this later).
  clidat <- merge(sample.unit.all, blocks, by=c('id.sample'), all=T) 
  clidat$block[is.na(clidat$block)] <- sample(1:10,1)
  clidat <- clidat[,c("id.sample","block")]
  
  t <- merge(sample.units.id, clidat,by="id.sample",all.x=TRUE)
  # head(t)
  # plot(t$lon, t$lat, col=t$block, pch=19, xlab="Longitude", ylab="Latitude", cex=0.1)
  # names(sample.units.id)
  write.csv(t, paste0(getwd(),"\\GitHub Version\\Data\\Explanatory layers\\Blocks\\",r,"_block.csv"))
  print(r)
  rm(list=c("sample.unit.realm","sample.units.id","sample.unit.mean", "sample.unit.var", "sample.unit.all", "blocks", "clidat", "t"))
  
}
