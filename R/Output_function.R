## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## Function for setting up file structure
## Adapted by C Howard from code provided by R. Bagchi and D. Baker
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

Output_files <- function(outDir,realm, taxa){
  
  ##create directory to store all ENSEMBLE data
  outData <- outDir
  if(file.exists(outData) == FALSE) {dir.create(outData)}
  ##create a directory to store taxa specific data 
  if(file.exists(paste(outData, taxa, sep="\\")) == FALSE){ #create directory for storing split data
    dir.create(paste(outData,taxa, sep="\\"))}
  ## Create file to create zone specific data
  if(file.exists(paste(outData, taxa, realm,sep="\\")) == FALSE){ #create directory for storing split data
    dir.create(paste(outData,taxa, realm,sep="\\"))}
}