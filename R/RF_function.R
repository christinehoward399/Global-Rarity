## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## Random Forest Function
## Adapted by C Howard from code provided by R. Bagchi and D. Baker
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

RF_funct<-function(source.data,outDir,realm,taxa,blocks, Response){
    
    #change to read in files correctly
    data.model <- source.data
    data.model<-data.model[complete.cases(data.model),]
    
    predNames <- sapply(c(CoVars),function(x,...) colnames(data.model)[which(colnames(data.model)==x)])#get predictor names
    formula <- as.formula(paste(Response," ~ ",paste(predNames, collapse=" + "),sep=''))  
    
    ##Set up parallel
    sfInit(parallel=TRUE, cpus=CP)
    sfLibrary(randomForest)#;sfLibrary(PresenceAbsence); sfLibrary(modEvA)
    sfExport(list=c("data.model","formula","realm","taxa","Response"))
    
    mod <- lapply(1:3,function(m){
      mod <- sfLapply(blocks,function(blockNo){
        
        cat('\n',"block = ",blockNo,"; mtry = ",m,'\n',sep=" ")
        fit.blocks <- subset(data.model, get(realm) != blockNo)
        test.block <- subset(data.model, get(realm) == blockNo)
        init.ntree <- 1000; init.Rsq <- 0.01; mod.improve <- "TRUE"
        
        while(mod.improve == "TRUE"){   
          
          model1 <- randomForest(formula, data=fit.blocks, ntree = init.ntree, mtry=m)   
          PRED <- predict(model1,newdata=test.block,type="response",se.fit=FALSE)
          PRED <- as.data.frame(PRED)
          test.block$id<-row.names(test.block)
          eval.data.rsq <- cbind(test.block[,c("id",Response)],PRED)
          
          ## Evaluate model fit
          SSE<-sum((eval.data.rsq[,Response]-eval.data.rsq$PRED)^2, na.rm=T)
          TSS<-sum((eval.data.rsq[,Response]-mean(fit.blocks[,Response]))^2, na.rm = T)
          Rsq<-1-(SSE/TSS)
          
          try(if(((Rsq/init.Rsq)) > 1.01){mod.improve <- "TRUE";init.Rsq <- Rsq;init.ntree <- init.ntree +500}else{mod.improve <- "FALSE"; best.ntree<-(init.ntree-500)}, silent=T)
          if(is.na(Rsq)){mod.improve <- "FALSE"; best.ntree<-(init.ntree-500); Rsq<-init.Rsq}
        }  
        return(list(block = blockNo, mtry=m, ntree=best.ntree, Rsq = init.Rsq))
        cat("\n")
      })
      return(mod)
    })
    
    ##Find optimum values
    evaluation.tab <- as.data.frame(matrix(unlist(mod),ncol=4,byrow=TRUE)); colnames(evaluation.tab) <- c("block","mtry","ntrees","Rsq")
    Rsq.mean.mtry <- aggregate(Rsq~mtry,FUN=mean,data=evaluation.tab)
    opt.mtry <-  Rsq.mean.mtry[which(Rsq.mean.mtry[,2]==max(Rsq.mean.mtry[,2])),1]
    opt.ntree <- max(evaluation.tab$ntrees)
    
    ##Set up parallel
    sfExport(list=c("data.model","formula","realm","opt.mtry","opt.ntree","Response"))
    
    #fit final models to data
    mod <- sfLapply(blocks,function(blockNo,m,nt){
      
      cat('\n', "Fitting final model to block",blockNo,'\n',sep=" ")
      fit.blocks <- subset(data.model, get(realm) != blockNo)
      test.block <- subset(data.model, get(realm) == blockNo)
      
      model1 <- randomForest(formula, data=fit.blocks, ntree = nt, mtry=m, importance = T)
      PRED <- predict(model1,newdata=test.block,type="response",se.fit=FALSE)
      PRED <- as.data.frame(PRED)
      test.block$id<-row.names(test.block)
      eval.data.Rsq <- cbind(test.block[,c("id",Response)],PRED)
      eval.data.Rsq<-eval.data.Rsq[complete.cases(eval.data.Rsq),]
      
      SSE<-sum((eval.data.Rsq[,Response]-eval.data.Rsq$PRED)^2)
      TSS<-sum((eval.data.Rsq[,Response]-mean(fit.blocks[,Response]))^2)
      Rsq<-1-(SSE/TSS)
      
      return(list(realm=realm, block=blockNo, Rsq=Rsq,mod=model1, mtry=m, ntree=nt))
    },m = opt.mtry, nt = opt.ntree)
    
    ##Save block models
    save(mod, file=paste(outDir,"/",taxa,"/", realm, "/",taxa,".",realm,".","model.output.RF.rda",sep=""),compress="bzip2")
    sfStop()
  }