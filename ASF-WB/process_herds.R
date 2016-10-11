#### PROCESS_HERDS library (V. 0.15.1)
####
#### Includes the following functions:
####
#### updateHerds
#### DIRinf3
#### INDflex
#### LASinf
#### constructAInfHerd
####



################################################################
### updateHerds adjusted for the Danish project (TH and LEC) ###
################################################################
updateHerds <- function () {
  
  ## Depopulating diagnosed and tagged herds

  if (dim(depopQueue)[1]>0){
   ## Remove duplets and those already being depopulated
    tmpIndex<-!duplicated(depopQueue[,1]) #& !(depopQueue[,1]%in%beingDepoped[,1])
    depopQueue<<-depopQueue[tmpIndex,,drop=FALSE]

    getCHR <- aHerd$chr[depopQueue[1]]
    extraCull <- which(aHerd$chr%in%getCHR)
    extraCull <- extraCull[ !extraCull%in%depopQueue[1] ]
    if(length(extraCull)>0) {
      depopQueue   <<- rbind(depopQueue,cbind(extraCull,gTime))
      InfEC        <- extraCull[aHerd$status[extraCull]%in%c(2,3,4,7)]
      if(length(InfEC)>0){
        aHerd$Diagnosed[InfEC]     <<- TRUE
        aHerd$diagnosisTime[InfEC] <<- gTime
        aInfHerd$setDiagnosed(InfEC)
     }
    }
   ## Move herds from queue to being depoped
    cullMat <-matrix(numeric(0),ncol=3)
    cullMat <- cbind(depopQueue,aHerd$herdSizeCull[depopQueue[,1]])
    
      if(sum(cullMat[,3])<= Capacity){
        beingDepoped<<-rbind(beingDepoped,
                                   cbind(cullMat[,1],1))
             toCull <- 1:dim(cullMat)[2]
        }
  
     else{
      cumNumAnimals<-cumsum(cullMat[,3])
      cullCat<- which(cumNumAnimals<=Capacity)

      partial <- which(!(cumNumAnimals<=Capacity))[1] 
      if(!is.null(partial)) {
         cullMat[partial,3] <- (cumNumAnimals[partial]-Capacity)
         aHerd$herdSizeCull[ cullMat[partial,1] ] <<- cullMat[partial,3]
        }
      if(length(cullCat)>0) 
               beingDepoped<<-rbind(beingDepoped, cbind(cullMat[cullCat,1],1))

      toCull<- cullCat
                    
     }## EndOf else    


   if(length(toCull)>0) depopQueue<<-depopQueue[-toCull,,drop=FALSE]

  if (dim(beingDepoped)[1]>0){
      tmpID<-beingDepoped[,1]
      chronicle$addEntry(itn=iteration, state=5, newInfection=0, time=gTime,
                         changedHerdsList=aHerd$ID[tmpID],
                         HerdSize=aHerd$herdSize[tmpID],
                         HerdTypes=aHerd$herdType[tmpID])

      aHerd$timeToRemovedPZ[ tmpID ] <<- gTime 
      aHerd$timeToRemovedSZ[ tmpID ] <<- gTime 
      aHerd$status[ tmpID ]          <<- 5
      aHerd$timeCulled[tmpID]        <<- gTime
      aInfHerd$delInf(tmpID,ignoreNotInfected=TRUE)
      
      gDepopulationCount[iteration]<<-
        gDepopulationCount[iteration]+length(tmpID) 
      beingDepoped<<-matrix(numeric(0),ncol=2)

  if(Detailed){
      DepopMatOut<<- rbind(DepopMatOut,cbind(iteration,gTime,tmpID))
      if(dim(DepopMatOut)[1]>= DumpData){
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
         NAMED <- paste(runID,"DepopHerds.txt",sep="-")
         write.table(DepopMatOut,NAMED,append=TRUE,col.names = F,row.names = F)
         DepopMatOut<<- matrix(numeric(0),ncol=3)
         }
        }
    }
  }##EndOf: if (dim(depopQueue)[1]>0)

  ## Simulate one day of intra herd dynamics:
  aInfHerd$simDay()

}##EndOF function updateHerds2

###################################
### DIRinf3 (Adjusted for the Danish ASF modeling project by TH)
###
### Selecting herds infected by direct contact
###################################

DIRinf3<- function(Lambda,distprob,DCtoRiskDist,RiskInf,MovMat,restMovedSize,label,tStart=0,tEnd=Inf){

if(verbose){
    cat("DIRinf3 envir ");print(environment())
    print(parent.env(environment()))
    cat("Dist exists? ");print(exists("Dist"))
 }
pMat<-ConstpMat(Dist,probList$DistCat,distprob)

  
list(  
       init= function(){}
       ,
       day=function(){ 
         if ( gTime>tStart & gTime<tEnd){
           if (verbose) cat("Entered DIRinf3$day()\n")
              
           ## first determine the number of outgoing contacts per infected herd. then determine the number of moved animals per batch
           ## and then determine the probability that at least one infected animal is infected and then determine the number of
           ## infectious contacts and then determine the number of infected animals in the batches and make sure that there is
           ## at least one infected animal in each infectious batch

           numDC        <- rpois(length(aHerd[[Lambda]][infHerdNums]),aHerd[[Lambda]][infHerdNums])*aHerd$relDC[infHerdNums]
           MovedAnimals <- sapply(aHerd$NumMovAnimal[infHerdNums],function(x) eval(x,list(n=1)))
           MovedAnimals[ MovedAnimals>=aHerd$herdSize[infHerdNums] ] <- ceiling(aHerd$herdSize[infHerdNums[MovedAnimals>=aHerd$herdSize[infHerdNums]]]/restMovedSize) # number moverd animals from a herds should not
                                                                                                               # exceed the herd size. if so, then it is restricted to the median number of moved animals based on the
                                                                                                              # the data from the movement datanase; 1/10 and 1/35 of the herd for weaners and sows, respectively
           InfnessHerds <- (1-(1-aInfHerd$getInfnessDC(infHerdNums))^MovedAnimals) # determine the probability that the contact is infectious
           numDC        <- RandContacts(numDC * InfnessHerds * aHerd[[RiskInf]][infHerdNums])
           
          if (!identical(infHerdNums,infHerdNumsLast)) {
             if(verbose)  print("Calculating newInfected distance matrix")
             aHerd[[DCtoRiskDist]] <<- t(t(pMat$get(infHerdNums)))
           }
           
           AllNewInfs = NULL
           AllNewAnimals = NULL
           for(i in 1:length(infHerdNums)) {
             if (numDC[i]>0) {
               curHerd <- infHerdNums[i] 
               curHerdTyp <- aHerd$herdType[curHerd] ## bring the herd type category of the infected herd "i" 
                
              ### Newest change 23/06/2015 (TH)     
               Prob <- MovMat[curHerdTyp,aHerd$herdType]
               Prob <- Prob * aHerd[[DCtoRiskDist]][,i]
               contactHerds <- sample(aHerd$ID[aHerd$herdSize>MovedAnimals[i]],numDC[i],rep=T,prob=Prob[aHerd$herdSize>MovedAnimals[i]])
               contactHerds <- unique(contactHerds)
             
 
               if(length(contactHerds>0)){
                  traceMatrix <<- rbind(traceMatrix,
                       cbind(gTime,contactHerds,curHerd,label))
                  }# EndOf if
               redTraceMatrix <- which((traceMatrix[,1]+ TracePeriod)<gTime)
               if(length(redTraceMatrix > 0)) 
                  traceMatrix <<- traceMatrix[-redTraceMatrix,,drop=FALSE]

               dProbInfection <- rep(1-.001*Tstoch,length(contactHerds))
               dProbInfection <- dProbInfection * aHerd$relDC[contactHerds] * (aHerd$status[contactHerds]==1)
                                              
           
               ## Call the appropriate Tinfs to determine newInfected animals
               newinfanimals <- Tinfs(contactHerds,dProbInfection,aHerd$Sus) ## determine if the herd will be infected
               newInf <- contactHerds[as.logical(newinfanimals)] ## get only new infected herds

               if(length(newInf)>0) {
                 AllNewInfs = c(AllNewInfs,newInf)
                 AllNewAnimals = c(AllNewAnimals,newinfanimals[as.logical(newinfanimals)])
                 aHerd$status[newInf] <<-3        # infected status
                 aHerd$timeInfected[newInf] <<-gTime
                 ## Add infected individual(s) to each new herd; no latent period
                 aHerd$infSource[newInf] <<- curHerd  # record infection source herd
               aHerd$infMode[newInf] <<-label
               }
             }
           }

         ##LEC Adding new herds to aInfHerd
           if(length(AllNewInfs)>0){
             if(any(duplicated(AllNewInfs)))
               print("Duplicated direct infs")
             ##LEC All newly introduced animals are subclinical - should change
             ##LEC but requires reconsidering Tinfs and sampling from the
             ##LEC originating herd instead of scaling risk with proportion
             ##LEC of infected animals.
             aInfHerd$addInf(AllNewInfs,cbind(0,AllNewAnimals,0),1)
           
             chronicle$addEntry(itn=iteration, state=3, newInfection=label,
                                time=gTime,changedHerdsList=aHerd$ID[AllNewInfs],
                                SourceHerds=aHerd$ID[aHerd$infSource[AllNewInfs]],
                                HerdSize=aHerd$herdSize[AllNewInfs],HerdTypes=aHerd$herdType[AllNewInfs])
            }
         }##EndOf if (gTime ...)
       }##EndOf DIRinf3$day()
       ,
       
       cleaniteration=function(){
         if (verbose) cat("Entered DIRinf3$inititeration()\n")
         pMat$wipe()
       }
       ) 
}##EndOf DIRinf()


###################################
### Indirect contact construction function
###
### Flexible function to determine the indirect infections
###################################
 INDflex<-function(lambda,distprob,relCont,pMatIDC,RiskInf,probMatrix=NULL,Reduction=1,Abattoir=FALSE,label,checking=FALSE,tStart=0,tEnd=Inf){ ## RateIND represent the rate of contacts out of a farm

   pMat<-ConstpMat(Dist,probList[['DistCat']],distprob)
 
   list(
       init= function(){}
       ,
       day=function(){
          if ( gTime>tStart & gTime<tEnd){
           if (verbose) cat("Entered INDflex$RandContact()\n")

          RandCont<-rpois(length(aHerd[[lambda]][infHerdNums]),aHerd[[lambda]][infHerdNums])*aHerd[[relCont]][infHerdNums]

          if(Abattoir){
            if(!outbreakDetected)
             RandCont<- rexp(length(RandCont),rateBefD) * RandCont 
             
            else 
             RandCont <- rexp(length(RandCont),rateAftD) * RandCont 
              
              }
                          
          dInfTmp <- aInfHerd$getInfnessIDC(infHerdNums)

            if(!outbreakDetected)
              NRC<-RandContacts(RandCont* aHerd[[RiskInf]][infHerdNums] * dInfTmp) 
 
            else
              NRC<-RandContacts(RandCont* aHerd[[RiskInf]][infHerdNums] * dInfTmp * Reduction)  

 
            
          if (!identical(infHerdNums,infHerdNumsLast)) {
             if(verbose)  print("Calculating newInfected distance dependent risk matrix") 
             aHerd[[pMatIDC]]<<- t(t(pMat$get(infHerdNums)))
             }

           AllNewInfs = NULL
           AllNewAnimals = NULL
           for(i in 1:length(infHerdNums)) {
             if (NRC[i]>0) {
               curHerd <- infHerdNums[i]

               ### Newest change 25/5/2011 (TH)     
                if (Abattoir){ 
                   Prob <- probMatrix[aHerd$herdType[curHerd],aHerd$herdType]
                   Prob <- Prob * aHerd[[pMatIDC]][,i]
                   }

                else{
                    Prob <- aHerd[[pMatIDC]][,i]

                    }
                contactHerds <- sample(aHerd$ID,NRC[i],rep=T,prob=Prob)              

                if(length(contactHerds>0)){
                  traceMatrix <<- rbind(traceMatrix,
                       cbind(gTime,contactHerds,curHerd,label))
                  }# EndOf if

                redTraceMatrix <- which((traceMatrix[,1]+ TracePeriod)<gTime)
                if(length(redTraceMatrix > 0)) 
                  traceMatrix <<- traceMatrix[-redTraceMatrix,,drop=FALSE]

                dProbInfection <- rep(1-.001*Tstoch,length(contactHerds))
                            
                dProbInfection <- dProbInfection*(aHerd[[relCont]][contactHerds])*(aHerd$status[contactHerds]==1)
                                  
                dProbInfection <- dProbInfection*aHerd$RelSusceptibility[contactHerds]

                ## Call the appropriate Tinfs to determine newInfected animals
                newinfanimals <- Tinfs(contactHerds,dProbInfection,aHerd$Sus)
                newInf <- contactHerds[newinfanimals>0]
                aHerd$status[newInf] <<-2  # latent infection status
                            
                if(length(newInf)>0) {
                 AllNewInfs = c(AllNewInfs,newInf)
                 AllNewAnimals = c(AllNewAnimals,newinfanimals[as.logical(newinfanimals)])
                 aHerd$infSource[newInf] <<- curHerd  # record infection source herd 
                 aHerd$timeInfected[newInf] <<-gTime
                 aHerd$infMode[newInf] <<-label
               
                 }# End if
               }# End if
             }# End for

              ## Adding new herds to aInfHerd
               if(length(AllNewInfs)>0){
                 if(any(duplicated(AllNewInfs)))
                   print("Duplicated indirect infs")
                    aInfHerd$addInf(AllNewInfs,cbind(AllNewAnimals,0,0),0)

               chronicle$addEntry(itn=iteration, state=2, newInfection=label,
                                time=gTime,changedHerdsList=aHerd$ID[AllNewInfs],
                                SourceHerds=aHerd$ID[aHerd$infSource[AllNewInfs]],
                                HerdSize=aHerd$herdSize[AllNewInfs],HerdTypes=aHerd$herdType[AllNewInfs])
                    }               
                 }# End of if
               } # End of function day
            ,
       cleaniteration=function(){
         if (verbose) cat("Entered INDflex$inititeration()\n")
         pMat$wipe()
               }
  
       )# End of list
}## End of INDflex

        

###################################
### LASinf as changed by TH
###
### Determine indirect infections through local area spread
###################################
LASinf <- function(localsize,label,tStart=0,tEnd=Inf) {

  ## Functions: 'day' and 'inititeration'
  list(
       init= function(){}
       ,
       day=function(){
         if ( gTime>tStart & gTime<tEnd){
           AllNewInfs = NULL
           AllNewAnimals = NULL

         for(i in 1:length(infHerdNums)) { 
           spreadProb <- Dist$get(infHerdNums[i])

           ##categorize the distances according to the predefined categories in the ASFoptions 
           ## file and then get the corresponding probabilities
           spreadProb <- colSums(sapply(spreadProb,function(x) x>=LocSpLim))+1
           spreadProb <- DistList$distProb[spreadProb]

           ## get infectiousness of the infected herd 
           dInfTmp <- aInfHerd$getInfnessIDC(infHerdNums[i]) 
           
             LAScontact <- (1:gMaxHerds)[spreadProb>0 & aHerd$status==1 ]

             ## susceptible herds that will be infected is dependent on the probability distance from the infected herd, the infectiousness of the infected herd
     
             newinfanimals <- Tinfs(LAScontact,(spreadProb[LAScontact] * dInfTmp ),aHerd$Sus)
             newInf <- LAScontact[newinfanimals>0]  # list of newly infected herds
             newinfanimals<-newinfanimals[newinfanimals>0]
           
             ## Some recieving herds are not fully susceptible:
             tmp<-aHerd$RelSusceptibility[newInf] > runif(length(newInf))
             newInf<-newInf[tmp]
             newinfanimals<-newinfanimals[tmp]

             if(sum(newInf)>0) {
               AllNewInfs = c(AllNewInfs,newInf)
               AllNewAnimals = c(AllNewAnimals,newinfanimals)
               aHerd$status[newInf]       <<- 2                 # latent infection status
               aHerd$infSource[newInf]    <<- infHerdNums[i]    # record infection source herd 
               aHerd$timeInfected[newInf] <<-gTime
               aHerd$infMode[newInf]      <<-label
             }
           }##EndOf for
           ## Adding new herds to aInfHerd
           if(length(AllNewInfs)>0){
             if(any( dupIndex<-duplicated(AllNewInfs) )){
               print("Duplicated LAS infs")
               AllNewInfs<-AllNewInfs[!dupIndex]
               AllNewAnimals<-AllNewAnimals[!dupIndex]
             }
             ##LEC All newly infected animals are latent:
             aInfHerd$addInf(AllNewInfs,cbind(AllNewAnimals,0,0),0)
             
             chronicle$addEntry(itn=iteration, state=2, newInfection=7, time=gTime,
                                changedHerdsList=aHerd$ID[AllNewInfs],
                                SourceHerds=aHerd$ID[aHerd$infSource[AllNewInfs]],
                                HerdSize=aHerd$herdSize[AllNewInfs],HerdTypes=aHerd$herdType[AllNewInfs])
           }
         }##EndOf if (gTime ...)
       }##EndOf LASinf$day()
       ,
       cleaniteration=function(){
         if (verbose) cat("Entered LASinf$inititeration()\n")
       }
       )
}##EndOf LASinf()

#######################################################
## Infection from contact to wild boar in a          ##
## situation where ASF is endemic in wild boar. (TH) ##
#######################################################
 WildBoar <- function(relCont,RiskVar,Reduction=1,ProbCont,PrevelanceWB=0.031,label,tStart=0,tEnd=Inf){

  list(
       init= function(){}
       ,
       day=function(){
         if ( gTime>tStart & gTime<tEnd){

     ## Here we determine the frequency of contact (lambda) to wild boar that is dependent on risk profile in relation to WB, which is
     ## dependent on the distance to the nearest wild boar population or actually habitat.
          FerqCont <- LambdaWB[aHerd[[RiskVar]]]

     ## Here we determine the number of time WB may have potential contact to the herds.  
           if(!outbreakDetected) RandCont <- rpois(length(FerqCont),FerqCont)
     ## If outbreak is detected, then there is a possibility that contacts to WB by all herds can be reduced by a "Reduction" factor.
           else                  RandCont <- rpois(length(FerqCont),FerqCont) * aHerd[[relCont]] * Reduction
     ## Determine randomly whether their will be today at least one infectious contact to any of the herds
            Prob <- RandCont * aHerd[[ProbCont]] * PrevelanceWB * as.numeric(aHerd$status==1) * aHerd$RelSusceptibility 
            NumCont<-runif(length(RandCont))<=Prob # the probability that the contact is by an infected WB

     ## Keep only the herds that are infected by WB
            if(sum(NumCont)>0) {
             NewInfs <- which(NumCont)
             if(length(NewInfs)>0){
## Update information about newly infected herds
                aHerd$status[NewInfs]       <<- 2  # latent infection status
                aHerd$infSource[NewInfs]    <<- -1 # record infection source herd (WB) 
                aHerd$timeInfected[NewInfs] <<- gTime
                aHerd$infMode[NewInfs]      <<- label   
                aInfHerd$addInf(NewInfs,cbind(rep(1,length(NewInfs)),0,0),0)

## Update the chronicale of new information
               chronicle$addEntry(itn=iteration, state=2, newInfection=label,
                                time=gTime,changedHerdsList=aHerd$ID[NewInfs],
                                SourceHerds=0,
                                HerdSize=aHerd$herdSize[NewInfs],HerdTypes=aHerd$herdType[NewInfs])

           }#End of if(length( 
          }#End of if(sum(
        }##EndOf if (gTime>tStart ...)
       }##EndOf day()
       ,
       cleaniteration=function(){
         if (verbose) cat("Entered wildBoar$inititeration()\n")
       }
       )
}##EndOf WildBoar ()

#######################################################
## Infection from contact to wild boar in a          ##
## situation where ASF is endemic in wild boar. (TH) ##
#######################################################
 LSWildBoar <- function(RiskVar,PrevelanceWB=0.031,label,tStart=0,tEnd=Inf){

  list(
       init= function(){}
       ,
       day=function(){
         if ( gTime>tStart & gTime<tEnd){


     ## Here we determine the probability of infectious contact through local spread to wild boar.
          ProbCont <- LocalSPWB[aHerd[[RiskVar]]]
     ## Determine randomly whether their will be new infected herds
            ProbCont <- ProbCont * as.numeric(aHerd$status==1) * aHerd$RelSusceptibility * PrevelanceWB
            Infect <- runif(length(ProbCont))<=ProbCont# the probability that the contact is by an infected WB

     ## Keep only the herds that are infected by WB
            if(sum(Infect)>0){ 
             NewInfs <- which(Infect)
             if(length(NewInfs)>0){
## Update information about newly infected herds
                aHerd$status[NewInfs]       <<- 2  # latent infection status
                aHerd$infSource[NewInfs]    <<- -1 # record infection source herd (WB) 
                aHerd$timeInfected[NewInfs] <<- gTime
                aHerd$infMode[NewInfs]      <<- label   
                aInfHerd$addInf(NewInfs,cbind(rep(1,length(NewInfs)),0,0),0)

## Update the chronicale of new information
               chronicle$addEntry(itn=iteration, state=2, newInfection=label,
                                time=gTime,changedHerdsList=aHerd$ID[NewInfs],
                                SourceHerds=-1,
                                HerdSize=aHerd$herdSize[NewInfs],HerdTypes=aHerd$herdType[NewInfs])

           }#End of if(length( 
          }#End of if(sum(
        }##EndOf if (gTime>tStart ...)
       }##EndOf day()
       ,
       cleaniteration=function(){
         if (verbose) cat("Entered wildBoar$inititeration()\n")
       }
       )
}##EndOf LSWildBoar()


######################
### constructAInfHerd
###
### Call once to construct aInfHerd object.
###
######################
constructAInfHerd<-function(){
  ## Internal variables:

  ## VarDef herds    Matrix with info on infected herds
  ## VarDes Column names are in herdsCol. A matrix without names is faster
  nColHerds<-19
  herds      <-matrix(numeric(0),ncol=nColHerds)
  DeadMat    <-matrix(numeric(0),ncol=DaysDead)
  DeadMatSur <-matrix(numeric(0),ncol=DaysSurDead)

  delHerds<-matrix(numeric(0),ncol=nColHerds) ##Keeping deleted herds
  herdsCol<-c("Sus1","latent2","SubC3","Clinic4","Immune5","Total6","status7",
              "ID8","p9","NA.latDur10","NA.SubDur11","TClic12","TDiag13",
              "Diagnosed14","InfByDC15","TInfected16","herdType17",
              "Vaccinated18","TLastAnCli19")
  dimnames(herds)[[2]]<-herdsCol
  
  ## VarDef latent,subclinical  Used to model transitions between stages
  ## VarDes Matrix with a row for each farm and columns for
  ## VarDes   each day latent or subclinical or clinical
  latent <- matrix(numeric(0),ncol=ncol(herdtypes$latDurFreq))
  subclinical <- matrix(numeric(0),ncol=ncol(herdtypes$SubDurFreq))
  clinical <- matrix(numeric(0),ncol=ncol(herdtypes$CliDurFreq))

  "rpoly2"<-function(np){
    tabulate(sample(1:(length(np)-1),size=np[1],replace=TRUE,
                    prob=np[-1]),nbins=(length(np)-1))

  }
  
  list( 
# List holding all information to be access from outside
#################################################################
       getIDs=function(){ return(herds[,8])  }
       ,       #####################
       getstatus=function(){ return(herds[,7]) }
       ,       #####################
       getTDiag=function(){ return(herds[,13]) }
       ,       #####################
       getDelHerds=function(){ return(delHerds) }
       ,       #####################
       getHerds=function(){ return(herds) }
       ,       #####################
       setHerds=function(newHerds){ herds<<-newHerds }
       ,       #####################
       getDiagnosed=function(){ return(herds[,14]) }
       ,       #####################
       getDiagnosedIDs=function(){ return(herds[herds[,14]==1,8]) }
       ,       #####################
       getTClic=function(IDs) { ### Time the herd showed clinical signs (TH)
         Herds<-herds
         index<-match(IDs,Herds[,8])
         if(any(is.na(index))){
           warning("ConstructAInfHerd: Some IDs not in aInfHerd!")
           index<-index[!is.na(index)]
            }
           tmp <- Herds[index,12] 
           return(tmp) 
        } 
       ,       #####################
       setDiagnosed=function(IDs){
         index<-match(IDs,herds[,8])
         if (any(is.na(index))){
           #warning("ConstructAInfHerd: Some IDs not in aInfHerd!") this warning is stopped, because some herds will be culled and hence treated as
           # diagnosed if they have the same chr as a detected herd.
           index<-index[!is.na(index)]
         }
         herds[index,14]<<-1
       }
       ,       #####################
       getDelIDs=function(){ return(delHerds[,8]) }
       ,       #####################
       getInfnessDC=function(IDs=herds[,8]){
         index<-match(IDs,herds[,8])
         if (any(is.na(index))){
           warning("ConstructAInfHerd: Some IDs not in aInfHerd!")
           index<-index[!is.na(index)]
         }
         tmp<-(herds[index,2]+herds[index,3]+herds[index,4])/(herds[index,6]-(herds[index,5]*PerDeadAnim))

         return(tmp)
       }
       ,       #####################
       getInfnessIDC=function(IDs=herds[,8]){
         index<-match(IDs,herds[,8])
         if (any(is.na(index))){
           warning("ConstructAInfHerd: Some IDs not in aInfHerd!")
           index<-index[!is.na(index)]
         }

         tmp<-((herds[index,3]*InfPropSub)+(herds[index,4]) + ((rowSums(t(exp((-dim(DeadMat[index,,drop=FALSE])[2]+1):0*ImpDeadTime)*t(DeadMat[index,,drop=FALSE]))))*PerDeadAnim*DeadImpact))/
               (herds[index,6]-(herds[index,5]*PerDeadAnim))

         return(tmp)
       }

       ,       #####################
       anyInf=function(){ return(nrow(herds)>0) }
       ,       #####################
       ##numInfAnimals: three column matrix for latent, subC, and Clin
       addInf=function(IDs,numInfAnimals,DC){
         tmpID<-(IDs%in%herds[,8]) ##Only those that were not already infected
         if (any(tmpID)){
           IDs<-IDs[!tmpID]
           numInfAnimals<-numInfAnimals[!tmpID,,drop=FALSE]
           DC<-DC[!tmpID]

         }
         if(length(IDs)>0){
           tmp<-cbind(aHerd$Sus[IDs],0,0,0,0,aHerd$Sus[IDs],1,IDs,
                      aHerd$p[IDs],0,0,
                      Inf,Inf,0,DC,gTime,aHerd$herdType[IDs],
                      0,0)

           tmp[,1:4]<-tmp[,1:4]+cbind(-rowSums(numInfAnimals),numInfAnimals)
           tmp[,7]<-apply(tmp[,1:4,drop=FALSE],1,function(x) max((1:4)[x>=1]))##status
           ## Updating diagnosis if initial status is 4 (Clinical)
           if(any(tmpIndex<-(tmp[,7]==4))){
             tmp[tmpIndex,12]<-gTime  ## Time herd showed clinical signs
           #  tmp[tmpIndex,13]<-gTime+tmp[tmpIndex,19]
           }         
           herds<<-rbind(herds,tmp)
           ## Distributing the animals on No. days being latent or subclincal or clinical
           latent<<-rbind(latent,t(apply(cbind(tmp[,2],herdtypes$latDurFreq[tmp[,17],,drop=FALSE]),1,rpoly2)))
           subclinical<<-rbind(subclinical,t(apply(cbind(tmp[,3],herdtypes$SubDurFreq[tmp[,17],,drop=FALSE]),1,rpoly2)))
           clinical<<-rbind(clinical,t(apply(cbind(tmp[,4],herdtypes$CliDurFreq[tmp[,17],,drop=FALSE]),1,rpoly2)))
           DeadMat<<-rbind(DeadMat,t(sapply(IDs,function(x)rep(0,(DaysDead)))))
           DeadMatSur<<-rbind(DeadMatSur,t(sapply(IDs,function(x)rep(0,(DaysSurDead)))))


         }##EndOf if(length(tmpID))

       }

       ,       #####################
       delInf=function(IDs,ignoreNotInfected=FALSE){
         index<-match(IDs,herds[,8])
         if (any(is.na(index))){
           if (!ignoreNotInfected)
             warning("ConstructAInfHerd: Cannot remove IDs not in aInfHerd!",
                     IDs,index)
           index<-index[!is.na(index)]
         }
         if(length(index)>0){
           delHerds<<-rbind(delHerds,herds[index,])
           herds<<-herds[-index,,drop=FALSE]
           latent<<-latent[-index,,drop=FALSE]
           subclinical<<-subclinical[-index,,drop=FALSE]
           clinical<<-clinical[-index,,drop=FALSE]
           DeadMat<<-DeadMat[-index,,drop=FALSE]
           DeadMatSur<<-DeadMatSur[-index,,drop=FALSE]
         }
       }
       ,       #####################
       simDay=function(){
         if (nrow(herds)>0){
           ## Sus -> latent
           rS2L<-rbinom(nrow(herds),herds[,1],
                 1-exp(- herds[,9]*( (herds[,3]*InfPropSub)+(herds[,4]) + 
                 ((rowSums(t(exp((-dim(DeadMat)[2]+1):0*ImpDeadTime)*t(DeadMat))))*PerDeadAnim*DeadImpact) ) / (herds[,6]-(herds[,5]*PerDeadAnim))))


           ## Updating totals
           herds[,1:5]<<-herds[,1:5]+cbind(-rS2L,rS2L-latent[,1],
                                           latent[,1]-subclinical[,1],
                                           subclinical[,1]-clinical[,1],
                                            clinical[,1])
           ## Distributing new latent infected
           newLat<-apply(cbind(rS2L,herdtypes$latDurFreq[herds[,17],,drop=FALSE]),1,rpoly2)
           ## Distributing new subclinical infected
           newSubC<-apply(cbind(latent[,1],herdtypes$SubDurFreq[herds[,17],,drop=FALSE]),1,rpoly2)
           ## Distributing new clinical infected
           newCli<-apply(cbind(subclinical[,1],herdtypes$CliDurFreq[herds[,17],,drop=FALSE]),1,rpoly2)
           ## Update intra latent and subclinical and clinical waiting states
           latent      <<- cbind(latent[,-1,drop=FALSE],0)+t(newLat)
           subclinical <<- cbind(subclinical[,-1,drop=FALSE],0)+t(newSubC)
           DeadMat     <<- cbind(DeadMat,clinical[,1])
           DeadMatSur  <<- cbind(DeadMatSur,clinical[,1])
           DeadMat     <<- DeadMat[,-1,drop=FALSE]
           DeadMatSur  <<- DeadMatSur[,-1,drop=FALSE]
           clinical    <<- cbind(clinical[,-1,drop=FALSE],0)+t(newCli)
           ## Checking for herds that changed status
           ## status 2 -> 3 : ## status latent to subclinical
           tmpIndex<-herds[,7]==2 & herds[,3]>=1
           if(sum(tmpIndex)>0){
             chronicle$addEntry(itn=iteration, state=3, newInfection=0,
                                time=gTime, changedHerdsList=herds[tmpIndex,8],
                                HerdSize=aHerd$herdSize[herds[tmpIndex,8]],
                                HerdTypes=aHerd$herdType[herds[tmpIndex,8]])
             herds[tmpIndex,7] <<-3      ## status: SubClinical
             aHerd$status[herds[tmpIndex,8]]<<-3
           }
           ## status 3 -> 4 : # status subclinical to clinical
           tmpIndex<-herds[,7]==3 & herds[,4]>=1
           if(sum(tmpIndex)>0){
             chronicle$addEntry(itn=iteration, state=4, newInfection=0,
                                time=gTime, changedHerdsList=herds[tmpIndex,8],
                                HerdSize=aHerd$herdSize[herds[tmpIndex,8]],
                                HerdTypes=aHerd$herdType[herds[tmpIndex,8]])

             herds[tmpIndex,7] <<-4      ## status: Clinical
             herds[tmpIndex,12]<<-gTime  ## Time herd showed clinical signs
             #herds[tmpIndex,13]<<-gTime+herds[tmpIndex,19]
             aHerd$status[herds[tmpIndex,8]]<<-4
           }
           ## status 4 -> 7 : status clinical to immune
           tmpIndex<-(herds[,7]==4 & herds[,5]==herds[,6]) | (herds[,7]==4 & herds[,2]==0&herds[,3]==0&herds[,4]==0&herds[,5]>0&rowSums(DeadMat)==0)
           if(sum(tmpIndex)>0){
             chronicle$addEntry(itn=iteration, state=7, newInfection=0,
                                time=gTime, changedHerdsList=herds[tmpIndex,8],
                                HerdSize=aHerd$herdSize[herds[tmpIndex,8]],
                                HerdTypes=aHerd$herdType[herds[tmpIndex,8]])

             herds[tmpIndex,7] <<- 7      ## status: Immune
             #herds[tmpIndex,12]<<-gTime  ## Time herd became Immune
             aHerd$status[herds[tmpIndex,8]] <<- 7
             aHerd$immuneTime[herds[tmpIndex,8]] <<- gTime
           }
           ## status 7 -> 1 : status immune to susceptible
           tmpIndex <- herds[,7]==7 & herds[,5]<herds[,6]
           if(sum(tmpIndex)>0){
             chronicle$addEntry(itn=iteration, state=1, newInfection=0,
                                time=gTime, changedHerdsList=herds[tmpIndex,8],
                                HerdSize=aHerd$herdSize[herds[tmpIndex,8]],
                                HerdTypes=aHerd$herdType[herds[tmpIndex,8]])

             herds[tmpIndex,7] <<- 1      ## status: Susceptible again
             aHerd$status[herds[tmpIndex,8]] <<- 1
             aHerd$SusAgain[herds[tmpIndex,8]] <<- gTime
             aInfHerd$delInf(herds[tmpIndex,8])
             ## Update the indexes of infectious herds
             tmpTagged<-which( aHerd$status==6 & aHerd$timeInfected<Inf & aHerd$SusAgain==0)
             infHerdNums <<- unique(c(tmpTagged,
                             (aInfHerd$getIDs())[aInfHerd$getstatus()%in%(3:4)]))


           }

         ## update the number sick and dead animals
         aHerd$Mortality[herds[,8]] <<- (herds[,4]+(herds[,5]*PerDeadAnim)) 
         aHerd$Survived[herds[,8]]  <<- (herds[,5]*(1-PerDeadAnim))

         }##EndOf if (nrow(herds)>0)
       }
       ,
        #############################
##### Get the dead animals animals during last week and sick animals today for the infected herds
       getDead=function(IDs){
         index<-match(IDs,herds[,8])
         if (any(is.na(index))){
           warning("ConstructAInfHerd: Some IDs in zones are not in aInfHerd")
           index<-index[!is.na(index)]
         }
         tmp <- round(rowSums(DeadMatSur)[index]*PerDeadAnim)
          return(tmp)
        }
,
               #####################
       getInfected=function(IDs){
         index<-match(IDs,herds[,8])
         if (any(is.na(index))){
           warning("ConstructAInfHerd: Some IDs in not in aInfHerd!")
           index<-index[!is.na(index)]
         }
         tmp<-round(herds[index,3]+herds[index,4]+(herds[index,5]*(1-PerDeadAnim)))##Include infected animals that survived death

         return(tmp)
        }
       ,       #####################
       wipe=function(){
         delHerds<<-matrix(numeric(0),ncol=nColHerds)
         herds<<-matrix(numeric(0),ncol=nColHerds)
         dimnames(herds)[[2]]<<-herdsCol
         latent <<-matrix(numeric(0),ncol=ncol(herdtypes$latDurFreq))
         subclinical <<- matrix(numeric(0),ncol=ncol(herdtypes$SubDurFreq))
         clinical <<- matrix(numeric(0),ncol=ncol(herdtypes$CliDurFreq))
         DeadMat<<-matrix(numeric(0),ncol=DaysDead)
         DeadMatSur<<-matrix(numeric(0),ncol=DaysSurDead)
       }
       )##EndOf list
}##EndOf constructAInfHerd
