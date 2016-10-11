#### CONTROL MEASURES library Version 0.15.1
#### Implementation of control measures including
#### circular  zones, tracing, control of diagnosed herds, 
#### national standstill and surveillance of herds (TH)
#############################################################

SurvZone<-function(size=10,size2=3,effectDC,effectIMC,effectILC,label){
  if (!is.expression(effectDC))
    effectDC<-parse(text=effectDC)
  if (!is.expression(effectIMC))
    effectIMC<-parse(text=effectIMC)
  if (!is.expression(effectILC))
    effectILC<-parse(text=effectILC)

  effDCSZ<-numeric(0)
  effIMCSZ<-numeric(0)
  effILCSZ<-numeric(0)
  SZcenters<-NULL
  
  list(
       init=function(){ ## Initialize next iteration
         effDCSZ<<-1-eval(effectDC,list(n=gMaxHerds))
         effIMCSZ<<-1-eval(effectIMC,list(n=gMaxHerds))
         effILCSZ<<-1-eval(effectILC,list(n=gMaxHerds))
         SZcenters<<-rep(FALSE,gMaxHerds)
       }
       ,
       day=function(){ ## Daily update of variables
         ## First: Check if inCirc should be updated

         if( SZoneDuration < maxTime ){ ## If centers may be removed.
           tmpCenters <- ( aHerd$Diagnosed &
                          !(gTime>(aHerd$timeToRemovedSZ + SZoneDuration)))
         }
         else
           tmpCenters <- aHerd$Diagnosed
         
         
         if(sum(SZcenters!=tmpCenters)>0){ ## Update inCirc if changes
             SZcenters<<-tmpCenters
            }##EndOf if (sum ...)

        aHerd$inSurZone    <<- rep(FALSE,gMaxHerds)
        aHerd$InPZToo      <<- rep(FALSE,gMaxHerds)

        Centers <- which(SZcenters)
        for(i in 1:length(Centers)){
            CurCenter <- Centers[i]
            Distance  <- calDist(CurCenter)
            indexZone <-  Distance < size 
            indexPZone <- Distance <= size2       
            aHerd$inSurZone[indexZone] <<- TRUE
            aHerd$InPZToo[indexPZone]  <<- TRUE
         }

           ## only herd on surveillance zones
           inCircle <- aHerd$inSurZone & !aHerd$InPZToo

           ## Second: Update aHerd$relDC, aHerd$relIMC and aHerd$relILC
            aHerd$relDC[inCircle] <<-pmin(aHerd$relDC[inCircle],effDCSZ[inCircle])
            aHerd$relIMC[inCircle]<<-pmin(aHerd$relIMC[inCircle],effIMCSZ[inCircle])
            aHerd$relILC[inCircle]<<-pmin(aHerd$relILC[inCircle],effILCSZ[inCircle])    

  ### here only herds within the surveillance zones are processed.
  ### herds will get a visit, if they become included in a zone while they were never visited,
  ### or when they have been visited earlier, but after inclusion in a surveillance zone following a period they were outside a surveillance zone.
  ### herds in overlapping zones will get assigned a new visit (SV2) every SecSurVisitOLSZ number of days (by default 0) as long as they are still in a SZ
  ### and have had not been visited within the last DelayVisit number of days (by default 7 days) this means a new visit every DelayVisit days  

             surv1 <- which((inCircle & !aHerd$Diagnosed & aHerd$timeToSV1 < gTime & (aHerd$timeToSV2 < gTime) & aHerd$visitCount==0)|
                            (inCircle & !aHerd$Diagnosed & aHerd$timeToSV1 < gTime & (aHerd$timeToSV2 < gTime) & aHerd$visitCount>0  & aHerd$NotInZoneYest))

             surv2 <- which((inCircle & !aHerd$Diagnosed & ((aHerd$timeToSV2 + DelayVisit) < gTime) & aHerd$visitCount>0 & !aHerd$NotInZoneYest)) 

          if(firstSurvVisit) aHerd$timeToSV1[surv1] <<- gTime + DelayStartVisit

          aHerd$timeToSV2[surv1] <<- gTime + SecSurVisit     + DelayStartVisit
          aHerd$timeToSV2[surv2] <<- gTime + SecSurVisitOLSZ + DelayStartVisit

          ## inform whether herds are today in surveillance zones or not. this info will be used tomorrow
          aHerd$NotInZoneYest[!inCircle] <<- TRUE
          aHerd$NotInZoneYest[inCircle]  <<- FALSE


  ### herds that had a scheduled visit because they were in overlapping SZ will not get the visit unless they are still in the SZ
          IndexVisStop <- (aHerd$timeToSV2>=gTime)&!inCircle
            if(sum(IndexVisStop)>0){
              if(firstSurvVisit) aHerd$timeToSV1[IndexVisStop] <<- 1
                                 aHerd$timeToSV2[IndexVisStop] <<- 1
              }

          
  #if(Detailed){
   #   if(sum(inCircle)>0) SurvZoneMatOut<<- rbind(SurvZoneMatOut,cbind(iteration,gTime,which(inCircle)))
    #  if(dim(SurvZoneMatOut)[1]>= DumpData){
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
     #    NAMES <- paste(runID,"SurvHerds.txt",sep="-")
      #   write.table(SurvZoneMatOut,NAMES,append=TRUE,col.names = F,row.names = F)
       #  SurvZoneMatOut<<- matrix(numeric(0),ncol=3)
        # }
        #}
       
          aHerd$daysInZones[aHerd$inSurZone] <<- aHerd$daysInZones[aHerd$inSurZone] + 1
          aHerd$daysInSZones[aHerd$inSurZone] <<- aHerd$daysInSZones[aHerd$inSurZone] + 1
          aHerd$timeSatInSZone[(aHerd$inSurZone)&aHerd$timeSatInSZone==0] <<- gTime
 
        }#EndOf Day
       ,
       getIn=function(){  return(inCirc) }
       ,
       getLabel=function(){ return(label) }
        )##EndOf list
}##EndOf controlCirc5

########################################################
##
## Simulate the protection zone
##
########################################################

ProtZone<-function(size=3,effectDC,effectIMC,effectILC,label){
  if (!is.expression(effectDC))
    effectDC<-parse(text=effectDC)
  if (!is.expression(effectIMC))
    effectIMC<-parse(text=effectIMC)
  if (!is.expression(effectILC))
    effectILC<-parse(text=effectILC)

  effDCPZ<-numeric(0)
  effIMCPZ<-numeric(0)
  effILCPZ<-numeric(0)
  PZcenters<-NULL
  
  list(
       init=function(){ ## Initialize next iteration
         effDCPZ<<-1-eval(effectDC,list(n=gMaxHerds))
         effIMCPZ<<-1-eval(effectIMC,list(n=gMaxHerds))
         effILCPZ<<-1-eval(effectILC,list(n=gMaxHerds))
         PZcenters<<-rep(FALSE,gMaxHerds)
       }
       ,
       day=function(){ ## Daily update of variables
         ## First: Check if inCircPZ should be updated

         if( PZoneDuration < maxTime ){ ## If centers may be removed.
           tmpCenters <- ( aHerd$Diagnosed &
                          !(gTime>(aHerd$timeToRemovedPZ + PZoneDuration)))
         }
         else
           tmpCenters <- aHerd$Diagnosed
         
         
         if(sum(PZcenters!=tmpCenters)>0){ 
             PZcenters<<-tmpCenters
            }##EndOf if (sum ...)

        aHerd$inProtZone     <<- rep(FALSE,gMaxHerds)
        Centers <- which(PZcenters)
        for(i in 1:length(Centers)){
            CurCenter <- Centers[i]
            Distance  <- calDist(CurCenter)
               indexPZone <- Distance <= size 
               aHerd$inProtZone[indexPZone]  <<- TRUE
           }


           ## Second: Update aHerd$relDC, aHerd$relIMC and aHerd$relILC
            aHerd$relDC[aHerd$inProtZone] <<-pmin(aHerd$relDC[aHerd$inProtZone],effDCPZ[aHerd$inProtZone])
            aHerd$relIMC[aHerd$inProtZone]<<-pmin(aHerd$relIMC[aHerd$inProtZone],effIMCPZ[aHerd$inProtZone])
            aHerd$relILC[aHerd$inProtZone]<<-pmin(aHerd$relILC[aHerd$inProtZone],effILCPZ[aHerd$inProtZone])    

          ### here herds in the protection zones are processed.
            ## herds that had not been visited before or had been out of a protection zone and then came again in a protection zone are selected for 2 visits (PV1 and PV2)
            survProt1 <- which((aHerd$inProtZone & !(aHerd$Diagnosed) & (aHerd$timeToPV2 < gTime) & aHerd$visitCount==0)|
                               (aHerd$inProtZone & !(aHerd$Diagnosed) & (aHerd$timeToPV2 < gTime) & aHerd$visitCount> 0 & aHerd$NotInPZYest))


            ## herds that are still in a protection zone despite that DelayVisit (by default 7) days had passed since they had recieved the first and second PZ visits are selected and given a visit (PV2) every 
            ## DelayStartVisitOLPZ number of days (by default 0 days). this will make sure that the herds are visited within DelayVisit (by default 7) days 
            ## before the last protection zone is lefted. this is important because a zone may be lefted only after surveillance of the herds. 
            survProt2 <- which((aHerd$inProtZone & !(aHerd$Diagnosed) & (aHerd$timeToPV1 < (gTime+DelayVisit)) & (aHerd$timeToPV2 < (gTime+DelayVisit)) & 
                                aHerd$visitCount > 0 & !aHerd$NotInPZYest))

            aHerd$timeToPV1[survProt1] <<- gTime + DelayStartVisit
            aHerd$timeToPV2[survProt1] <<- gTime + SecProtVisit + DelayStartVisit
            aHerd$timeToPV2[survProt2] <<- gTime + DelayStartVisitOLPZ + DelayStartVisit

          IndexVisStopPV <- aHerd$timeToPV2>=gTime&!aHerd$inProtZone
            if(sum(IndexVisStopPV)>0){
              aHerd$timeToPV2[IndexVisStopPV] <<- 1
              aHerd$timeToPV1[IndexVisStopPV] <<- 1
            }

           aHerd$NotInPZYest[!aHerd$inProtZone] <<- TRUE
           aHerd$NotInPZYest[aHerd$inProtZone]  <<- FALSE
          ## herds that are in protection zone and had a planned SV visit will not get the SV visits as they will get 2 visits for being in PZ. 
          if(firstSurvVisit) aHerd$timeToSV1[aHerd$inProtZone] <<- 1
                             aHerd$timeToSV2[aHerd$inProtZone] <<- 1
        

  #if(Detailed){
  #    if(sum(aHerd$inProtZone)>0) ProtZoneMatOut<<- rbind(ProtZoneMatOut,cbind(iteration,gTime,which(aHerd$inProtZone)))
  #    if(dim(ProtZoneMatOut)[1]>= DumpData){
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
  #       NAMEP <- paste(runID,"ProtHerds.txt",sep="-")
  #       write.table(ProtZoneMatOut,NAMEP,append=TRUE,col.names = F,row.names = F)
  #       ProtZoneMatOut<<- matrix(numeric(0),ncol=3)
  #      }
  #     }
          aHerd$daysInZones[aHerd$inProtZone] <<- aHerd$daysInZones[aHerd$inProtZone] + 1
          aHerd$daysInPZones[aHerd$inProtZone] <<- aHerd$daysInPZones[aHerd$inProtZone] + 1
          aHerd$timeSatInPZone[(aHerd$inProtZone)&aHerd$timeSatInPZone==0] <<- gTime
 
        }#EndOf Day
       ,
       getIn=function(){ }
       ,
       getLabel=function(){ return(label) }
        )##EndOf list
}##EndOf ProtZone

##################################################################################
## Set herds for queue to surveillance and survay them based on the available   ##
## resources for Surveillance. Surveillance is based on number of herds/category##
## and is carried out based on herd level. (TH)                                 ##
##################################################################################

SurvZonesHerds <- function(){


list(
       init=function(){

       }
       ,
       day=function(){ 

### Herds that should be included in Queuing for surveillance and Select herds that will be serologically tested
# first select the herds for PV2 and SV2, where the herds will be tested.
 setQueue1 <- (aHerd$timeToPV2==gTime) & !(aHerd$Diagnosed) & !(aHerd$status%in%c(5,6))

## Find herds that have to get the PV1
   setQueue2 <-  aHerd$timeToPV1==gTime & !(aHerd$Diagnosed) & !(aHerd$status%in%c(5,6))

## If SV1 should be excuted then allow select the herds that should get this visit, otherwise do not excute the visit
 if(firstSurvVisit){
 setQueue3 <-  aHerd$timeToSV1==gTime & !(aHerd$Diagnosed) & !(aHerd$status%in%c(5,6))
   }
 else{
   setQueue3 <- FALSE
   }

#print(sum(setQueue3))

 ## Find herds that have to get the SV2
   setQueue4 <- aHerd$timeToSV2==gTime & !(aHerd$Diagnosed) & !(aHerd$status%in%c(5,6))

 ## Find herds that have to get the tracing visit of indirect contacts (IDC)
   setQueue5 <- aHerd$timeToVisitTraceIDC==gTime & !(aHerd$Diagnosed) & !(aHerd$status%in%c(5,6))

 ## Find herds that have to get the tracing visit of direct contacts (DC)
   setQueue6 <- aHerd$timeToVisitTraceDC==gTime & !(aHerd$Diagnosed) & !(aHerd$status%in%c(5,6))  


## Randomly select a proportion of the herds for PV1 and SV1 to be tested
 if(sum(setQueue2)>0) SerSetQueue2 <- rbinom(sum(setQueue2), size=1, prob=ProbSelPV1)
 if(sum(setQueue3)>0) SerSetQueue3 <- rbinom(sum(setQueue3), size=1, prob=ProbSelSV1)
 if(sum(setQueue4)>0) SerSetQueue4 <- rbinom(sum(setQueue4), size=1, prob=ProbSelSV2) 
 if(sum(setQueue5)>0) SerSetQueue5 <- rbinom(sum(setQueue5), size=1, prob=ProbSelTIDC)             
 
## Arrange the selected herds for visiting and those that have to be tested for testing
setQueue1.1 <- matrix(numeric(0),ncol=4)
setQueue2.1 <- matrix(numeric(0),ncol=4)
setQueue3.1 <- matrix(numeric(0),ncol=4)  
setQueue4.1 <- matrix(numeric(0),ncol=4)
setQueue5.1 <- matrix(numeric(0),ncol=4)     
setQueue6.1 <- matrix(numeric(0),ncol=4) 
## the categories in the last line of each matrix represent the type of visit 1=PV2, 2=PV1, 3=SV1, 4=SV2, 5=IDC, 6=DC.
## these categories are used to determine which test is used for the economic analysis of the alternative scenarios.
 if(sum(setQueue1)>0) setQueue1.1 <- cbind(which(setQueue1),gTime,1,1)
 if(sum(setQueue2)>0) setQueue2.1 <- cbind(which(setQueue2),gTime,SerSetQueue2,ifelse(SerSetQueue2==1,2,0))
 if(sum(setQueue3)>0) setQueue3.1 <- cbind(which(setQueue3),gTime,SerSetQueue3,ifelse(SerSetQueue3==1,3,0))
 if(sum(setQueue4)>0) setQueue4.1 <- cbind(which(setQueue4),gTime,SerSetQueue4,ifelse(SerSetQueue4==1,4,0))
 if(sum(setQueue5)>0) setQueue5.1 <- cbind(which(setQueue5),gTime,SerSetQueue5,ifelse(SerSetQueue5==1,5,0))
 if(sum(setQueue6)>0) setQueue6.1 <- cbind(which(setQueue6),gTime,1,6)



### Set herds that are in the zones or traced for queuing to surveillence and remove duplicates
## here we priortise the ones that will have PCR and Serology for sure.
         if(dim(setQueue1.1)[1]>0) zoneQueue <<- rbind(zoneQueue,setQueue1.1)
         if(dim(setQueue6.1)[1]>0) zoneQueue <<- rbind(zoneQueue,setQueue6.1)
         if(dim(setQueue5.1)[1]>0) zoneQueue <<- rbind(zoneQueue,setQueue5.1)
         if(dim(setQueue2.1)[1]>0) zoneQueue <<- rbind(zoneQueue,setQueue2.1)
         if(dim(setQueue3.1)[1]>0) zoneQueue <<- rbind(zoneQueue,setQueue3.1)
         if(dim(setQueue4.1)[1]>0) zoneQueue <<- rbind(zoneQueue,setQueue4.1)

         if(dim(zoneQueue)[1]>0){
           zoneQueue<<-zoneQueue[!duplicated(zoneQueue[,1]),,drop=FALSE]        
### Remove herds that have been or to be depopulated (in ring depopulation) from queuing to surveillance.
           tmp <- which(zoneQueue[,1]%in%depopQueue[,1] | aHerd$Diagnosed[zoneQueue[,1]] | aHerd$status[zoneQueue[,1]]==6)
           if(length(tmp)>0) zoneQueue<<-zoneQueue[-tmp,,drop=FALSE]  
           }

###### Calculate the total number of herds in queue today.
TotNumQueue <<- 0
TotNumQueue <<- dim(zoneQueue)[1]

### Move herds from queue to being Surveyed
   if(dim(zoneQueue)[1]>0){
     SurvMat <- zoneQueue
#if(iteration==2) print(SurvMat)
       
     survInd <- 1:dim(SurvMat)[1]  
       if(length(survInd) <= (CapSurvay)) 
          toSurv <- survInd
       else         
          toSurv <- survInd[1:sum(1:length(survInd)<=CapSurvay)]

     SurvMat2 <- SurvMat[toSurv,,drop=FALSE]
#if(iteration==2) print(SurvMat2)

### Survey herds: following clinical surveillance, herds maybe tested (serology+PCR) once the herd has doubled mortality (including sick animals) 
    IndexCl<-SurvMat2[SurvMat2[,3]==0,1]

    aHerd$VisitClInspect[IndexCl] <<- aHerd$VisitClInspect[IndexCl] + 1

    toBeCulled <- which((aHerd$ID%in%IndexCl) & !(aHerd$Diagnosed) & (aHerd$status==4))
                 if(length(toBeCulled)>0){
                    SickTime   <- aInfHerd$getTClic(toBeCulled)
                    ExpectMort <- (gTime-SickTime) * aHerd$ExpMortality[toBeCulled]
                    tmp        <- aHerd$Mortality[toBeCulled] >= (ExpectMort * MortalityIncreaseZone) & aHerd$Mortality[toBeCulled] >= NumDeadAnimSurv
                    toBeCulled <- toBeCulled[tmp]
                    if(length(toBeCulled)>0){
                    ### herds will be sampled and tested with PCR and serology
                      aHerd$sampVisitSer[toBeCulled] <<- aHerd$sampVisitSer[toBeCulled] + 1
                      aHerd$sampVisitPCR[toBeCulled] <<- aHerd$sampVisitPCR[toBeCulled] + 1
                      TMP1 <- aInfHerd$getInfected(toBeCulled)
                      TMP2 <- 1-(1-(aHerd$NumSamp[toBeCulled]/(aHerd$herdSize[toBeCulled]-((TMP1-1)/2))))^TMP1
                      TMP2[aHerd$NumSamp[toBeCulled]==aHerd$herdSize[toBeCulled]&TMP1>0] <-  1
                      TMP2[TMP2>1] <- 1               
                      toBeCulled <- toBeCulled[runif(length(toBeCulled))<=TMP2] 
                      if(length(toBeCulled)>0){
                        depopQueue <<-rbind(depopQueue,cbind(toBeCulled,gTime))
                        aHerd$Diagnosed[toBeCulled]     <<- TRUE
                        aHerd$DiagSurv[toBeCulled]      <<- TRUE
                        aHerd$diagnosisTime[toBeCulled] <<- gTime
                        aInfHerd$setDiagnosed(toBeCulled)

                      }# End of if(length....
                   }#End of if
                 }#End of if 

# This part includes detection of herds based on serology testing. Suspesion and testing is dependent on whether the herd is clinical or not.
     IndexSer <-SurvMat2[SurvMat2[,4]%in%SerologyTesting ,1]
     aHerd$sampVisitSer[IndexSer] <<- aHerd$sampVisitSer[IndexSer] + 1 
     toBeCulledSer   <- which((aHerd$ID%in%IndexSer) & !(aHerd$Diagnosed) & (aHerd$status!=4))
                 if(length(toBeCulledSer)>0){
                   TMP1 <- floor(aHerd$Survived[toBeCulledSer])
                   TMP2 <-  1-(1-(aHerd$NumSamp[toBeCulledSer]/(aHerd$herdSize[toBeCulledSer]-((TMP1-1)/2))))^TMP1
                   TMP2[aHerd$NumSamp[toBeCulledSer]==aHerd$herdSize[toBeCulledSer]&TMP1>0] <-  1
                   TMP2[TMP2>1] <- 1                   
                   toBeCulledSer<- toBeCulledSer[runif(length(toBeCulledSer))<=TMP2] 
                   if(length(toBeCulledSer)>0){
                    depopQueue <<-rbind(depopQueue,cbind(toBeCulledSer,gTime))
                    aHerd$Diagnosed[toBeCulledSer]     <<- TRUE
                    aHerd$DiagSurv[toBeCulledSer]      <<- TRUE
                    aHerd$diagnosisTime[toBeCulledSer] <<- gTime
                    aInfHerd$setDiagnosed(toBeCulledSer)
                   }# End of if
                  }#End of if
#### detection of herds following surveillance visit (not including tracing. it comes later). If the herd is clinical and fits the criteria for suspesion for ASF, the herd will be tested by PCR and serology.
#### If there is no suspesion, the herd is tested only by serology and hence probability of detection is based only on survived animals following infection.
     IndexSer2 <-SurvMat2[SurvMat2[,4]%in%SerologyTesting ,1]
     toBeCulledSerCl <- which((aHerd$ID%in%IndexSer2) & !(aHerd$Diagnosed) & (aHerd$status==4))
     TestedSerClPCR <- integer(0)
                 if(length(toBeCulledSerCl)>0){
                    SickTime   <- aInfHerd$getTClic(toBeCulledSerCl)
                    ExpectMort <- (gTime-SickTime) * aHerd$ExpMortality[toBeCulledSerCl]
                    tmp        <- aHerd$Mortality[toBeCulledSerCl] >= (ExpectMort * MortalityIncreaseZone) & aHerd$Mortality[toBeCulledSerCl] >= NumDeadAnimSurv
                    toBeCulledSerClPCR<- toBeCulledSerCl[tmp]
                    ## The herds that will be tested by PCR.
                    TestedSerClPCR <- toBeCulledSerClPCR
                    if(length(toBeCulledSerClPCR)>0){
                      aHerd$sampVisitPCR[toBeCulledSerClPCR] <<- aHerd$sampVisitPCR[toBeCulledSerClPCR] + 1 
                      TMP1 <- aInfHerd$getInfected(toBeCulledSerClPCR)
                      TMP2 <- 1-(1-(aHerd$NumSamp[toBeCulledSerClPCR]/(aHerd$herdSize[toBeCulledSerClPCR]-((TMP1-1)/2))))^TMP1
                      TMP2[aHerd$NumSamp[toBeCulledSerClPCR]==aHerd$herdSize[toBeCulledSerClPCR]&TMP1>0] <-  1
                      TMP2[TMP2>1] <- 1               
                      toBeCulledSerClPCR<- toBeCulledSerClPCR[runif(length(toBeCulledSerClPCR))<=TMP2]
                     }
                    toBeCulledSerCl2<- toBeCulledSerCl[!tmp]
                    if(length(toBeCulledSerCl2)>0){
                      TMP1 <- floor(aHerd$Survived[toBeCulledSerCl2])
                      TMP2 <-  1-(1-(aHerd$NumSamp[toBeCulledSerCl2]/(aHerd$herdSize[toBeCulledSerCl2]-((TMP1-1)/2))))^TMP1
                      TMP2[aHerd$NumSamp[toBeCulledSerCl2]==aHerd$herdSize[toBeCulledSerCl2]&TMP1>0] <-  1
                      TMP2[TMP2>1] <- 1                   
                      toBeCulledSerCl2<- toBeCulledSerCl2[runif(length(toBeCulledSerCl2))<=TMP2]
                      }
                   toBeCulledSerCl3 <- unique(c(toBeCulledSerClPCR,toBeCulledSerCl2))
                   if(length(toBeCulledSerCl3)>0){
                    depopQueue <<-rbind(depopQueue,cbind(toBeCulledSerCl3,gTime))
                    aHerd$Diagnosed[toBeCulledSerCl3]     <<- TRUE
                    aHerd$DiagSurv[toBeCulledSerCl3]      <<- TRUE
                    aHerd$diagnosisTime[toBeCulledSerCl3] <<- gTime
                    aInfHerd$setDiagnosed(toBeCulledSerCl3)
                   }# End of if
                  }#End of if

# This part includes detection of herds based on PCR testing. Here we still allow that the herd can be tested by PCR, if not tested
# during the previous step. IF the herd has been tested within a day, then the herd should not be tested again during the same day
# by the same test regardless the reason for testing.
    IndexPCR <-SurvMat2[SurvMat2[,4]%in%PCRTesting & !SurvMat2[,1]%in%TestedSerClPCR,1]
    aHerd$sampVisitPCR[IndexPCR] <<- aHerd$sampVisitPCR[IndexPCR] + 1 
    toBeCulledPCR <- which((aHerd$ID%in%IndexPCR) & !(aHerd$Diagnosed)  & (aHerd$status%in%c(3,4))) 
                 if(length(toBeCulledPCR)>0){
                   TMP1 <- aInfHerd$getInfected(toBeCulledPCR)
                   TMP2 <-  1-(1-(aHerd$NumSamp[toBeCulledPCR]/(aHerd$herdSize[toBeCulledPCR]-((TMP1-1)/2))))^TMP1
                   TMP2[aHerd$NumSamp[toBeCulledPCR]==aHerd$herdSize[toBeCulledPCR]&TMP1>0] <-  1
                   TMP2[TMP2>1] <- 1                   
                   toBeCulledPCR<- toBeCulledPCR[runif(length(toBeCulledPCR))<=TMP2] 
                   if(length(toBeCulledPCR)>0){
                    depopQueue <<-rbind(depopQueue,cbind(toBeCulledPCR,gTime))
                    aHerd$Diagnosed[toBeCulledPCR]     <<- TRUE
                    aHerd$DiagSurv[toBeCulledPCR]      <<- TRUE
                    aHerd$diagnosisTime[toBeCulledPCR] <<- gTime
                    aInfHerd$setDiagnosed(toBeCulledPCR)
                   }# End of if
                  }#End of if


####################################################################################
### this is a matrix to output the day the herd was sat for surveillance and the ###
### day the herd was actually surveyed (gTime in this case)                      ###
####################################################################################
  if(Detailed){
     if(sum(SurvMat2[,3]==0)>0)                 ClSurvMatOut  <<- rbind(ClSurvMatOut,cbind(iteration,gTime,SurvMat2[SurvMat2[,3]==0,1]))
     if(sum(SurvMat2[,4]%in%SerologyTesting)>0) SerSurvMatOut <<- rbind(SerSurvMatOut,cbind(iteration,gTime,SurvMat2[ SurvMat2[,4]%in%SerologyTesting ,1]))
     if(sum(SurvMat2[,4]%in%PCRTesting)>0)      PCRSurvMatOut <<- rbind(PCRSurvMatOut,cbind(iteration,gTime,unique(c(TestedSerClPCR,IndexPCR))))

      if(dim(ClSurvMatOut)[1]>= DumpData){
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
         NAMEClSH <- paste(runID,"ClSurvayedHerds.txt",sep="-")
         write.table(ClSurvMatOut,NAMEClSH,append=TRUE,col.names = F,row.names = F)
         ClSurvMatOut <<- matrix(numeric(0),ncol=3)
         }
      if(dim(SerSurvMatOut)[1]>= DumpData){
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
         NAMESerSH <- paste(runID,"SerSurvayedHerds.txt",sep="-")
         write.table(SerSurvMatOut,NAMESerSH,append=TRUE,col.names = F,row.names = F)
         SerSurvMatOut <<- matrix(numeric(0),ncol=3)
         }
     if(dim(PCRSurvMatOut)[1]>= DumpData){
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
         NAMEPCRSH <- paste(runID,"PCRSurvayedHerds.txt",sep="-")
         write.table(PCRSurvMatOut,NAMEPCRSH,append=TRUE,col.names = F,row.names = F)
         PCRSurvMatOut <<- matrix(numeric(0),ncol=3)
         }
        }

         aHerd$visitCount[SurvMat2[,1]]  <<- aHerd$visitCount[SurvMat2[,1]] + 1
         aHerd$timeVisited[SurvMat2[,1]] <<- gTime
         zoneQueue                       <<- zoneQueue[-toSurv,,drop=FALSE]

 
## extend the duration of the protection and surveillance zones in case herds are queuing for longer than the duration of each of the zones 
     if( any((gTime-zoneQueue[,2]) > (PZoneDuration-(SecProtVisit+DelayStartVisit))) ) PZoneDuration <<- PZoneDuration+ max(gTime-zoneQueue[,2])
           
     if( any((gTime-zoneQueue[,2]) > (SZoneDuration-(SecSurVisit+DelayStartVisit))) )  SZoneDuration <<- SZoneDuration+ max(gTime-zoneQueue[,2])
           

     }#End of if(dim(zoneQueue)[1]>0) 
 
    } # end of day
     ,
       getLabel=function(){ }
       )##EndOf list

  } #end of


######################################################
## Implementation of global control measures (Standstill)
## .
##
## The relative effects for DC and IC are stored in aHerd$relDC
## and aHerd$relIC respectively.
## Function adjusted for the Danish Project (LEC and TH)
###################################################

controlAll<-function(effectDC,
                              #effectIMC
                                        #,effectILC,
                                                    label,duration=3){
  if (!is.expression(effectDC))
    effectDC<-parse(text=effectDC)
  #if (!is.expression(effectIMC))
  #  effectIMC<-parse(text=effectIMC)
#if (!is.expression(effectILC))
  #  effectILC<-parse(text=effectILC)
  effDC<-numeric(0)
  #effIMC<-numeric(0)
 # effILC<-numeric(0)
  
  list(
       init=function(){ ## Initialize next iteration
         effDC<<-1-eval(effectDC,list(n=gMaxHerds))
        # effIMC<<-1-eval(effectIMC,list(n=gMaxHerds))
        # effILC<<-1-eval(effectILC,list(n=gMaxHerds))
       }
       ,
       day=function(){ ## Daily update of variables
	   if(gDaysUntilBaseline+duration>gTime){
           ## Update aHerd$relDC, aHerd$relIMC and aHerd$relILC
           aHerd$relDC <<- pmin(aHerd$relDC,effDC)
          # aHerd$relIMC <<- pmin(aHerd$relIMC,effIMC)
          # aHerd$relILC <<- pmin(aHerd$relILC,effILC)
         }

       }
       ,
       getLabel=function(){ return(label) }
       )##EndOf list
}##EndOf controlCirc

##################################
##
## A function to trace herds that have been contacted 
## by infected herds with animal contact for serveillance
##
########################################

traceDC <- function(prob,probdetect=1,delay,tracetime,duration,label){ ### prob her includes the probability to report a contact and to actually conduct a visit.###
 # if (!is.expression(prob))
 #   prob<-parse(text=prob)
  if (!is.expression(delay))
   delay<-parse(text=delay)
  if (!is.expression(tracetime))
   tracetime<-parse(text=tracetime)
  
  delayVisit<-numeric(0)
  forwardTraceTime<-numeric(0)
  backwardTraceTime<-numeric(0)

  list(
    
      init=function(){ ## Initialize next iteration
         delayVisit <<- eval(delay,list(n=gMaxHerds))
         forwardTraceTime <<- eval(tracetime,list(n=gMaxHerds))
         backwardTraceTime <<- eval(tracetime,list(n=gMaxHerds))
       }
       ,
      
      day= function(){

     
### Forward tracing (trace herds that had recieved animals from infected herds (TH))
  
     forwardToday <- ((aHerd$diagnosisTime+forwardTraceTime)==gTime) 
      if(sum(forwardToday)>0){
         traced <- which(traceMatrix[,3]%in%which(forwardToday) & (traceMatrix[,4]==1) & ((traceMatrix[,1]+duration) > gTime) )
	   if (length(traced)>0)
            traced <- traced[runif(length(traced))<prob]
          if(length(traced)>0){
	      traced <- traced[ aHerd$timeToVisitTraceDC[ traceMatrix[traced,2] ] < gTime ] ## Not revisiting already scheduled visits.                                                                                                                                                ## do not include tracing IDC, because now they should be tested.                       
            if(length(traced)>0){
              aHerd$timeToVisitTraceDC[ traceMatrix[traced,2] ] <<- delayVisit[ traceMatrix[traced,2] ] + gTime
              aHerd$traceDC[ traceMatrix[traced,2] ]            <<- aHerd$traceDC[ traceMatrix[traced,2] ]+ 1 
    #if(Detailed){
    #  TraceDCMatOut<<- rbind(TraceDCMatOut,cbind(iteration,gTime,traceMatrix[traced,2]))
    #  if(dim(TraceDCMatOut)[1]>= DumpData){
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
    #     NAMETDC <- paste(runID,"TDCHerds.txt",sep="-")
    #     write.table(TraceDCMatOut,NAMETDC,append=TRUE,col.names = F,row.names = F)
    #     TraceDCMatOut<<- matrix(numeric(0),ncol=3)
    #    }
    #   }
         }#End of if(length..
       }#End of if(sum...
       }

### Backward tracing ( traced herds that had sent animals to diagnosed herds)
          
     backwardToday<-(aHerd$diagnosisTime+backwardTraceTime)==gTime
      
      if(sum(backwardToday)>0){
	  traced <- which(traceMatrix[,2]%in%which(backwardToday) & (traceMatrix[,4]==1) & ((traceMatrix[,1]+duration) > gTime) )
         if(length(traced)>0)
           traced <- traced[runif(length(traced))<prob]
          if(length(traced)>0){
   	     traced <- traced[ aHerd$timeToVisitTraceDC[ traceMatrix[traced,3] ] < gTime ] ## Not revisiting already scheduled visits.                                                                                                                                               ## do not include tracing IDC, because now they should be tested.                       
            if(length(traced)>0){
           aHerd$timeToVisitTraceDC[ traceMatrix[traced,3] ] <<- delayVisit[ traceMatrix[traced,3] ] + gTime
           aHerd$traceDC[ traceMatrix[traced,3] ] <<- aHerd$traceDC[ traceMatrix[traced,3] ]+ 1
      #if(Detailed){
      #TraceDCMatOut<<- rbind(TraceDCMatOut,cbind(iteration,gTime,traceMatrix[traced,3]))
      #if(dim(TraceDCMatOut)[1]>= DumpData){
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
      #   NAMETDC <- paste(runID,"TDCHerds.txt",sep="-")
      #   write.table(TraceDCMatOut,NAMETDC,append=TRUE,col.names = F,row.names = F)
      #   TraceDCMatOut<<- matrix(numeric(0),ncol=3)
      #  }
      # }
          }#End of if(length..
         }#End of if(sum... 
        }         

      } #End of day
      ,
       getLabel=function(){ return(label) }
       )##EndOf list
     }##EndOf traceDC    

    
     
##################################
##
## A function to trace herds that have been contacted 
## by infected herds for serveillance
##
########################################

traceIDC <- function(timetotrace,delayvisitMed,delayvisitLow,duration,label){
 if (!is.expression(timetotrace))
   timetotrace<-parse(text=timetotrace)
 if (!is.expression(delayvisitMed))

   delayvisitMed<-parse(text=delayvisitMed)
 if (!is.expression(delayvisitLow))
   delayvisitLow<-parse(text=delayvisitLow)


  delay <- matrix(0,ncol=5)
  forwardTraceTime<-numeric(0)
  backwardTraceTime<-numeric(0)

  list(
    
      init=function(){ ## Initialize next iteration
         #probSelectTIDC <<- c(0,0.838,0.2,0.125)# the probability that a movement will not be forgetten and it will be traced and visited,
         forwardTraceTime<<- eval(timetotrace,list(n=gMaxHerds))
         backwardTraceTime<<- eval(timetotrace,list(n=gMaxHerds))
         delay <<- cbind(0,eval(delayvisitMed,list(n=gMaxHerds)),eval(delayvisitMed,list(n=gMaxHerds)),eval(delayvisitMed,list(n=gMaxHerds)),eval(delayvisitLow,list(n=gMaxHerds)))
       }
       ,
      day= function(){


### Forward tracing
          
     forwardToday <- (aHerd$diagnosisTime+forwardTraceTime)==gTime

      if(sum(forwardToday)>0){
	    traced <- which(traceMatrix[,3]%in%which(forwardToday) & (traceMatrix[,4]%in%ToTracedIDC) & ((traceMatrix[,1]+duration) > gTime) )
        if(length(traced)>0)
          traced <- traced[runif(length(traced))<probSelectTIDC[traceMatrix[traced,4]]]
         if(length(traced)>0){
   	     traced <- traced[ aHerd$timeToVisitTraceDC[ traceMatrix[traced,2] ] < gTime & aHerd$timeToVisitTraceIDC[ traceMatrix[traced,2] ] < gTime ] ## Not revisiting already scheduled visits.            
          if(length(traced)>0){
             tmp<-numeric(0)
             for(i in traced){
                tmp <- c(tmp,delay[traceMatrix[i,2],traceMatrix[i,4]])
                }
             aHerd$timeToVisitTraceIDC[ traceMatrix[traced,2] ] <<- tmp + gTime
             aHerd$traceIDC[ traceMatrix[traced,2] ] <<- aHerd$traceIDC[ traceMatrix[traced,2] ] + 1
    #if(Detailed){
    #  TraceIDCMatOut<<- rbind(TraceIDCMatOut,cbind(iteration,gTime,traceMatrix[traced,2]))
    #  if(dim(TraceIDCMatOut)[1]>= DumpData){
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
    #     NAMETIDC <- paste(runID,"TIDCHerds.txt",sep="-")
    #     write.table(TraceIDCMatOut,NAMETIDC,append=TRUE,col.names = F,row.names = F)
    #     TraceIDCMatOut<<- matrix(numeric(0),ncol=3)
    #    }
    #   }
          }#End of if(length
         }#End of if  
        }
   
### Backward tracing   
          
     backwardToday<-(aHerd$diagnosisTime+backwardTraceTime)==gTime
       
      if(sum(backwardToday)>0){
          traced <- which(traceMatrix[,2]%in%which(backwardToday) & (traceMatrix[,4]%in%ToTracedIDC) & ((traceMatrix[,1]+duration) > gTime) )
	  if(length(traced)>0)
   	    traced <- traced[runif(length(traced))<probSelectTIDC[traceMatrix[traced,4]]]
         if(length(traced)>0){
   	     traced <- traced[ aHerd$timeToVisitTraceDC[ traceMatrix[traced,3] ] < gTime & aHerd$timeToVisitTraceIDC[ traceMatrix[traced,3] ] < gTime ]## Not revisiting already scheduled visits.                         
          if(length(traced)>0){
             tmp<-numeric(0)
             for(i in traced){
                tmp <- c(tmp,delay[traceMatrix[i,3],traceMatrix[i,4]])
                }
            aHerd$timeToVisitTraceIDC[ traceMatrix[traced,3] ] <<- tmp + gTime 
            aHerd$traceIDC[ traceMatrix[traced,3] ] <<- aHerd$traceIDC[ traceMatrix[traced,3] ] + 1
   # if(Detailed){
   #   TraceIDCMatOut<<- rbind(TraceIDCMatOut,cbind(iteration,gTime,traceMatrix[traced,3]))
   #   if(dim(TraceIDCMatOut)[1]>= DumpData){
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
   #      NAMETIDC <- paste(runID,"TIDCHerds.txt",sep="-")
   #      write.table(TraceIDCMatOut,NAMETIDC,append=TRUE,col.names = F,row.names = F)
   #      TraceIDCMatOut <<- matrix(numeric(0),ncol=3)
   #     }
   #    }
            }#End of if(length...
           }#End of if(sum...
          }
         } #End of day
      ,
       getLabel=function(){ return(label) }
       )##EndOf list
     }##EndOf traceCont        

###########################################################
##
## A function to control movements from diagnosed herds
##
###########################################################

controlDiag<-function(effectDC,effectIMC,effectILC,label,duration=Inf){
  
  if (!is.expression(effectIMC))
    effectIMC<-parse(text=effectIMC)
  if (!is.expression(effectILC))
    effectILC<-parse(text=effectILC)
  effDC<-numeric(0)
  effIMC<-numeric(0)
  effILC<-numeric(0)
  
  list(
       init=function(){ ## Initialize next iteration
         effDC<<-1-rep(effectDC,gMaxHerds)
         effIMC<<-1-eval(effectIMC,list(n=gMaxHerds))
         effILC<<-1-eval(effectILC,list(n=gMaxHerds))
       }
       ,
       day=function(){ 
         ## First determine which herds were diagnosed

         diag<-aHerd$Diagnosed
           
         ## Second: Update aHerd$relDC and aHerd$relIC
         aHerd$relDC[diag]<<-pmin(aHerd$relDC[diag],effDC[diag])
         aHerd$relIMC[diag]<<-pmin(aHerd$relIMC[diag],effIMC[diag])
         aHerd$relILC[diag]<<-pmin(aHerd$relILC[diag],effILC[diag])

       }
       ,
       getIn=function(){  return(diag) }
       ,
       getLabel=function(){ return(label) }
       )##EndOf list
}##EndOf controlCirc



