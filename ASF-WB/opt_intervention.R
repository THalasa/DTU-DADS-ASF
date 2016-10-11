##### OPT_INTERVENTION library version 0.15.1 
##### Includes the following functions:
#### DummyInter, SurvDead and CullRing (TH)

#################################################
#################################################

## Dummy intervention
DummyInter<-function(){
  list(
       day=function(){
       },
       cleaniteration=function(){
       }
       )##EndOf list of functions
}##EndOf DummyInter


################################
### Survey Dead animals (TH) ###
################################
 SurvDead <- function(Zone=1){

    list(
       init=function(){}
       ,
       day=
        switch(Zone,
          {function(){
          ### On daily basis, check which herds to be set for surveillance of dead animals.
          Temp <- (aHerd$inProtZone) & !aHerd$Diagnosed & !aHerd$status%in%c(5,6) & aHerd$ToSurvDead<gTime
          if(sum(Temp)>0) aHerd$ToSurvDead[Temp] <<- gTime + DaysSurDead
          
          toBeCulled <- which(aHerd$SubDeadSamp == gTime)
          if(length(toBeCulled)>0){
            depopQueue <<-rbind(depopQueue,cbind(toBeCulled,gTime))
            aHerd$Diagnosed[toBeCulled]     <<- TRUE
            aHerd$DiagSurvDead[toBeCulled]  <<- TRUE
            aHerd$diagnosisTime[toBeCulled] <<- gTime
            aInfHerd$setDiagnosed(toBeCulled)
           }#End of if 

          ### Check which herds to be survayed today, given that they still are in a zone.
          SurvToday    <- which((aHerd$inProtZone) & !aHerd$Diagnosed & (aHerd$status<=4) & aHerd$ToSurvDead==gTime)

          if(length(SurvToday)>0){

 # if(Detailed){
 #     SurDeadMatOut<<- rbind(SurDeadMatOut,cbind(iteration,gTime,SurvToday))
 #     if(dim(SurDeadMatOut)[1]>= DumpData){
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
 #        NAMESD <- paste(runID,"SurDeadHerds.txt",sep="-")
 #        write.table(SurDeadMatOut,NAMESD,append=TRUE,col.names = F,row.names = F)
 #        SurDeadMatOut<<- matrix(numeric(0),ncol=3)
 #        }
 #       }
           ## Update tested herds
           aHerd$SampledDead[SurvToday]  <<- aHerd$SampledDead[SurvToday]  + 1

           Deaths      <- round(aHerd$ExpMortality[SurvToday] * DaysSurDead)
           Tmp         <- match(infHerdNums,SurvToday)

           Tmp         <- Tmp[!is.na(Tmp)]
           ### Get the infected herds that will be survayed today
           InfToSurvay <- infHerdNums[infHerdNums%in%SurvToday]
           if(length(InfToSurvay)>0){
             DeadOfASF    <-  aInfHerd$getDead(InfToSurvay)
             Deaths[Tmp]  <-  Deaths[Tmp] + DeadOfASF
             numOfSamples <-  ifelse(Deaths<=aHerd$NumSampDead[SurvToday],Deaths,aHerd$NumSampDead[SurvToday])
             probDet      <-  (1-(1-(numOfSamples[Tmp]/(Deaths[Tmp]-((DeadOfASF-1)/2))))^DeadOfASF)
             probDet[numOfSamples[Tmp]==Deaths[Tmp]&DeadOfASF>0] <- 1
             probDet[probDet>1] <- 1
             Submitted     <-  InfToSurvay[runif(length(InfToSurvay))<=probDet]
             aHerd$SubDeadSamp[Submitted]  <<- gTime + DelaySubDeadSamp
             aHerd$DeadSampTest[SurvToday] <<- aHerd$DeadSampTest[SurvToday] + numOfSamples
             }#End of if
             ### Update the number of tested animals
             if(length(InfToSurvay)==0){
                numOfSamples <-  ifelse(Deaths<=aHerd$NumSampDead[SurvToday],Deaths,aHerd$NumSampDead[SurvToday])
                aHerd$DeadSampTest[SurvToday] <<- aHerd$DeadSampTest[SurvToday] + numOfSamples
             }#End of if
            }#End of if      
           }#End of function
          }#End of first part
           ,
          {function(){
          ### On daily basis, check which herds to be set for surveillance of dead animals.
          Temp <- (aHerd$inProtZone|aHerd$inSurZone) & !aHerd$Diagnosed & aHerd$status<=4 & aHerd$ToSurvDead<gTime
          if(sum(Temp)>0) aHerd$ToSurvDead[Temp] <<- gTime + DaysSurDead

          toBeCulled <- which(aHerd$SubDeadSamp == gTime)
          if(length(toBeCulled)>0){
            depopQueue <<-rbind(depopQueue,cbind(toBeCulled,gTime))
            aHerd$Diagnosed[toBeCulled]     <<- TRUE
            aHerd$DiagSurvDead[toBeCulled]  <<- TRUE
            aHerd$diagnosisTime[toBeCulled] <<- gTime
            aInfHerd$setDiagnosed(toBeCulled)
           }#End of if 

          ### Check which herds to be survayed today, given that they still are in a zone.
          SurvToday    <- which((aHerd$inProtZone|aHerd$inSurZone) & !aHerd$Diagnosed & !aHerd$status%in%c(5,6) & aHerd$ToSurvDead==gTime)
           if(length(SurvToday)>0){

 # if(Detailed){
 #     SurDeadMatOut<<- rbind(SurDeadMatOut,cbind(iteration,gTime,SurvToday))
 #     if(dim(SurDeadMatOut)[1]>= DumpData){
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
 #        NAMESD <- paste(runID,"SurDeadHerds.txt",sep="-")
 #        write.table(SurDeadMatOut,NAMESD,append=TRUE,col.names = F,row.names = F)
 #        SurDeadMatOut<<- matrix(numeric(0),ncol=3)
 #        }
 #       }
           ## Update tested herds
           aHerd$SampledDead[SurvToday]  <<- aHerd$SampledDead[SurvToday]  + 1

           ### Number of tested animals 
           Deaths      <- round(aHerd$ExpMortality[SurvToday] * DaysSurDead)
           Tmp         <- match(infHerdNums,SurvToday)
           Tmp         <- Tmp[!is.na(Tmp)]
           ### Get the infected herds that will be survayed today
           InfToSurvay <- infHerdNums[infHerdNums%in%SurvToday]
           if(length(InfToSurvay)>0){
             DeadOfASF    <-  aInfHerd$getDead(InfToSurvay)
             Deaths[Tmp]  <-  Deaths[Tmp] + DeadOfASF
             numOfSamples <-  ifelse(Deaths<=aHerd$NumSampDead[SurvToday],Deaths,aHerd$NumSampDead[SurvToday])
             probDet      <-  (1-(1-(numOfSamples[Tmp]/(Deaths[Tmp]-((DeadOfASF-1)/2))))^DeadOfASF)
             probDet[numOfSamples[Tmp]==Deaths[Tmp]&DeadOfASF>0] <- 1
             probDet[probDet>1] <- 1
             Submitted     <-  InfToSurvay[runif(length(InfToSurvay))<=probDet]
             aHerd$SubDeadSamp[Submitted] <<- gTime + DelaySubDeadSamp
             aHerd$DeadSampTest[SurvToday] <<- aHerd$DeadSampTest[SurvToday] + numOfSamples
             }#End of if
             ### Update the number of tested animals
             if(length(InfToSurvay)==0){
                numOfSamples <-  ifelse(Deaths<=aHerd$NumSampDead[SurvToday],Deaths,aHerd$NumSampDead[SurvToday])
                aHerd$DeadSampTest[SurvToday] <<- aHerd$DeadSampTest[SurvToday] + numOfSamples
             }#End of if
            }#End of if      
           }#End of function
          }#End of second part

         )#End of switch
       ,
       cleaniteration=function(){
       }
       )#End of list
    }#End of SurvDead

######################################
## Ring Culling Danish project (TH) ##
######################################
CullRing <- function(size,startCulling,type=1){
  CullRingSize <- size
    list(
       init=function(){}
       ,
       day=
        switch(type,
        { function(){
         if(gTime >= (startCulling+gDaysUntilBaseline)){
          if (length(infHerdNums)>0) {
           tmpID<-aInfHerd$getDiagnosedIDs()
           tmpR<-rowSums(as.matrix(Dist$get(tmpID)<= CullRingSize))
           Rdepops <-tmpR & !aHerd$status%in%c(5,6) & !(aHerd$Diagnosed) & aHerd$cullEligible 
           if (sum(Rdepops)>0){
            aHerd$status[Rdepops] <<-6
            aHerd$timeToTaggedForDepop[Rdepops] <<- gTime
            aHerd$timeToPV1[Rdepops] <<- 0
            aHerd$timeToPV2[Rdepops] <<- 0
            aHerd$timeToSV1[Rdepops] <<- 0
            aHerd$timeToSV2[Rdepops] <<- 0
            aHerd$timeToVisitTraceIDC[Rdepops] <<- 0
            aHerd$timeToVisitTraceDC[Rdepops]  <<- 0
            depopQueue<<-rbind(depopQueue,cbind(which(Rdepops),gTime))
  if(Detailed){
      PreEmpMatOut<<- rbind(PreEmpMatOut,cbind(iteration,gTime,which(Rdepops)))
      if(dim(PreEmpMatOut)[1]>= DumpData){
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
         NAMED <- paste(runID,"PreEmpHerds.txt",sep="-")
         write.table(PreEmpMatOut,NAMED,append=TRUE,col.names = F,row.names = F)
         PreEmpMatOut<<- matrix(numeric(0),ncol=3)
         }
        }
       }
            }  ##EndOf if
           }
          }
         },
        { function(){
          if(sum(aHerd$Diagnosed) >= startCulling){
          if (length(infHerdNums)>0) {
           tmpID<-aInfHerd$getDiagnosedIDs()
           tmpR<-rowSums(as.matrix(Dist$get(tmpID)<= CullRingSize))
           Rdepops <-tmpR & !aHerd$status%in%c(5,6) & !(aHerd$Diagnosed) & aHerd$cullEligible 
           if (sum(Rdepops)>0){
            aHerd$status[Rdepops] <<-6
            aHerd$timeToTaggedForDepop[Rdepops] <<- gTime
            aHerd$timeToPV1[Rdepops] <<- 0
            aHerd$timeToPV2[Rdepops] <<- 0
            aHerd$timeToSV1[Rdepops] <<- 0
            aHerd$timeToSV2[Rdepops] <<- 0
            aHerd$timeToVisitTraceIDC[Rdepops] <<- 0
            aHerd$timeToVisitTraceDC[Rdepops]  <<- 0
            depopQueue<<-rbind(depopQueue,cbind(which(Rdepops),gTime))
  if(Detailed){
      PreEmpMatOut<<- rbind(PreEmpMatOut,cbind(iteration,gTime,which(Rdepops)))
      if(dim(PreEmpMatOut)[1]>= DumpData){
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
         NAMED <- paste(runID,"PreEmpHerds.txt",sep="-")
         write.table(PreEmpMatOut,NAMED,append=TRUE,col.names = F,row.names = F)
         PreEmpMatOut<<- matrix(numeric(0),ncol=3)
         }
        }
       }
            }
           }  
          } 
         })
       ,
       cleaniteration=function(){
       }
       )##EndOf list of functions
} ##EndOf  ring culling
