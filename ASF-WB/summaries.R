#### summary library version 0.15.1
#### 
####
#### 
sumTh<-function(step){
  switch(step,
         init={if (verbose) cat("Initializing...\n")
               ## initialize distance calculations
               eval.parent(expression(aHerd$dDist <- NULL))
               ## global count of herds depopulated
               eval.parent(expression(gDepopulationCount <- rep(0,n)))
               ## global count of vaccinated herds
               eval.parent(expression(gVaccineCount <- rep(0,n)))
               ## number of times each herd was infected
               eval.parent(expression(gEpiDur <- rep(NA,n)))
               eval.parent(expression(ObsEpiDur <- rep(NA,n)))
               eval.parent(expression(NumInf <- rep(0,n)))
               eval.parent(expression(NumDet <- rep(0,n)))
               eval.parent(expression(NumCulled <- rep(0,n)))
               eval.parent(expression(Recovered<- rep(0,n)))
               eval.parent(expression(NumCullAnim <- rep(0,n)))
               eval.parent(expression(VisEpi <- rep(0,n)))
               eval.parent(expression(ClVisCount <- rep(NA,n)))
               eval.parent(expression(PCRVisHerdCount <- rep(NA,n)))
               eval.parent(expression(SerVisHerdCount <- rep(NA,n)))
               eval.parent(expression(PCRVisAnimCount <- rep(NA,n)))
               eval.parent(expression(SerVisAnimCount <- rep(NA,n)))
               eval.parent(expression(traceVisCount <- rep(0,n)))
               eval.parent(expression(VisitCostsTes <- rep(NA,n)))
               eval.parent(expression(VisitCostsCl <- rep(NA,n)))
               eval.parent(expression(logCostsSw <- rep(NA,n)))
               eval.parent(expression(CCostsPigFS <- rep(NA,n)))
               eval.parent(expression(CCostsPigInd <- rep(NA,n)))
               eval.parent(expression(ESCosts <- rep(NA,n)))
               eval.parent(expression(CosSS <- rep(NA,n)))
               eval.parent(expression(CompCosts <- rep(NA,n)))
               eval.parent(expression(LoExpLivSwEU<- rep(NA,n)))
               eval.parent(expression(LoExpSwProdEU<- rep(NA,n)))
               eval.parent(expression(LoExpOutZoneEU<- rep(NA,n)))
               eval.parent(expression(LoExpSwNEU<- rep(NA,n)))
               eval.parent(expression(cosWSFin <- rep(NA,n)))
               eval.parent(expression(cosWSWean <- rep(NA,n)))
               eval.parent(expression(MeetDelay <- rep(0,n)))
               eval.parent(expression(EUblock <- rep(0,n)))
               eval.parent(expression(redPrEU <- rep(0,n)))
               eval.parent(expression(VisAftEpEnd <- rep(0,n)))
               eval.parent(expression(ZoneDur <- rep(0,n)))
               eval.parent(expression(FirstEpiDet<- rep(0,n)))
               eval.parent(expression(diagHerdsSurv<- rep(0,n)))
               eval.parent(expression(SurvCosDead<- rep(0,n)))
               eval.parent(expression(NumSurvDeadAnim<- rep(0,n)))
               eval.parent(expression(NumSurvDeadHerd<- rep(0,n)))
               eval.parent(expression(DetFromSurvDead<- rep(0,n)))
               eval.parent(expression(InfectedFromWB<- rep(0,n)))


             },
         day={
             

         },
         iter={
           if(verbose) cat("\nRecording iteration Values\n")
           gEpiDur[iteration]         <<- gTime - gDaysUntilBaseline
           ObsEpiDur[iteration]       <<- max(aHerd$timeCulled) - gDaysUntilBaseline
           NumDet[iteration]          <<- sum(aHerd$Diagnosed)
           NumInf[iteration]          <<- sum(aHerd$Diagnosed | aHerd$SusAgain>0 | aHerd$infMode>0 | aHerd$status==7)
           NumCulled[iteration]       <<- sum(aHerd$status%in%c(5,6)) 
           Recovered[iteration]       <<- sum(aHerd$SusAgain>0 | aHerd$status==7)
           NumCullAnim[iteration]     <<- sum(aHerd$herdSize[aHerd$status%in%c(5,6)])
           NumSurvDeadAnim[iteration] <<- sum(aHerd$DeadSampTest)
           NumSurvDeadHerd[iteration] <<- sum(aHerd$SampledDead)
           DetFromSurvDead[iteration] <<- sum(aHerd$DiagSurvDead)
           InfectedFromWB[iteration]  <<- sum(aHerd$infMode%in%6:7)
        ### Herds that are in the zones and must be visited after the model run is finished (last detected herd is culled).
        ### Notice that we do not include the possibility of testing after with PCR for herds assigned to serology because of suspecion, 
        ### as we do in the survHerds function during the model run, because at this stage we know that there are no clinical herds. 
        ### The outbreak is finished and hence no suspesion will happen.
             if(2%in%SerologyTesting) aHerd$sampVisitSer[aHerd$timeToPV1>=gTime] <- aHerd$sampVisitSer[aHerd$timeToPV1>=gTime] + 1
             if(3%in%SerologyTesting) aHerd$sampVisitSer[aHerd$timeToSV1>=gTime] <- aHerd$sampVisitSer[aHerd$timeToSV1>=gTime] + 1
             if(4%in%SerologyTesting) aHerd$sampVisitSer[aHerd$timeToSV2>=gTime] <- aHerd$sampVisitSer[aHerd$timeToSV2>=gTime] + 1

             if(1%in%PCRTesting) aHerd$sampVisitPCR[aHerd$timeToPV2>=gTime] <- aHerd$sampVisitPCR[aHerd$timeToPV2>=gTime] + 1
             if(2%in%PCRTesting) aHerd$sampVisitPCR[aHerd$timeToPV1>=gTime] <- aHerd$sampVisitPCR[aHerd$timeToPV1>=gTime] + 1
             if(3%in%PCRTesting) aHerd$sampVisitPCR[aHerd$timeToSV1>=gTime] <- aHerd$sampVisitPCR[aHerd$timeToSV1>=gTime] + 1
             if(4%in%PCRTesting) aHerd$sampVisitPCR[aHerd$timeToSV2>=gTime] <- aHerd$sampVisitPCR[aHerd$timeToSV2>=gTime] + 1  

             if(!(4%in%SerologyTesting)) aHerd$VisitClInspect[aHerd$timeToSV2>=gTime] <- aHerd$VisitClInspect[aHerd$timeToSV2>=gTime] + 1
             if(!(2%in%PCRTesting)) aHerd$VisitClInspect[aHerd$timeToPV1>=gTime]      <- aHerd$VisitClInspect[aHerd$timeToPV1>=gTime] + 1
             if(firstSurvVisit&!(3%in%PCRTesting)) aHerd$VisitClInspect[aHerd$timeToSV1>=gTime] <- aHerd$VisitClInspect[aHerd$timeToSV1>=gTime] + 1

             VisCountCl  <- sum(aHerd$VisitClInspect)
             ### assign the herds that will get the IDC tracing visit and sampled for serology and PCR
             TEMP <- rep(FALSE,length(aHerd$timeToVisitTraceIDC))
             Tmp <- which(aHerd$timeToVisitTraceIDC>=gTime)
             if(length(Tmp)>0){
                Tmp <- Tmp[runif(length(Tmp)) <= ProbSelTIDC]
                TEMP[Tmp] <- TRUE
             }
             
             VisCountSer <- aHerd$sampVisitSer + (aHerd$timeToPV2>=gTime) + (aHerd$timeToVisitTraceDC>=gTime) + TEMP
             VisCountPCR <- aHerd$sampVisitPCR + (aHerd$timeToVisitTraceDC>=gTime) + TEMP

             ClVisCount[iteration]       <<- VisCountCl 
             SerVisHerdCount[iteration]  <<- sum(VisCountSer) 
             PCRVisHerdCount[iteration]  <<- sum(VisCountPCR) 
             SerVisAnimCount[iteration]  <<- sum(VisCountSer * aHerd$NumSamp)
             PCRVisAnimCount[iteration]  <<- sum(VisCountPCR * aHerd$NumSamp) 

           diagHerdsSurv[iteration] <<- sum(aHerd$DiagSurv)     
           FirstEpiDet[iteration]   <<- gDaysUntilBaseline
           MeetDelay[iteration]     <<- meetDelay           

           ### in case of detailed information is demanded, then check the herds that will still have to be surveyed and att then to the
       ### relevant surveillance output data file
        if (Detailed){      
             ClinicalInsp<- rep(FALSE,gMaxHerds)
             if(!(4%in%SerologyTesting)) ClinicalInsp[aHerd$VisitClInspect>0&aHerd$timeToSV2>=gTime] <- TRUE
             if(!(2%in%PCRTesting))      ClinicalInsp[aHerd$VisitClInspect>0&aHerd$timeToPV1>=gTime] <- TRUE
             if(firstSurvVisit&!(3%in%PCRTesting)) ClinicalInsp[aHerd$VisitClInspect>0&aHerd$timeToSV1>=gTime] <- TRUE
     if(sum(ClinicalInsp)>0){
      ClSurvMatOut  <<- rbind(ClSurvMatOut,cbind(iteration,gTime,which(ClinicalInsp)))
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
         NAMEClSH <- paste(runID,"ClSurvayedHerds.txt",sep="-")
         write.table(ClSurvMatOut,NAMEClSH,append=TRUE,col.names = F,row.names = F)
         ClSurvMatOut <<- matrix(numeric(0),ncol=3)
         }

             SerTestingPV1<-rep(FALSE,gMaxHerds)
             SerTestingSV1<-rep(FALSE,gMaxHerds)
             SerTestingSV2<-rep(FALSE,gMaxHerds)
             SerTestingPV2<-rep(FALSE,gMaxHerds)
             SerTestingTDC<-rep(FALSE,gMaxHerds)
             SerTestingTIDC<-rep(FALSE,gMaxHerds)
             if(2%in%SerologyTesting) SerTestingPV1[aHerd$timeToPV1>=gTime] <- TRUE
             if(3%in%SerologyTesting) SerTestingSV1[aHerd$timeToSV1>=gTime] <- TRUE
             if(4%in%SerologyTesting) SerTestingSV2[aHerd$timeToSV2>=gTime] <- TRUE
             SerTestingPV2[aHerd$timeToPV2>=gTime] <- TRUE
             SerTestingTDC[aHerd$timeToVisitTraceDC>=gTime] <- TRUE
             SerTestingTIDC[TEMP] <- TRUE ## TEMP is the herds that will be traced through IDC. determined above.
      if(sum(SerTestingPV1|SerTestingSV1|SerTestingSV2|SerTestingPV2|SerTestingTDC|SerTestingTIDC)>0){
       SerTesting <- c(which(SerTestingPV1),which(SerTestingSV1),which(SerTestingSV2),which(SerTestingPV2),which(SerTestingTDC),which(SerTestingTIDC))       
       SerSurvMatOut <<- rbind(SerSurvMatOut,cbind(iteration,gTime,SerTesting))
### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
         NAMESerSH <- paste(runID,"SerSurvayedHerds.txt",sep="-")
         write.table(SerSurvMatOut,NAMESerSH,append=TRUE,col.names = F,row.names = F)
         SerSurvMatOut <<- matrix(numeric(0),ncol=3)
         }

             TestingWPCRPV1<-rep(FALSE,gMaxHerds)
             TestingWPCRSV1<-rep(FALSE,gMaxHerds)
             TestingWPCRPV2<-rep(FALSE,gMaxHerds)
             TestingWPCRSV2<-rep(FALSE,gMaxHerds)
             TestingWPCRTDC<-rep(FALSE,gMaxHerds)
             TestingWPCRTIDC<-rep(FALSE,gMaxHerds)
             if(2%in%PCRTesting) TestingWPCRPV1[aHerd$timeToPV1>=gTime] <- TRUE
             if(3%in%PCRTesting) TestingWPCRSV1[aHerd$timeToSV1>=gTime] <- TRUE
             if(4%in%PCRTesting) TestingWPCRSV2[aHerd$timeToSV2>=gTime] <- TRUE
             if(1%in%PCRTesting) TestingWPCRPV2[aHerd$timeToPV2>=gTime] <- TRUE
             TestingWPCRTDC[aHerd$timeToVisitTraceDC>=gTime] <- TRUE
             TestingWPCRTIDC[TEMP] <- TRUE

     if(sum(TestingWPCRPV1|TestingWPCRSV1|TestingWPCRSV2|TestingWPCRPV2|TestingWPCRTDC|TestingWPCRTIDC)>0){
      TestingWPCR <- c(which(TestingWPCRPV1),which(TestingWPCRSV1),which(TestingWPCRSV2),which(TestingWPCRPV2),which(TestingWPCRTDC),which(TestingWPCRTIDC))
      PCRSurvMatOut <<- rbind(PCRSurvMatOut,cbind(iteration,gTime,TestingWPCR))

### NAME here will be exactly the same as that in the initialization file, 
### so no worries; no overwriting will happen ;-) (TH)
         NAMEPCRSH <- paste(runID,"PCRSurvayedHerds.txt",sep="-")
         write.table(PCRSurvMatOut,NAMEPCRSH,append=TRUE,col.names = F,row.names = F)
         PCRSurvMatOut <<- matrix(numeric(0),ncol=3)
         }
        }#end of if(Detalied){                  

        ###########################################################
        ###########################################################
        #################### ECONOMIC EFFECTS #####################
        ###########################################################
        ###########################################################

            #################################################
            #################################################
            ################## DIRECT COSTS #################
            #################################################
            #################################################
      

##########Calculate costs of visits  
### we include as well herds that are queuing for visit in calculating the costs. the costs include the costs of clinical inspection for herds that will only be clinicaly inspected and the herds that
### will be samples, the costs of serology sampling and testing and the costs of PCR.        

            VisCosCl    <- sum(VisCountCl  * CostsVisFarSw) 
            VisCosSer   <- sum(VisCountSer * aHerd$NumSamp * CostsSamTesSer) + sum(VisCountSer * CostsVisFarSw)
            VisCosPCR   <- sum(VisCountPCR * aHerd$NumSamp * CostsSamTesPCR) + sum(VisCountPCR * CostsVisFarSw)

            #VisCosSer30 <- sum(VisCountSer * CostsSamTesSer30) + sum(VisCountSer * CostsVisFarSw) 
            #VisCosPCR30 <- sum(VisCountPCR * CostsSamTesPCR30) + sum(VisCountPCR * CostsVisFarSw)
            #VisCosSer60 <- sum(VisCountSer * CostsSamTesSer60) + sum(VisCountSer * CostsVisFarSw)
            #VisCosPCR60 <- sum(VisCountPCR * CostsSamTesPCR60) + sum(VisCountPCR * CostsVisFarSw)
            

           VisitCostsCl[iteration] <<- VisCosCl 
           VisitCostsTes[iteration] <<- VisCosSer + VisCosPCR
           SurvCosDead[iteration]  <<- (sum(aHerd$DeadSampTest) * (CostAnimPCR + CostAnimSer)) + (sum(aHerd$SampledDead)*(CosIndsendelse)) 

           #VisitCosts30[iteration] <<- VisCosSer30 + VisCosPCR30
           #VisitCosts60[iteration] <<- VisCosSer60 + VisCosPCR60
 ############# matrix to estimate several direct cost

            CCDMat <- cbind(aHerd$sows[aHerd$status%in%c(5,6)],aHerd$finishers[aHerd$status%in%c(5,6)],aHerd$daysInZones[aHerd$status%in%c(5,6)])                  

##########Calculate costs of logistics (including testing and distruction) 
            
           logCostsSw[iteration]    <<- dim(CCDMat)[1] * logCosSw 
 
##########Calculate costs of cleaning and disinfection

            CDMat <- cbind(aHerd$sows[aHerd$Diagnosed],aHerd$finishers[aHerd$Diagnosed])

            CDCosPig <- dim(CDMat)[1]  * clCosPig
           CCostsPigFS[iteration]  <<- CDCosPig 


##########Calculate costs of empty stables
            
            CCDMat[CCDMat[,3] < max(PZoneDuration,SZoneDuration),3] <- max(PZoneDuration,SZoneDuration)

            EmStCosSow  <- sum(CCDMat[,1] * CCDMat[,3] * CESPSow)
            EmStCosFin  <- sum(CCDMat[,2] * CCDMat[,3] * CESPFin)
        
           ESCosts[iteration] <<- EmStCosSow + EmStCosFin 
                                   

##########Calculate compensations
 
            CompSw     <- sum(CCDMat[,1] * comSows) + sum(CCDMat[,2]  * comFin)
            CompCosts[iteration] <<- CompSw 

#########Calculate the costs of welfare slaughter
          index <- aHerd$daysInZones>0
          WelSl <- cbind(aHerd$sows[index],aHerd$daysInZones[index])
          WelSl[WelSl[,2] < max(PZoneDuration,SZoneDuration),2] <- max(PZoneDuration,SZoneDuration)

          cosWSFin[iteration]  <<- sum((WelSl[,1] * FactWSFinishers) * (WelSl[,2]-14) * CosWSFinishers )
          cosWSWean[iteration] <<- sum((WelSl[,1] * FactWSWeaners) * (WelSl[,2]-14) * CosWSWeaners)

##########Calculate costs of 3 days standstill
           
           if(any(aHerd$timeCulled>0)) CosSS[iteration] <<- rpert(1,30000000,40000000,50000000)
           else                        CosSS[iteration] <<- 0

       
            ####################################################
            ####################################################
            ################### EXPORT LOSSES ##################
            ####################################################
            ####################################################

########## Export loss from swine production systems #############################################################################################
## Ban on export of live animals from the whole country and starting from detection of the first outbreak until optional days                   ##
##   after culling last infected herd in the country and lifting the zones. Meat product ban on farms within the zones during and a waiting     ##
##   period after the outbreak.                                                                                                                 ##
## Ban on export of live animals and swine products from the whole country from the detection of the first infection until the culling of the   ##
##   last infected herd and lifiting the zone including an optional delay and a delay until the meeting to regain the free status. Live animals  ##
##   and swine products (from outside the zones) are assumed to be sold in the EU market with a 25% reduced price.                               ##
##################################################################################################################################################
 
          SwinExMat <- cbind(aHerd$herdSize,aHerd$daysInZones)
          SwinExMat[SwinExMat[,2]>0 & SwinExMat[,2] < max(PZoneDuration,SZoneDuration),2] <- max(PZoneDuration,SZoneDuration)


            if(any(aHerd$timeCulled>0)){
               LoExpLivSwEU[iteration]  <<- ((expBanAftCEU + max(aHerd$timeCulled)+ max(PZoneDuration,SZoneDuration)- gDaysUntilBaseline + 1) * totExpLivSwEU)    # ZoneDuration - gDaysUntilBaseline + 1, ZoneDuration because it will take extra ZoneDuration 
               LoExpSwProdEU[iteration] <<- sum((SwinExMat[,1] / sum(SwinExMat[,1])) * (SwinExMat[,2]) * totExpSwProdEU) ## export loss from inside the zones                                                                                   
               LoExpSwNEU[iteration]    <<-  (expBanAftCNEU + max(aHerd$timeCulled) - gDaysUntilBaseline + 1 + meetDelay) * totExpSwNEU * reducedPrice 
              }   
                                                                                                      
            else {
               LoExpLivSwEU[iteration]  <<- 0  
               LoExpSwProdEU[iteration] <<- 0
               LoExpSwNEU[iteration]    <<- 0  
              } 

                #########################
                ### Producing results ###
                #########################
### Exporting the summary results ###
      SumResOut <<- cbind(FirstEpiDet[iteration],gEpiDur[iteration],ObsEpiDur[iteration],NumInf[iteration],InfectedFromWB[iteration],NumDet[iteration],NumCulled[iteration],diagHerdsSurv[iteration],
                DetFromSurvDead[iteration],Recovered[iteration],NumCullAnim[iteration],ClVisCount[iteration],SerVisHerdCount[iteration],SerVisAnimCount[iteration],
                PCRVisHerdCount[iteration],PCRVisAnimCount[iteration],NumSurvDeadHerd[iteration],NumSurvDeadAnim[iteration],VisitCostsCl[iteration],
                VisitCostsTes[iteration],SurvCosDead[iteration],logCostsSw[iteration],CCostsPigFS[iteration],ESCosts[iteration],cosWSFin[iteration],
                cosWSWean[iteration],CompCosts[iteration],CosSS[iteration],LoExpSwProdEU[iteration],LoExpLivSwEU[iteration],LoExpSwNEU[iteration],MeetDelay[iteration]) 
      NAME <- paste(runID,"ASF.txt",sep="-") 
      write.table(SumResOut,NAME,append=T,sep=" ",col.names = F,row.names=F)                                                     
      SumResOut  <<- matrix(numeric(0),ncol=32)

### Exporting the Infected herds matrix ### 
      IndexAIH <- which(aHerd$Diagnosed | aHerd$SusAgain>0 | aHerd$infMode>0 | aHerd$status==7)
      AllInfHerds <<- cbind(iteration,aHerd$diagnosisTime[IndexAIH],aHerd$immuneTime[IndexAIH],IndexAIH,
                      aHerd$timeInfected[IndexAIH],aHerd$infMode[IndexAIH],aHerd$infSource[IndexAIH]) 
      NAMEInf <- paste(runID,"AllInfHerds.txt",sep="-")
      write.table(AllInfHerds,NAMEInf,append=T,sep=" ",col.names = F,row.names=F)
      AllInfHerds <<- matrix(numeric(0),ncol=7)
                                                                                                                      
         if(Detailed){ 

           #NAMES    <- paste(runID,"SurvHerds.txt",sep="-")
           #NAMEP    <- paste(runID,"ProtHerds.txt",sep="-")
           #NAMETDC  <- paste(runID,"TDCHerds.txt",sep="-")
           #NAMETIDC <- paste(runID,"TIDCHerds.txt",sep="-")
           NAMED    <- paste(runID,"DepopHerds.txt",sep="-")
           NAMEPE   <- paste(runID,"PreEmpHerds.txt",sep="-")
          # NAMESD   <- paste(runID,"SurDeadHerds.txt",sep="-")
           #write.table(SurvZoneMatOut,NAMES,append=T,col.names = F,row.names=F)
           #write.table(ProtZoneMatOut,NAMEP,append=T,col.names = F,row.names=F)
           #write.table(TraceDCMatOut,NAMETDC,append=T,col.names = F,row.names=F)
           #write.table(TraceIDCMatOut,NAMETIDC,append=T,col.names = F,row.names=F)
           write.table(DepopMatOut,NAMED,append=T,col.names = F,row.names=F)
           write.table(PreEmpMatOut,NAMEPE,append=T,col.names = F,row.names=F)
           #write.table(SurDeadMatOut,NAMESD,append=T,col.names = F,row.names=F)

# SurvZoneMatOut   <<- matrix(numeric(0),ncol=3)
# ProtZoneMatOut   <<- matrix(numeric(0),ncol=3)
# TraceDCMatOut    <<- matrix(numeric(0),ncol=3)
# TraceIDCMatOut   <<- matrix(numeric(0),ncol=3)
 DepopMatOut      <<- matrix(numeric(0),ncol=3)
 PreEmpMatOut     <<- matrix(numeric(0),ncol=3)
 #SurDeadMatOut    <<- matrix(numeric(0),ncol=3)
         }                                                                                                                                                                
     cat(iteration," ")
         },
         final={
   

return(aHerd)                                                                                                    
          
        },
         ## This will be executed if nothing above matches!
         cat("Step should be either init, day, iter, or final")
       )
} ## End of sumTH


