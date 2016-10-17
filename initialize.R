### INITIALIZE library V. 0.15.1
### Includes the following functions:
###
### selectIndexHerd
### createASFvars
### initializeASFvars

######################
### selectIndexHerd
######################
## The options are used to select rows in aHerd.
## This function is called in the beginning of each simulated outbreak
## See help file in package for further details.
selectIndexHerd<-function(args=NULL,...){
  if (is.null(args))
    args<-list(...)

  #sel<-rep(TRUE,length(aHerd[[1]]))
  #Freq <- table(args)
  ## Allowing more than one index case in each outbreak:
  if (nTmp<-match("nIndex",names(args),nomatch=0)){
    nIndex<-args[[nTmp]]
    args<-args[-nTmp] ## removing "nIndex" from the list of arguments
  }
  else
    nIndex<-1
  ## Matching the remaining arguments:
  if (length(args)>0){
    indexAHerd<-match(names(args),names(aHerd))


    if (any(is.na(indexAHerd)))
      stop(paste("The argument(s):",
                 paste(names(args)[is.na(indexAHerd)],collapse=" and "),
                 "was/were not recognized.\n"))
  }

tmp <- args[[1]]
  tmp[iteration]
 }


######################
### selectIndexHerd 2
######################
## Start the infection in more than 1 index herds
selectIndHerdMore<-function(args=NULL){

#### args should be a data frame containing the index herds
#### Each row in args includes the multiple index herds (History herds)
#### that will be used to initiate the epidemic
    tmp <- args
    tmp <- tmp[iteration,]
    tmp <- unlist(unname(tmp))
    tmp <- tmp[!is.na(tmp)]
 }


##################################################################
### createASFvars
###
### This funtion is called once to read data and initialize variables
##################################################################
createASFvars <- function() {
  ## Load Data
  if (verbose) cat("Loading Data\n")

  ## read in herd information from comma delimited infofile
  eval.parent(expression(aHerd<-NULL))
  if (!is.null(infofile))
    if ((infofile=="pkg3county")&&("package:ASF"%in%search())){
      data(pkg3county,envir=environment())
      aHerd <<- as.list(pkg.three.county)
    }
    else
      aHerd <<- as.list(read.table(infofile,sep=";",header=T))

  names(aHerd)[1]<<-"ID" ## Making sure the ID column is named "ID"
  ## Check aHerd for non-unique IDs
  if(length(unique(aHerd$ID))<length(aHerd$ID))
    stop("Herd ID numbers are not unique. Simulations Fails",
         paste(unique(aHerd$ID[duplicated(aHerd$ID)]),collapse=", "))
  ## read in herd types from comma delimited typesfile
  eval.parent(expression(herdtypes<-NULL))
  if (!is.null(typesfile))
    if ((typesfile=="pkgtypes")&&("package:ASF"%in%search())){
      data(pkgtypes,envir=environment())
      herdtypes<<-as.list(pkgtypes)
    }
    else
      herdtypes<<-as.list(read.table(typesfile,sep=";",as.is=TRUE,header=TRUE))

  ## modify one column of herdtypes table if tornado plot option is selected
  if (!is.null(tornCol)) {
    if (is.numeric(herdtypes[[tornCol]]))
      herdtypes[[tornCol]] <<- 0.9 * herdtypes[[tornCol]]
    else
      herdtypes[[tornCol]] <<-
        paste(as.character(tornMult),"*",herdtypes[[tornCol]])
  }
  ## Parsing 'herdtypes' once and for all:
  for (i in 3:length(herdtypes))
    herdtypes[[i]]<<-parse(text=herdtypes[[i]])
  ##################################################
  ##
  ## New Changes for the Danish Model, TH and LEC
  ##
  ##################################################
  ### Read the movement and catrgory matrises for All and for weaners
   eval.parent(expression(MovMatAll<-NULL))
   eval.parent(expression(MovMatWean<-NULL))
   eval.parent(expression(MovAb<-NULL))

   MovMatAll<<-as.matrix(read.table(fileMovMatAll,sep=";",dec=",")) 
   MovMatWean<<-as.matrix(read.table(fileMovMatWean,sep=";",dec=","))
   MovAb<<-as.matrix(read.table(fileMovAb,sep=";",dec=",")) 

eval.parent(expression(MovSwProb<-NULL))
   MovSwProb<<-as.list(read.table(fileMovSwProb,sep=";",dec=","))
for(i in 1:length(MovSwProb)) {
  MovSwProb[[i]] <<- MovSwProb[[i]]/diff(c(0,probList$DistCat)^2) ## Normalizing with relative area
  MovSwProb[[i]][is.na(MovSwProb[[i]])] <<- 0
  MovSwProb[[i]] <<- MovSwProb[[i]]/sum(MovSwProb[[i]])
}
eval.parent(expression(MovWeProb<-NULL))
   MovWeProb<<-as.list(read.table(fileMovWeProb,sep=";",dec=","))
for(i in 1:length(MovWeProb)) {
  MovWeProb[[i]] <<- MovWeProb[[i]]/diff(c(0,probList$DistCat)^2) ## Normalizing with relative area
  MovWeProb[[i]][is.na(MovWeProb[[i]])] <<- 0
  MovWeProb[[i]] <<- MovWeProb[[i]]/sum(MovWeProb[[i]])
}
eval.parent(expression(SwMovAbProb<-NULL))
   SwMovAbProb<<-as.list(read.table(fileSwMovAbProb,sep=";",dec=",")) 
for(i in 1:length(SwMovAbProb)) {
  SwMovAbProb[[i]] <<- SwMovAbProb[[i]]/diff(c(0,probList$DistCat)^2) ## Normalizing with relative area
  SwMovAbProb[[i]][is.na(SwMovAbProb[[i]])] <<- 0
  SwMovAbProb[[i]] <<- SwMovAbProb[[i]]/sum(SwMovAbProb[[i]])
}
eval.parent(expression(MedRiskMovProb<-NULL))
   MedRiskMovProb<<-as.list(read.table(fileMedRiskMovProb,sep=";",dec=",")) 
for(i in 1:length(MedRiskMovProb)) {
  MedRiskMovProb[[i]] <<- MedRiskMovProb[[i]]/diff(c(0,probList$DistCat)^2) ## Normalizing with relative area
  MedRiskMovProb[[i]][is.na(MedRiskMovProb[[i]])] <<- 0
  MedRiskMovProb[[i]] <<- MedRiskMovProb[[i]]/sum(MedRiskMovProb[[i]])
}

eval.parent(expression(LowRiskMovProb<-NULL))
   LowRiskMovProb<<-as.list(read.table(fileLowRiskMovProb,sep=";",dec=",")) 
for(i in 1:length(LowRiskMovProb)) {
  LowRiskMovProb[[i]] <<- LowRiskMovProb[[i]]/diff(c(0,probList$DistCat)^2) ## Normalizing with relative area
  LowRiskMovProb[[i]][is.na(LowRiskMovProb[[i]])] <<- 0
  LowRiskMovProb[[i]] <<- LowRiskMovProb[[i]]/sum(LowRiskMovProb[[i]])
}

eval.parent(expression(toPlot<- NULL))#just to test
eval.parent(expression(toPlot2<- NULL))#just to test

## Further transformation of latent and subclinical duration frequencies
  herdtypes$latDurFreq<<-t(sapply(herdtypes$latDurFreq,function(x) if (is.character(x)) eval(parse(text=x)) else eval(x) ))
  herdtypes$SubDurFreq<<-t(sapply(herdtypes$SubDurFreq,function(x) if (is.character(x)) eval(parse(text=x)) else eval(x) ))
  herdtypes$CliDurFreq<<-t(sapply(herdtypes$CliDurFreq,function(x) if (is.character(x)) eval(parse(text=x)) else eval(x) ))

  ## Declaring some variables in the parent scope.(Accessible to all functions)
  eval.parent(expression(outbreakDetectedLast<-outbreakDetected<-FALSE))
  eval.parent(expression(depopQueue <- NULL))
  eval.parent(expression(beingDepoped <- NULL))
  eval.parent(expression(indexHerd<-NA))
  eval.parent(expression(iteration <- NULL))
  eval.parent(expression(infHerdNums <- NULL))
  eval.parent(expression(gTime<-0))
  eval.parent(expression(initHideMap <- hideMap))
   if(is.null(initHideMap))
    eval.parent(expression(hideMap <- T))

  eval.parent(expression(gMaxHerds <- length(aHerd$herdType)))
  aHerd$Sus <<- rep(NA,gMaxHerds)
  aHerd$K <<- rep(NA,gMaxHerds)
  aHerd$cullEligible <<- rep(T,gMaxHerds)

  if (!identical(cullTypes,"all"))
    aHerd$cullEligible <<- aHerd$herdType %in% cullTypes

  ## choose functions
  eval.parent(expression(
      RFtrans <-switch(RFstoch+1,
                    function(z,nn,pp){ nn*pp }, # if False
                    rbinom                      # if True
                    )
      ))
  eval.parent(expression(
      Tinfs<-switch(Tstoch+1,
                 function(contactHerds,dProbInfection,Sus) {
                   rbinom(length(contactHerds),1,dProbInfection)
                 },
                 function(contactHerds,dProbInfection,Sus) {
                   rbinom(length(contactHerds),aHerd$Sus[contactHerds],
                          1-(1-dProbInfection)^(1/Sus[contactHerds]))
                 })
           ))

  ## Set initial herd status to 1 for every herd unless values taken from file
  if (ignoreStatus) {
    eval.parent(expression(initStatus <- rep(1,gMaxHerds)))
    eval.parent(expression(initTimeInfected <- rep(Inf,gMaxHerds)))
  }
  else {
    eval.parent(expression(initStatus <- aHerd$status))
    eval.parent(expression(initTimeInfected <- aHerd$timeInfected))
  }

  ## Initiate distributions of animals moved and distribute them over the different herdsize categories
   HerdSizeCategories <- sort(unique(aHerd$herdSizeCat)) 
   HerdSizeCatDistAO  <- parse(text=c('round(rpert(n,1,1,2))','round(rpert(n,1,61,196))','round(rpert(n,1,85,415))','round(rpert(n,1,110,433))','round(rpert(n,1,81,371))'))
   HerdSizeCatDist    <- parse(text=c('round(rpert(n,1,1,2))','round(rpert(n,1,25,290))','round(rpert(n,1,304,800))','round(rpert(n,1,300,910))','round(rpert(n,1,290,800))'))
   
   aHerd$NumMovAnimal <<- rep(0,length(aHerd$herdType))
   for(i in HerdSizeCategories){
    Index <- aHerd$herdType==1&aHerd$herdSizeCat==HerdSizeCategories[i]
    aHerd$NumMovAnimal[Index] <<- rep(HerdSizeCatDistAO[i],sum(Index))
   }
   for(i in HerdSizeCategories){
    for(j in 2:10){
    Index <- aHerd$herdType==j&aHerd$herdSizeCat==HerdSizeCategories[i]
    aHerd$NumMovAnimal[Index] <<- rep(HerdSizeCatDist[i],sum(Index))
    }
   }

 eval.parent(expression(AllInfHerds <- NULL))
 AllInfHerds <<- matrix(numeric(0),ncol=7)
 NAMEInf <- paste(runID,"AllInfHerds.txt",sep="-")
 write.table(AllInfHerds,NAMEInf,sep=" ")

 eval.parent(expression(SumResOut  <- NULL))
 SumResOut  <<- matrix(numeric(0),ncol=32)
 NAME <- paste(runID,"ASF.txt",sep="-") 
 write.table(SumResOut,NAME,sep=" ")
######################################################################
### initiate the Matrix that includes surveyed herds (TH) Oct 2015 ###
######################################################################
if(Detailed){
 eval.parent(expression(ClSurvMatOut  <- NULL))
 eval.parent(expression(SerSurvMatOut <- NULL))
 eval.parent(expression(PCRSurvMatOut <- NULL))
 #eval.parent(expression(SurvZoneMatOut<- NULL))
 #eval.parent(expression(ProtZoneMatOut<- NULL))
 #eval.parent(expression(TraceDCMatOut <- NULL))
 #eval.parent(expression(TraceIDCMatOut<- NULL))
 eval.parent(expression(DepopMatOut   <- NULL))
 eval.parent(expression(PreEmpMatOut  <- NULL))
# eval.parent(expression(SurDeadMatOut <- NULL))


 ClSurvMatOut     <<- matrix(numeric(0),ncol=3)
 SerSurvMatOut    <<- matrix(numeric(0),ncol=3)
 PCRSurvMatOut    <<- matrix(numeric(0),ncol=3)
# SurvZoneMatOut   <<- matrix(numeric(0),ncol=3)
# ProtZoneMatOut   <<- matrix(numeric(0),ncol=3)
# TraceDCMatOut    <<- matrix(numeric(0),ncol=3)
# TraceIDCMatOut   <<- matrix(numeric(0),ncol=3)
 DepopMatOut      <<- matrix(numeric(0),ncol=3)
 PreEmpMatOut     <<- matrix(numeric(0),ncol=3)
# SurDeadMatOut    <<- matrix(numeric(0),ncol=3)

# NAMES    <- paste(runID,"SurvHerds.txt",sep="-")
# NAMEP    <- paste(runID,"ProtHerds.txt",sep="-")
 NAMEClSH <- paste(runID,"ClSurvayedHerds.txt",sep="-")
 NAMESerSH<- paste(runID,"SerSurvayedHerds.txt",sep="-")
 NAMEPCRSH<- paste(runID,"PCRSurvayedHerds.txt",sep="-")
# NAMETDC  <- paste(runID,"TDCHerds.txt",sep="-")
# NAMETIDC <- paste(runID,"TIDCHerds.txt",sep="-")
 NAMED    <- paste(runID,"DepopHerds.txt",sep="-")
 NAMEPE   <- paste(runID,"PreEmpHerds.txt",sep="-")
# NAMESD   <- paste(runID,"SurDeadHerds.txt",sep="-")

# write.table(SurvZoneMatOut,NAMES,col.names = F,row.names=F)
# write.table(ProtZoneMatOut,NAMEP,col.names = F,row.names=F)
 write.table(ClSurvMatOut,NAMEClSH,col.names = F,row.names=F)
 write.table(SerSurvMatOut,NAMESerSH,col.names = F,row.names=F)
 write.table(PCRSurvMatOut,NAMEPCRSH,col.names = F,row.names=F)
# write.table(TraceDCMatOut,NAMETDC,col.names = F,row.names=F)
# write.table(TraceIDCMatOut,NAMETIDC,col.names = F,row.names=F)
 write.table(DepopMatOut,NAMED,col.names = F,row.names=F)
 write.table(PreEmpMatOut,NAMEPE,col.names = F,row.names=F)
# write.table(SurDeadMatOut,NAMESD,col.names = F,row.names=F)
}
 }

##################################################################
### initializeASFvars
###
### This function is called at the beginning of each simulated epidemic
##################################################################
initializeASFvars <- function() {
  if (verbose) cat("iteration",iteration,"\n")

  ## If seed is negative the seed is set at the beginning of each iteration:
  if (!is.null(seed))
    if(seed<0)
      set.seed(iteration+abs(seed))
  
  ## Determine the day the first detection will happen
  gDaysUntilBaseline <<- 0

 
  ## output variable to calculate the number of herds queuing for surveillance per day. does not affect anything (TH)
  TotNumQueue <<- 0

   ## delay on meeting of OIE to give free status. 
     meetDelay <<- round(rpert(1,30,60,90))

#Additions for economic calculations (TH) 10-06-2014
  ## reduction on the price for products to EU
  reducedPriceEU <<- rpert(1,0.05,0.06,0.1)


  ## Capacity for surveillance. it is initiated here with 50 herds per day and further
  ## manipulated in the ASFEngine file.
    CapSurvay <<- 50

  ## Initialize waiting periods and transmission probabilities
  aHerd$taggedDur <<-rep(0,gMaxHerds)
  aHerd$DCFromNumEffCons <<- rep(0,gMaxHerds)
  aHerd$DCtoEffCons <<- rep(0,gMaxHerds)
  aHerd$LRICEffCons <<- rep(0,gMaxHerds)
  aHerd$relDC <<- rep(1,gMaxHerds)
  aHerd$relIMC <<- rep(1,gMaxHerds)
  aHerd$relILC <<- rep(1,gMaxHerds)
    ## New additions TH ###
  aHerd$timeToVisitTraceIDC  <<- rep(0,gMaxHerds)
  aHerd$timeToVisitTraceDC  <<- rep(0,gMaxHerds)
  aHerd$timeToPV1     <<- rep(0,gMaxHerds)
  aHerd$timeToPV2     <<- rep(0,gMaxHerds)
  aHerd$timeToSV1     <<- rep(0,gMaxHerds)
  aHerd$timeToSV2     <<- rep(0,gMaxHerds)
  aHerd$visitCount    <<- rep(0,gMaxHerds)
  aHerd$daysInSZones  <<- rep(0,gMaxHerds)
  aHerd$daysInPZones  <<- rep(0,gMaxHerds)
  aHerd$daysInZones   <<- rep(0,gMaxHerds)
  aHerd$traceIDC      <<- rep(0,gMaxHerds)
  aHerd$traceBDC      <<- rep(0,gMaxHerds)
  aHerd$allInZones    <<- rep(0,gMaxHerds)
  aHerd$DetTaggedDep  <<- rep(FALSE,gMaxHerds)
  aHerd$timeVisited   <<- rep(0,gMaxHerds)
  aHerd$timeSatInSZone<<- rep(0,gMaxHerds)
  aHerd$timeSatInPZone<<- rep(0,gMaxHerds)
  aHerd$infMode       <<- rep(0,gMaxHerds)
  aHerd$sampVisitPCR  <<- rep(0,gMaxHerds)
  aHerd$sampVisitSer  <<- rep(0,gMaxHerds)
  aHerd$Mortality     <<- rep(0,gMaxHerds)
  aHerd$NotInZoneYest <<- rep(TRUE,gMaxHerds)
  aHerd$NotInPZYest   <<- rep(TRUE,gMaxHerds)
  aHerd$VisitClInspect<<- rep(0,gMaxHerds)
  aHerd$immuneTime    <<- rep(0,gMaxHerds)
  aHerd$timeCulled    <<- rep(0,gMaxHerds)
  aHerd$ToSurvDead    <<- rep(0,gMaxHerds)
  aHerd$DeadSampTest  <<- rep(0,gMaxHerds)
  aHerd$SampledDead   <<- rep(0,gMaxHerds)
  aHerd$DiagSurvDead  <<- rep(FALSE,gMaxHerds)
  aHerd$SubDeadSamp   <<- rep(0,gMaxHerds)
  aHerd$NumSampDead   <<- ifelse(aHerd$herdSize>=numTestDead,numTestDead,aHerd$herdSize)
  aHerd$SusAgain      <<- rep(0,gMaxHerds)
  aHerd$Survived      <<- rep(0,gMaxHerds)


  ## This part works like list; it evaluates the text entries from the input
  ## table, whether those entries are numeric or R functions (like random
  ## variate selection), in order to set the waiting periods and transmission
  ## probabilities by herd type.
  for (i in herdtypes$herdTypeID) {
    herdIndex <- aHerd$herdType==i
    nType <- sum(herdIndex)
    typeIndex <- herdtypes$herdTypeID==i
    aHerd$K[herdIndex] <<- eval(herdtypes$K[typeIndex],list(n=nType)) ## Intra-herd interaction rate
   #############################################
   ### newly add for the Danish ASF project (TH)
   #############################################
   aHerd$RiskDC[herdIndex] <<-
      eval(herdtypes$RiskDC[typeIndex],list(n=nType))
   aHerd$RiskAb[herdIndex] <<-
      eval(herdtypes$RiskAb[typeIndex],list(n=nType))
    aHerd$LamMRC[herdIndex] <<-
      eval(herdtypes$LamMRC[typeIndex],list(n=nType))
   aHerd$RiskMRC[herdIndex] <<-
      eval(herdtypes$RiskMRC[typeIndex],list(n=nType))
   aHerd$LamLRC[herdIndex] <<-
      eval(herdtypes$LamLRC[typeIndex],list(n=nType))
   aHerd$RiskLRC[herdIndex] <<-
      eval(herdtypes$RiskLRC[typeIndex],list(n=nType))
   aHerd$RelSusceptibility[herdIndex] <<- 
      eval(herdtypes$RelSusceptibility[typeIndex],list(n=nType))
   aHerd$LocSpProfile[herdIndex] <<- 
      eval(herdtypes$LocSpProfile[typeIndex],list(n=nType))
   aHerd$RiskWB[herdIndex] <<- 
      eval(herdtypes$RiskWB[typeIndex],list(n=nType))
  }##EndOf for (i in ...)


  ## Making sure that both 'Sus' and 'herdSize' are in 'aHerd'
  ##LEC: This should be changed so that only one name is used!
  if("herdSize"%in%names(aHerd))
    aHerd$Sus<<-aHerd$herdSize
  else
    aHerd$herdSize<<-aHerd$Sus
 

  ## Making sure that the above distributed numbers make sense in the model:
  aHerd$Sus[aHerd$Sus<2] <<- 2
  aHerd$K[aHerd$K<0] <<- 0
  aHerd$taggedDur[aHerd$taggedDur<0] <<- 0

  aHerd$p <<- aHerd$K    # Reed-Frost probability
  aHerd$p[aHerd$p>1]<<-1
  
  for (i in 1:length(newInfMethods)) ## Initializing infection methods
    newInfMethods[[i]]$init()

  for (i in 1:length(controlMethods)) ## Initializing control measures
    controlMethods[[i]]$init()

#### Reset Values of other parameters
  gTime<<-0

  aHerd$timeToTaggedForDepop<<-rep(Inf,gMaxHerds)
  aHerd$Diagnosed<<-rep(F,gMaxHerds)
  aHerd$diagnosisTime <<- rep(Inf,gMaxHerds)
  aHerd$DiagSurv<<-rep(F,gMaxHerds)

  outbreakDetectedLast<<-outbreakDetected<<-FALSE
  
  depopQueue<<- matrix(numeric(0),ncol=2) ## rows: [ID, day when entered queue]
  beingDepoped<<- matrix(numeric(0),ncol=2) ## rows: [ID, days to depoped]
  traceMatrix <<- matrix(numeric(0),ncol=4)
  zoneQueue<<- matrix(numeric(0),ncol=4)

  aHerd$timeToRemovedSZ<<-rep(Inf,gMaxHerds)
  aHerd$timeToRemovedPZ<<-rep(Inf,gMaxHerds)
  aHerd$infSource<<-rep(NA,gMaxHerds)

#### set initial values for index herd(s)

  ## SelectindexHerd(s) for next simulation
  aHerd$timeInfected<<-initTimeInfected
  aHerd$status<<-initStatus
  if (ignoreStatus) {
    if(is.null(stepInFile)){
      indexHerd<<-indexHerdFunction(args=indexHerdSelect)
      aHerd$status[indexHerd] <<-2+indexDirect
      aHerd$timeInfected[indexHerd] <<-0
      
      if (indexDirect)
        aInfHerd$addInf(indexHerd,matrix(c(0,1,0),byrow=TRUE,ncol=3,nrow=length(indexHerd)),1)
      else
        aInfHerd$addInf(indexHerd,matrix(c(1,0,0),byrow=TRUE,ncol=3,nrow=length(indexHerd)),0)
    }
    else{ ## Using stepInFile
      tmp<-read.table(stepInFile,sep=",",header=T)
      names(tmp)[1]<-"ID" #First column must be unique ID
      aHerdIndex<-match(tmp$ID,aHerd$ID)
      assign("stepIn",cbind(aHerdIndex,tmp),envir=parent.env(environment()))
      if(any(is.na(aHerdIndex)))
        stop("stepInFile: ",stepInFile," contains IDs not in herdfile. \n Missing IDs are:",
             paste(stepIn$ID[is.na(aHerdIndex)],collapse=", "))
      indexHerd<-aHerdIndex[stepIn$infectionDay == min(stepIn$infectionDay)]
      gDaysUntilBaseline<<-min(stepIn$diagnosedDay)
      if("depopedDay"%in%names(stepIn)){
        traced<-stepIn$diagnosedDay>=stepIn$depopedDay
        ## "diagnosing" the herds time to diagnosis before depop:
        stepIn$diagnosedDay<<-stepIn$depopedDay-aHerd$timeToDiagnosis[stepIn$aHerdIndex]
      }
      ##DK for (day in min(stepIn$infectionDay):(gDaysUntilBaseline-1)){
      for (day in min(stepIn$infectionDay):max(stepIn$infectionDay)){
        gTime<<-day
        newInf<-which(stepIn$infectionDay==day)
        if(length(newInf)>0){## Add these
          for (j in 1:length(newInf)){ ## One at a time
            if(stepIn$infType[ newInf[j] ]==1){ ## Direct => subClin
              aInfHerd$addInf(aHerdIndex[ newInf[j] ],cbind(0,1,0),1)
              aHerd$status[aHerdIndex[newInf[j]]] <<- 3      # subclin infection status
              ##DK aInfHerd$addInf(aHerdIndex[ newInf[j] ],cbind(1,0,0),1)
              ##DK Ups:  the scenario dictates that the moved animal is latent
            }
            else{
              aInfHerd$addInf(aHerdIndex[ newInf[j] ],cbind(1,0,0),0)
              aHerd$status[aHerdIndex[newInf[j]]] <<- 2      # latent infection status
            }
          }
        }## Added newInfected herds
        aInfHerd$simDay()
      }
      aHerd$timeInfected[stepIn$aHerdIndex]<<-stepIn$infectionDay
      aHerd$infSource[stepIn$aHerdIndex]<<-stepIn$infectionSource
      ##infHerdNums <<- (aInfHerd$getIDs())[aInfHerd$getstatus()%in%(3:4)]
### Time line!!!
### Set Time Diagnosed as to make sure depop happens as it should!!!      
    }
  }
  else{ ## What to do for infected herds given their time of diagnosis
    

  }
  
  vaccToday <<- 0
  Lock <<- FALSE

  if(is.null(initHideMap))
    hideMap <<- TRUE

#### chronicle index herds:
  if(sum(aHerd$status==2)>0)
    chronicle$addEntry(itn=iteration, state=2, newInfection=1, time=1,
                       changedHerdsList=aHerd$ID[aHerd$status==2],
                       HerdSize=aHerd$herdSize[aHerd$status==2],
                       HerdTypes=aHerd$herdType[aHerd$status==2])
  if(sum(aHerd$status==3)>0)
    chronicle$addEntry(itn=iteration, state=3, newInfection=1, time=1,
                       changedHerdsList=aHerd$ID[aHerd$status==3],
                       HerdSize=aHerd$herdSize[aHerd$status==3],
                       HerdTypes=aHerd$herdType[aHerd$status==3])
  if(sum(aHerd$status==4)>0)
    chronicle$addEntry(itn=iteration, state=4, newInfection=1, time=1,
                       changedHerdsList=aHerd$ID[aHerd$status==4],
                       HerdSize=aHerd$herdSize[aHerd$status==4],
                       HerdTypes=aHerd$herdType[aHerd$status==4])
}

