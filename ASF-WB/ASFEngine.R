#### ASfEngine library version 0.15.1
#### Includes the following functions:
####
#### ASFEngine (TH)
####

ASFEngine <- function() {
  if(verbose) cat("Started engine\n")

  for(iteration in 1:n) { ## Looping trough n simulated epidemics
    iteration <<- iteration
    initializeASFvars()             # function in INITIALIZE library

     while (any(aHerd$status!=1 & aHerd$status!=5 & aHerd$status!=6 & aHerd$status!=7) & gTime < maxTime) {

      gTime<<-gTime+1
      outbreakDetectedLast<<-outbreakDetected
      infHerdNumsLast <<- infHerdNums ## Vector of infectious IDs 

   ### Increase the Surveillance capacity 
       if(gTime>=(gDaysUntilBaseline+7))  {CapSurvay <<- 100}
       if(gTime>=(gDaysUntilBaseline+14)) {CapSurvay <<- 200}
      

      if(verbose) cat("Time:",gTime,"\n")
      updateHerds( )                # function in PROCESS_HERDS library
      limitMovements( )             # function in LIMIT_MOVEMENTS library
      
    ## Update relative effect of measures on DC, IMC and ILC
    aHerd$relDC<<-aHerd$relIMC<<-aHerd$relILC<<-rep(1,gMaxHerds) ## All start with one
      if (outbreakDetected){
         for (i in 1:length(controlMethods)){ ## The controls may lower some.
          controlMethods[[i]]$day()
        }
      }
   
      if (outbreakDetected) {
        for (i in 1:length(InterMethods)){
          InterMethods[[i]]$day() ## Calling the intervention functions daily
        }      
      }  
      ## Finding the indexes of infectious herds
      tmpTagged<-which( aHerd$status==6 & aHerd$timeInfected<Inf & aHerd$SusAgain==0)
      infHerdNums <<- unique(c(tmpTagged,
                             (aInfHerd$getIDs())[aInfHerd$getstatus()%in%(3:4)]))
      if (verbose) print(c("aInfHerd$getIDs: ",aInfHerd$getIDs()))
   
      ## Making newInfections if there still are some infected herds.
      ## newInfFunctions is given as an input to ASFoptions() and translated
      ## into newInfMethods
      if ( length(infHerdNums)>0 & length(aInfHerd$getIDs())>0 ){
          for (i in 1:length(newInfMethods)){

            newInfMethods[[i]]$day() ## Calling the newInfFunctions daily
 
          }##EndOf for
        }##EndOf if
      summaryFunction("day")
print(c(iteration,gTime))
#print(table(aHerd$status))

    }##EndOf while
  
 summaryFunction("iter")
    ## Reseting the distance function and the infected herd function
    Dist$wipe()
    aInfHerd$wipe()
    for (i in 1:length(newInfMethods))
      newInfMethods[[i]]$cleaniteration() ## Calling the newInfFunctions to clean up
    for (i in 1:length(InterMethods))
      InterMethods[[i]]$cleaniteration() ## Calling the intervention functions to clean up
    gc(verbose=F,reset=T) #Garbage collecting to free memory.
  }##EndOf for(iteration in 1:n)
}##EnfOf function()
