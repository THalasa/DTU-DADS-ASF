#### ASF model interface version 0.15.1
#### Includes the following functions:
####
#### FMD (TH)
####

#### List of optional function input variables:
## args         list of options
## test         if TRUE FMD() exists without doing calculations.
##              This feature is only included for use .Rd files
##              At some point a small default dataset should be included in
##              the package.

ASF <- function(args=ASFoptions() ,test=FALSE ) {
  if (test){
    print("ASF: Just a test - didn't calculate anything!   ;-)")
    return(0)
  }

  if (class(args)=="list"){
    argnames<-names(args)

    for (i in 1:length(args)){
      assign(argnames[i],args[[i]])
    }
  }
  else
    stop("args must be a list containing all required parameter values!")
  if(verbose) cat("Assigned all options\n")
 
  ## Parsing textual parameters:
  if (!is.expression(delayTimes))
    delayTimes<-parse(text=delayTimes)
  if (!is.expression(newInfFunctions))
    newInfFunctions<-parse(text=newInfFunctions)
  if (!is.expression(controlFunctions))
    controlFunctions<-parse(text=controlFunctions)
  if (!is.expression(interventionFunctions))
    interventionFunctions<-parse(text=interventionFunctions)
  if(verbose) cat("Parsed textual parameters\n")

  ## Prespecified seeds can be used
  if (!is.null(seed))
    if(seed<0)
      set.seed(1) ## This is further handled when initializing each iteration
    else
      set.seed(seed)

  ## In order not to mess around in the global environment all functions
  ## that need access to internal variables need to be re-scoped.
  ## Re-scoping the internal functions:
  environment(ASFEngine)<-environment()
  environment(createASFvars)<-environment()
  environment(initializeASFvars)<-environment()
  environment(limitMovements)<-environment()
  environment(updateHerds)<-environment()
  environment(constructChronicle)<-environment()
  environment(constructDist)<-environment()
  environment(ConstpMat)<-environment()
  environment(constructAInfHerd)<-environment()
  environment(calDist)<-environment()

  ## Re-scoping graphics functions:
  environment(fmdplot)<-environment()
  environment(cbands)<-environment()
  environment(points1)<-environment()
  environment(points2)<-environment()
  environment(updateCharts)<-environment()
  environment(epimap)<-environment()
  environment(filled3)<-environment()

  ## Re-scoping user supplied functions:
  environment(indexHerdFunction)<-environment()
  environment(summaryFunction)<-environment()
 

  newInfList<-list()
  newInfFun<-list()
  for (i in 1:length(newInfFunctions)){
    newInfList[[i]]<-as.list(newInfFunctions[[i]]) # Splits into name and arg
    newInfFun[[i]]<-eval(newInfList[[i]][[1]])  # Get the function
    environment(newInfFun[[i]])<-environment()  # Rescope the function
  }


  controlList<-list()
  controlFun<-list()
  for (i in 1:length(controlFunctions)){
    controlList[[i]]<-as.list(controlFunctions[[i]]) # Splits into name and arg
    controlFun[[i]]<-eval(controlList[[i]][[1]])  # Get the function
    environment(controlFun[[i]])<-environment()  # Rescope the function
  }


  InterList<-list()
  InterFun<-list()
  for (i in 1:length(interventionFunctions)){
    InterList[[i]]<-as.list(interventionFunctions[[i]]) # Splits into name and arg
    InterFun[[i]]<-eval(InterList[[i]][[1]])  # Get the function
    environment(InterFun[[i]])<-environment()  # Rescope the function
  }

  if(verbose) cat("Rescoped functions internally\n")
  
  ## Constructing the chronicle, that logs all changes of status
  ## including newInfections and the mode of infection.
  chronicle<-constructChronicle(herdfile=infofile, typesfile=typesfile)

  createASFvars()          # function in INITIALIZE library
  summaryFunction("init")

  ## Initializing distance related objects
  aInfHerd<-constructAInfHerd()
  if ("long"%in%names(aHerd)) ## latlong or UTM ?
    Dist<-constructDist(aHerd$lat,aHerd$long,norm=2)
  else                        
    Dist<-constructDist(aHerd$north,aHerd$east,norm=3)

  ## Creating and initializing the newInfFunctions
  ## each of which returns a list of functions: day and cleaniteration
  newInfMethods<-list()
  for (i in 1:length(newInfFunctions)){
    newInfMethods[[i]]<-do.call(newInfFun[[i]],newInfList[[i]][-1])
  }
  ## Same thing for the control measures:
  controlMethods<-list()
  for (i in 1:length(controlFunctions)){
    controlMethods[[i]]<-do.call(controlFun[[i]],controlList[[i]][-1])
  }
  ## Same thing for the optional interventions:
  InterMethods<-list()
  for (i in 1:length(interventionFunctions)){
    InterMethods[[i]]<-do.call(InterFun[[i]],InterList[[i]][-1])
  }
  
  ## Start and run the model now
  ASFEngine()

  ## Wrapping up
  chronicle$wrapUp(savefile=chroniclefile)

  ## Making final summary including what is returned by ASF()
  summaryFunction("final")  ## This should be the last line !!!
}
