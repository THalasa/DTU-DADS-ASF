#### CHRONICLE library version 0.05
#### Includes the following functions:
####
#### constructChronicle
####

## chronicle is a ASF function that sets up chronicle, the log of herd changes
## during a run, and functions associated with it.  
## Most of these functions are called outside of the main loops, so speed is
## not important; the exception is $addentry 
## Variables of chronicle are documented as they are declared.

constructChronicle<-function(herdfile=warning("No herdfile specified when chronicle Initialized"), typesfile=warning("No typesfile specified when Chronicle Initialized"), dumpTrigger=1e5) {

 ######################## INITIALIZE AND ADD DATA:

  number<-0    ## number of dumps-to-file that have taken place
  prevEntry<-0 ## keeps track of what line the list is at
  
  ## Choose an ID number for the run; this implements a random number, you
  ## could use a date/time; or could use a globally scoped runid from elsewhere
  if(exists("runID"))
    runID<<-get("runID")
  if(is.null(runID)||(runID==""))
    runID<<-round(runif(1,10000,99999))
  
  ## Check to make sure the id number is not being used by the system.
  ## while(runID %in% read.table(System.RunID)
  ## runID<<-round(runif(1,10000,99999)) 
  cat("\nFYI, chronicle has assigned ID",runID,"to this run\n")

  ## chronicle is set up to save locally and zero itself on a regular basis for
  ## several reasons: 
  ## 1. If the node or run dies during simulation, you have a backup of what
  ## happened 
  ## 2. Data structures beyond a certain size just work faster -- saving,
  ## loading, analyzing 
  ## 3. Memory will always be limited.
  ## So chronicle_dump runs when the entries are greater than Trigger, after an
  ## iteration begins (ie if finishes it's current iteration and then dumps). 
  ## The following sets up the trigger (which is also used for initialization)
  ## and the boolean that controls dumps 
  trigger   <- dumpTrigger
  dumpSaves <- TRUE
  if(!dumpTrigger) {
    trigger   <- 1000
    dumpSaves <- FALSE
  }

  ## chronicle keeps track of herd by an ID number specified in the file with
  ## info about herds, so to be analyzed we'll have to go back 
  ## to that same file.  The same is generally true of a types file: we may
  ## eventually have typesfiles tied to herdsfiles (e.g. if Texas 
  ## has different herd-types than California, the herdfile for each state
  ## would reference a different typesfile). 
  ##
  ## Since we're initially using only one herdfile, this will be a single
  ## value, but it will eventually be a list of filenames to open up  
  ## with herd info in them.
  
  ## These two values are going to be lists of vectors
  stateHerds  <- list()
  sourceHerds <- list()
  
  ##Creating the remaining parameters
  iteration<-NULL
  day<-NULL
  newState<-NULL
  newInfMode<-NULL
  dumpedFiles<-NULL
  herdsizes<-NULL
  herdTyps<-NULL
  

  ## Variables used in relation to traceBack
  ## Each line in matrices contains: Time, ChangedHerds, and mode of Infection
  traceList<-list() # Should be cleared by the end of each iteration

  
############### SECOND, ADD INTERNAL FUNCTIONS:
  
  zero<-function(){
    ## This function is called both during initialization and by dump() 
    
    ## R seems to take as much time assigning a single value into a vector as
    ## it takes to assign a vector of multiple values. 
    ## So we can save time by making initial assignments of a reasonable
    ## length.
    
    ## Keep track of how many times chronicle has been Dumped    
    number <<- number+1 
    prevEntry <<- 0 ## Reset the number of lines in the chronicle

    ## Empty out the herd lists but keep the length
    stateHerds <<- list()
    stateHerds[[trigger]] <<- array(dim=1)
    sourceHerds <<- stateHerds
    herdsizes <<- stateHerds
    herdTyps <<- stateHerds
    iteration <<- array(dim=(trigger)) ## Tracks iteration
    day       <<- array(dim=(trigger)) ## Tracks day
    
    ## Tracks the new status of the disease in the stateHerds
    newState    <<- array(dim=(trigger))
    
    ## Tracks how the herds got infected: 
    newInfMode  <<- array(dim=(trigger))  
    ## Legend to newInfMode / mode of infection:
    ## 0 - not a newInfection. 
    ## 1 - index case. 
    ## 2 - INDinf.
    ## 3 - DIRinf.
    ## 4 - LASinf.
    ## 5 - airborne (5 is not implemented)
  }##EndOf zero
  
  ## Function called from addEntry() when chronicle gets bigger than
  ##  Trigger   
  dump<-function(){
    savefile<-paste(number,runID,"ASF.tmp",sep="-")
    save(stateHerds, sourceHerds, herdsizes, herdTyps,
         iteration, day, newState, newInfMode, file=savefile)
    zero()
    dumpedFiles<<-c(savefile, dumpedFiles)
    if (verbose) cat(paste("Dumped: ",savefile,"\n",sep=""))
  }##EndOf dump

  
######################## END OF INTERNAL FUNCTIONS

  ## Finally, call the function that sets the list values to zero
  zero()

  ## Error messages for external methods: addEntry and wrapUp:
  noItr<-"Error in chronicle$addEntry, iteration not specified"
  noState<-"Error in chronicle$addEntry, state not specified"
  noHerdsList<-"Error in chronicle$addEntry, herdlist not specified"
  noHerdSize<-"Error in chronicle$addEntry, herdSize not specified"

  

  ## List containing the functions to be called from the outside.
  list(
       ## Function to add a new line to the chronicle
       ## This function is called routinely during engine operation
       ## - so keep it fast!   
       addEntry=function(itn=stop(noItr), state=stop(noState), newInfection=0,
         time=gTime, changedHerdsList=stop(noHerdsList), SourceHerds=NA, 
         HerdSize=stop(noHerdSize),HerdTypes=stop(noHerdTypes)){
         if(is.na(changedHerdsList[1]))
           stop("CHL problem in chronicle")

         ## dump chronicle if the number of lines is greater than trigger
         if(dumpSaves && prevEntry>=trigger)
           dump()

         thisentry<<-prevEntry+1    
         iteration[thisentry] <<- itn
         day[thisentry] <<- time ## default=gTime: scoped globally
         newState[thisentry] <<- state
         newInfMode[thisentry] <<- newInfection
         stateHerds[[thisentry]] <<- changedHerdsList
         sourceHerds[[thisentry]] <<- SourceHerds
         herdsizes [[thisentry]] <<- HerdSize
         herdTyps [[thisentry]] <<- HerdTypes
         prevEntry<<-thisentry

         if(traceBack){ ## Record if doing trace backs
           if (newInfection>0){ ## but only if a newInfection
             if (newInfection==1){
               if (verbose){ print("TraceList:"); print(traceList)}
               traceList<<-list() ## reset when index case
             }
             charSource<-as.character(SourceHerds)
             for ( i in 1:length(SourceHerds)){
               traceList[[ charSource[i] ]]<<-
                 rbind(traceList[[ charSource[i] ]],
                       matrix(c(time,changedHerdsList[i],newInfection),ncol=3))
             }##EndOf for ( i in ...)
           } ##EndOf if(newInfection>0)
         }##EndOf if (traceBack)

       }##EndOf addEntry
       ,
       wrapUp=function(savefile=FALSE){
         ## Making local copies of variables to wrap up to avoid load()
         ## overwriting the values. This could also have been done by
         ## wrapping in another name before dumping. 
         WUiteration    <- iteration   
         WUday          <- day         
         WUnewState     <- newState    
         WUnewInfMode <- newInfMode
         WUstateHerds   <- stateHerds  
         WUsourceHerds  <- sourceHerds 
         WUherdsizes <- herdsizes 
         WUherdTyps <- herdTyps 

         ## first check and see if this is necessary
         if(verbose) cat("Gathering up chronicle data\n")
         if(length(dumpedFiles)>0){
           if(verbose) print(dumpedFiles)
           for(filei in dumpedFiles){
             load(filei)
             WUiteration    <- c(iteration,WUiteration)
             WUday          <- c(day,WUday)
             WUnewState     <- c(newState,WUnewState)
             WUnewInfMode <- c(newInfMode,WUnewInfMode)
             WUstateHerds   <- c(stateHerds,WUstateHerds) 
             WUsourceHerds  <- c(sourceHerds,WUsourceHerds)
             WUherdsizes <- c(herdsizes,WUherdsizes)
             WUherdTyps <- c(herdTyps,WUherdTyps)
            }
         }
         if(verbose) cat("chronicle data gathered, writing into idx file")

         ## Create the bix file; And removing NAs
         daysNotNA<-!is.na(WUday)
         idx<-list(iteration    = WUiteration[daysNotNA],
                   day          = WUday[daysNotNA],         
                   newState     = WUnewState[daysNotNA],    
                   newInfMode   = WUnewInfMode[daysNotNA],
                   herdsizes    = WUherdsizes[daysNotNA],
                   herdTyps     = WUherdTyps[daysNotNA],
                   stateHerds   = WUstateHerds[daysNotNA],  
                   sourceHerds  = WUsourceHerds[daysNotNA] )
         idx$runID     <- runID
         idx$herdfile  <- herdfile
         idx$typesfile <- typesfile
         
         if(!(savefile))
           savefile<-paste(runID,"ASF.bix",sep="-")
         save(idx,file=savefile);
         cat("\nSaved idx in ",savefile,"\n")
         if(length(dumpedFiles)>0){
           file.remove(dumpedFiles)
           if (verbose) cat("Removed temporary files\n")
         }
       }##EndOf wrapUp
       )##EndOf list being returned
}##EndOf constructChronicle function
