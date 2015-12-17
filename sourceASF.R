##
## sourceASF V. 0.15.1
##
## This is a utility function for sourcing the R code for the ASF model.
## This file/function should NOT be included in the package.
sourceASF<-function(verbose=FALSE){
  if(verbose) cat("Reading in Dist library...",fill=T)
  source("Dist.R")
  #if(verbose) cat("Reading in ASFEngine library...",fill=T)
  source("ASFEngine.R")
  if(verbose) cat("Reading in initialize library...",fill=T)
  source("initialize.R")
  if(verbose) cat("Reading in limitMovements library...",fill=T)
  source("limitMovements.R")
  if(verbose) cat("Reading in process herds library...",fill=T)
  source("process_herds.R")
  if(verbose) cat("Reading in opt_intervention library...",fill=T)
  source("opt_intervention.R")
  if(verbose) cat("Reading in graphics library...",fill=T)
  source("graphics.R")
  if(verbose) cat("Reading in chronicle library...",fill=T)
  source("chronicle.R")
  if(verbose) cat("Reading in misc library...",fill=T)
  source("misc.R")
  if(verbose) cat("Reading in ASFoptions library...",fill=T)
  source("ASFoptions.R")
  if(verbose) cat("Reading in ASF library...",fill=T)
  source("ASF.R")
  if(verbose) cat("Reading in control measures library...",fill=T)
  source("controlMeasures.R")
  if(verbose) cat("Reading in ConstpMat library...",fill=T)
  source("ConstpMat.R")
  if(verbose) cat("Reading in summaries library...",fill=T)
  source("summaries.R")
  cat("Library sourcing complete.\n")
}

Tidy<-function(){
	files<-dir(pattern=".tmp")
	files<-c(files,dir(pattern=".bix"))
	
	if(length(files)>0){
		for(i in 1:length(files)) cat(files[i],"\n")
		makesure<-readline(prompt="Are you sure you want to delete all these? (y/n) ")
		if(makesure=="y") {
			unlink(files)
			cat("All done\n")
		}
		else
		  cat("No files deleted\n")
	}
	else
	  cat("No .tmp or .bix files found\n")
}

