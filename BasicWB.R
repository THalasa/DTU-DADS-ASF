source("sourceASF.R")
sourceASF()
indexHerdFile<-read.table("IndexHerdsWB.csv",sep=";")
optC<-ASFoptions(maxTime=365,n=1160,indexHerdSelect=list(ID=indexHerdFile[[1]]),runID="BasicWB",seed=-10)
a<-ASF(optC)

