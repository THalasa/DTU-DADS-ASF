source("../DTU-DADS-ASF/sourceASF.R")
sourceASF()
indexHerdFile<-read.table("../ASFinputs/IndexHerds.csv",sep=";")
optC<-ASFoptions(maxTime=16,n=1,indexHerdSelect=list(ID=indexHerdFile[[2]]),runID="TestSows1",Detailed=F,PCRTesting=c(1,5,6),seed=-10)
a<-ASF(optC)


source("../DTU-DADS-ASF/sourceASF.R")
sourceASF()
indexHerdFile<-read.table("../ASFinputs/IndexHerds.csv",sep=";")
optC<-ASFoptions(maxTime=365,n=5,indexHerdSelect=list(ID=indexHerdFile[[2]]),Detailed=T,runID="TestSows3",
                 interventionFunctions=c("SurvDead(Zone=2)"),seed=-10)
a<-ASF(optC)

source("../DTU-DADS-ASF/sourceASF.R")
sourceASF()
indexHerdFile<-read.table("../ASFinputs/IndexHerds.csv",sep=";")
optC<-ASFoptions(maxTime=365,n=5,indexHerdSelect=list(ID=indexHerdFile[[2]]),Detailed=T,runID="TestSows2",
                 interventionFunctions=c("CullRing(size=1,startCulling=1)"),seed=-10)
a<-ASF(optC)

source("../DTU-DADS-ASF/sourceASF.R")
sourceASF()
indexHerdFile<-read.table("../ASFinputs/TestIndHerds.csv",sep=";")
optC<-ASFoptions(maxTime=365,n=10,indexHerdSelect=indexHerdFile,indexHerdFunction="selectIndHerdMore",Detailed=T,runID="Test",seed=-10)
a<-ASF(optC)