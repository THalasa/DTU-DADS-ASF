source("sourceASF.R")
sourceASF()
indexHerdFile<-read.table("IndexHerdsWB.csv",sep=";")
optC<-ASFoptions(maxTime=365,n=1160,indexHerdSelect=list(ID=indexHerdFile[[1]]),
newInfFunctions=c( # Vector of functions used to make new infections (including parameters).
                   "DIRinf3('LamAll',MovSwProb,'pMatAll','RiskDC',MovMatAll,restMovedSize=35,label=1)", 
                   "DIRinf3('LambdaWeaners',MovWeProb,'pMatWea','RiskDC',MovMatWean,restMovedSize=10,label=1)", 
                   "INDflex('LamAb',SwMovAbProb,'relDC','pMatMovAb','RiskAb',probMatrix=MovAb,Reduction=0.5,Abattoir=TRUE,label=2)",
                   "INDflex('LamMRC',MedRiskMovProb,'relIMC','pMatMRC','RiskMRC',Reduction=1,label=3)",
                   "INDflex('LamLRC',LowRiskMovProb,'relILC','pMatLRC','RiskLRC',Reduction=1,label=4)",
                   "LASinf(localsize=2,label=5)"),
runID="Basic",seed=-10)
a<-ASF(optC)

