

rpert <- function(n,a,l,b) {
  mu <- (a+4*l+b)/6
  if (mu==l) v <- w <- 3 else {
    v <- (mu-a)*(2*l-a-b)/(l-mu)/(b-a)
    w <- v*(b-mu)/(mu-a)
    }
  a+(b-a)*rbeta(n,v,w)
}

set.seed(10)
#Outcome <- matrix(numeric(0),ncol=4)
 MatTTC <- matrix(numeric(0),ncol=7)
#write.table(MatTTC,"TimeToClearB6b.txt", sep=" ",col.names = F,row.names=F)
InfProbSubPar <- 0.5 #c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
DeadImpactPar <- 0.5 #c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

for(b in InfProbSubPar){
 for(j in DeadImpactPar){ 

InfProbSub  <- b
DeadImpact  <- j
PerDeadAnim <- 0.95  # Percentage of the animals that would die
#SelectHerds <- 10000 # Number of herds to model
PersTime    <- 5    # Number of days following death, the left over of a dead animal can be infectious.
Gamma       <- 1
# InfPropSub is the Relative risk of infectiousness of subclinically infected animals
# DeadImpact  is the Impact of materials of dead animals on infection process
# Gamma is A parameter to address uncertainty regarding virus survival time in dead materials

herdtypes <- as.list(read.table("typesfile.csv",header=T,sep=";",as.is=TRUE)) 
for (i in 3:length(herdtypes))
   herdtypes[[i]]<-parse(text=herdtypes[[i]])

  nColHerds<-20
  herds<-matrix(numeric(0),ncol=nColHerds)
  delHerds<-matrix(numeric(0),ncol=nColHerds) ##Keeping deleted herds
  herdsCol<-c("Sus1","latent2","SubC3","Clinic4","Immune5","Total6","status7",
              "ID8","p9","NA.latDur10","NA.SubDur11","TClic12","TDiag13",
              "Diagnosed14","InfByDC15","TInfected16","herdType17",
              "Vaccinated18","TLastAnCli19","Tested")
  dimnames(herds)[[2]]<-herdsCol
  
  ## VarDef latent,subclinical  Used to model transitions between stages
  ## VarDes Matrix with a row for each farm and columns for
  ## VarDes   each day latent or subclinical or clinical
  herdtypes$latDurFreq<-t(sapply(herdtypes$latDurFreq,function(x) if (is.character(x)) eval(parse(text=x)) else eval(x) ))
  herdtypes$SubDurFreq<-t(sapply(herdtypes$SubDurFreq,function(x) if (is.character(x)) eval(parse(text=x)) else eval(x) ))
  herdtypes$CliDurFreq<-t(sapply(herdtypes$CliDurFreq,function(x) if (is.character(x)) eval(parse(text=x)) else eval(x) ))

  latent <- matrix(numeric(0),ncol=ncol(herdtypes$latDurFreq))
  subclinical <- matrix(numeric(0),ncol=ncol(herdtypes$SubDurFreq))
  clinical <- matrix(numeric(0),ncol=ncol(herdtypes$CliDurFreq))

  "rpoly2"<-function(np){
    tabulate(sample(1:(length(np)-1),size=np[1],replace=TRUE,
                    prob=np[-1]),nbins=(length(np)-1))
  }
      gTime <- 1
      Days  <- 365
       Sus <-  read.table("HerdsWM.csv",sep=";")#round(runif(SelectHerds,min=2,max=10000))#k$herdSize#round(rexp(SelectHerds,rate=0.0018))#round(runif(SelectHerds,min=6,max=300))### round(runif(SelectHerds,min=1201,max=5000))## ###Number of animals within each herd. reflect the distribution of herd sizes from my data
       Sus <- Sus[,1]
       SelectHerds<-length(Sus)
       IDs <- 1:SelectHerds                                             # and corrected to an average herd size of about 620 animals, but it varies a bit.
      
       K <- 0.3 #rpert(SelectHerds,0.1,0.36,0.5)        # number of infectious contacts per day per infectious animal (median value)

       p <- K # probability of infectious contact within the specific-herd per day per infectious animal
       

      tmp<-cbind(Sus,1,0,0,0,Sus+1,2,IDs,
                 p,0,0,
                 Inf,Inf,0,0,gTime,1,
                 0,0,0)
      herds<-rbind(herds,tmp)

## Distributing the animals on No. days being latent or subclincal or clinical
           latent<-rbind(latent,t(apply(cbind(herds[,2],herdtypes$latDurFreq[herds[,17],,drop=FALSE]),1,rpoly2)))
           subclinical<-rbind(subclinical,t(apply(cbind(herds[,3],herdtypes$SubDurFreq[herds[,17],,drop=FALSE]),1,rpoly2)))
           clinical<-rbind(clinical,t(apply(cbind(herds[,4],herdtypes$CliDurFreq[herds[,17],,drop=FALSE]),1,rpoly2)))

Clinicals   <- matrix(numeric(0),nrow=SelectHerds,ncol=Days)
Infectness  <- matrix(numeric(0),nrow=SelectHerds,ncol=Days)
DiedWithin  <- matrix(numeric(0),nrow=SelectHerds,ncol=Days)
ImmHerds    <- matrix(numeric(0),ncol=3)
DeadMat     <- matrix(numeric(0),nrow=SelectHerds)
DEADYest    <- 0

      while(gTime<Days){
 
          gTime <- gTime+1
           ## Sus -> latent
           rS2L<-rbinom(nrow(herds),herds[,1],
                 1-exp(- herds[,9]*( (herds[,3]*InfProbSub)+(herds[,4]) + 
                 ((rowSums(t(exp((-dim(DeadMat)[2]+1):0*Gamma)*t(DeadMat))))*PerDeadAnim*DeadImpact) ) / (herds[,6]-(herds[,5]*PerDeadAnim))))

           ## Updating totals
           herds[,1:5]<-herds[,1:5]+cbind(-rS2L,rS2L-latent[,1],
                                           latent[,1]-subclinical[,1],
                                           subclinical[,1]-clinical[,1],
                                            clinical[,1])

Infectness[,gTime] <- rS2L
           ## Distributing new latent infected
           newLat<-apply(cbind(rS2L,herdtypes$latDurFreq[herds[,17],,drop=FALSE]),1,rpoly2)
           ## Distributing new subclinical infected
           newSubC<-apply(cbind(latent[,1],herdtypes$SubDurFreq[herds[,17],,drop=FALSE]),1,rpoly2)
           ## Distributing new clinical infected
           newCli<-apply(cbind(subclinical[,1],herdtypes$CliDurFreq[herds[,17],,drop=FALSE]),1,rpoly2)
           ## Update intra latent and subclinical and clinical waiting states
           latent<-cbind(latent[,-1,drop=FALSE],0)+t(newLat)
           subclinical<-cbind(subclinical[,-1,drop=FALSE],0)+t(newSubC)          
           DeadMat <- cbind(DeadMat,clinical[,1])
           if(gTime>(PersTime+1))
             DeadMat  <- DeadMat[,-1,drop=FALSE]

DiedWithin[,gTime]  <- clinical[,1]          
DEADYest <- clinical[,1] 
Clinicals[,gTime]  <- herds[,4]

           clinical<-cbind(clinical[,-1,drop=FALSE],0)+t(newCli)
           ## Checking for herds that changed status
           ## status 2 -> 3 :
           tmpIndex<-herds[,7]==2 & herds[,3]>=1
           if(sum(tmpIndex)>0){
              herds[tmpIndex,7] <- 3
           }
           ## status 3 -> 4 :
           tmpIndex<-herds[,7]==3 & herds[,4]>=1
           if(sum(tmpIndex)>0){
             herds[tmpIndex,7] <- 4
             herds[tmpIndex,12] <- gTime
           }
           ## status 4 -> 7 :
           tmpIndex<-(herds[,7]==4 & herds[,5]==herds[,6]) | (herds[,7]==4 & herds[,2]==0&herds[,3]==0&herds[,4]==0
                      &herds[,5]>0&rowSums(DeadMat)==0) #&DEADYest==0
           if(sum(tmpIndex)>0){
             ImmHerds <- rbind(ImmHerds,cbind(herds[tmpIndex,8],gTime,herds[tmpIndex,6]))
             herds[tmpIndex,7] <- 7
           }#End if
        }#End while

InfectnessD <- Infectness
InfectnessD[is.na(InfectnessD)] <- 0
InfnessD <- apply(InfectnessD,1,cumsum)
InfnessD <- t(InfnessD)/Sus
timeToClearance <- apply(InfnessD,1,function(x)min(which(x==max(x))))
#reg1 <- lm(timeToClearance~log(Sus))
#Outcome <- rbind(Outcome,c(InfProbSub,DeadImpact,reg1[[1]][1],reg1[[1]][2]))

timeToClearMat <- cbind(InfProbSub,DeadImpact,timeToClearance,Sus,herds[,5],ImmHerds[,2],herds[,12])

#write.table(timeToClearMat,"TimeToClearB6b.txt",append=TRUE,col.names = F,row.names = F)
 }
}

write.table(Clinicals,"ClinicalsM5E5B3.txt",sep=" ")

Clinicals[is.na(Clinicals)] <- 0

tiff(file="FigureClinicals.tiff",width=1000,height=1000,res=300,family="serif")
par(mar=c(4,4,1,1))
matplot(t(Clinicals),type="l", col="black",xlab="Days following infection",ylab="Number of clinical cases")
dev.off()











