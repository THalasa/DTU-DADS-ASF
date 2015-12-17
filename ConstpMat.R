

        ConstpMat<- function(dist,DistCat,distprobA,initSize=100,addSize=25,reduce=TRUE){

          pMat <- matrix(0,nrow=dist$getNrow(),ncol=dist$getNcol())
          pMatIDs<-rep(NA,initSize) 

 
          numRow<-nrow(pMat)    # Number of rows in pMat 
          numCol<-ncol(pMat)    # Number of columns allocated in pMat 
          usedCol<-0            # Last column used in pMat 
          xCol <- 0

           fillCat<-function(distVec,distCat){
              distVec<-ceiling(distVec)-1
              dTmp<-numeric(length(distVec))+1.0
              for (i in 1:length(distCat))
                dTmp <- dTmp + as.numeric(distVec>=distCat[i])
              return(dTmp)
           }#End of fillCat

       list(

         fillProb <<- function(pMatTmp,distprobA,tmpID){
           pMatTmp[]<-fillCat(pMatTmp[1:prod(dim(pMatTmp))],probList$DistCat)
           temp1<-pMatTmp>length(probList$DistCat)
           pMatTmp[temp1]<-length(probList$DistCat)
           for(i in unique(aHerd$region)){
            tmp <- aHerd$region==i
              pMatTmp[tmp,] <- distprobA[[i]][pMatTmp[tmp,]]
            }
           pMatTmp[is.na(pMatTmp)] <-0
           ## I do not normalize, because in such a case some shipments from Bornholm will end up in the sea :-) TH
           return(pMatTmp)
          }
          ,
         addIDs<-function(IDs){
          n<-length(IDs)
          ## First check if there are enough columns
           if(n > numCol-usedCol){
           ## The following line makes sure that enough columns are added
           addN<-( (usedCol + n - numCol) %/% addSize+1 )*addSize
           pMat<<-cbind(pMat,matrix(nrow=numRow,ncol=addN))
           pMatIDs<<-c(pMatIDs,rep(NA,addN))
           numCol<<-numCol+addN 
          }
         ## Then fill columns for the requested IDs
         newCol <- usedCol+(1:n)    
         pMat[,newCol]<<- dist$get(IDs)
         pMatIDs[newCol]<<-IDs
         pMat[,newCol]<<-fillProb(pMat[,newCol,drop=FALSE],distprobA,pMatIDs[newCol])            
         usedCol<<-usedCol+n
         return(newCol)
         }# endOf addIDs
        ,
        get=function(origins){
         found<-match(origins,pMatIDs)
         miss<-is.na(found)
         if (any(miss)){
           found[miss]<-addIDs(origins[miss])
         }
         if (reduce){
           if (verbose) cat("In pMat$get: Reducing to:",length(found)," of",origins,"\n")
           pMat[,1:length(found)]<<-pMat[,found,drop=FALSE]
           pMatIDs<<-c(pMatIDs[found],rep(NA,numCol-length(found)))
           return(pMat[,1:length(found),drop=FALSE])

             }

         else
           return(pMat[,found,drop=FALSE])
        }# endOf get
         ,
        
        wipe=function(){
         pMat<<-matrix(nrow=numRow,ncol=initSize)
         pMatIDs<<-rep(NA,initSize) 
         usedCol<<-0
         numCol<<-initSize
        }# endOf wipe
       )# End of list
      }## endOf ConstpMat

