### DIST library 
###
### Used to handle distance calculations
###
### Is basically a list containing functions
### and the needed variables.
###
### As part of initialization call: constructDist(coords)

######################
### constructDist
###
### Call once to construct distance object.
###
### north :  Latitude or UTM northings  (See 'norm')
### east  :  longitude or UTM eastings  (See 'norm')
### initSize : Initial number of columns in distance matrix
### addSize  : Number of columns to add when growing distance matrix
### norm     : Choose norm to be used for distance calculation
###              1) True spherical
###              2) Euclidian on lat long
###              3) Pure Euclidian (To be used with UTM coordinates)
######################
constructDist<-function(north,east,initSize=200,addSize=100,norm=3){

  ## Internal variables:
  ##
  ## VarDef distMat Distance matrix
  ## VarDes Matrix with one row per location and
  ## VarDes adjustable number of columns.
  distMat<-matrix(nrow=length(north),ncol=initSize)
  ## VarDef distIDs Distance IDs
  ## VarDes Vector of IDs for what is in the columns of distMat.
  distIDs<-rep(NA,initSize) #Mayby not NA

  numRow<-nrow(distMat) # Number of rows in distMat
  numCol<-ncol(distMat) # Number of columns allocated in distMat
  usedCol<-0            # Last column used in distMat

  rads <- 57.2958 # Conversion from degrees to radians
  meanNorth<-NA   # Only used for Eucludian using lat long
  rad2km<-1.609334 * 3958.75  # miles2km * rad2miles
  
  ## Preprocessing for calculation of distances:
  switch(norm,{
    ##1: Spherical  north=latitude, east=longitude
    north<-north/rads
    east<-east/rads
  },{
    ##2: Euclidean using lat long
    north<-north/rads
    east<-east/rads
    meanNorth<- mean(north)
  },{
    ##3: Eucludian using UTM (rectangular coordinates)
  })

  ## Defining a distance function (Tariq changed the name of this function to allow the use of the same function as adjusted beneath)
  calcD<-
    switch(norm,
           ##1: Spherical  north=latitude, east=longitude
           function(Index) {
             n <- length(Index)
             north2 <- t(north[Index])
             east2 <- east[Index]
             x <- sin(north) %*% sin(north2) + (cos(north) %*% cos(north2)) *
               cos((matrix(east,numRow,n,byrow=FALSE)-
                    matrix(east2,numRow,n,byrow=TRUE)))
             temp <- 1-x^2
             distance <- rad2km * atan(sqrt(temp*(temp>0)))/x
            
             return(distance)
           },
           ##2: Euclidean using lat long
           function(Index) {
             n <- length(Index)
             distance <- matrix(0,numRow,n)
             for (i in 1:n) {
               north2 <- north[Index[i]]
               east2 <- east[Index[i]]
               distance[,i] <- rad2km *
                 sqrt((north-north2)^2 +(cos(meanNorth)*(east-east2))^2) 
             }
             return(distance)
           },
           ##3: Eucludian using UTM (rectangular coordinates)
           function(Index){
             n <- length(Index)
             distance <- matrix(0,numRow,n)
             for (i in 1:n) {
               distance[,i] <-
                 sqrt((north-north[Index[i]])^2 +(east-east[Index[i]])^2)/1000
             }
             return(distance)
           }
           )
  ## Calculate the distances from IDs to all other locations
  ## and put them into distMat and distIDs as approiate
  addIDs<-function(IDs){
    n<-length(IDs)
    ## First check if there is enough columns
    if(n > numCol-usedCol){
      ## The following line makes sure that enough columns are added
      addN<-( (usedCol + n - numCol) %/% addSize+1 )*addSize
      distMat<<-cbind(distMat,matrix(nrow=numRow,ncol=addN))
      distIDs<<-c(distIDs,rep(NA,addN))
      numCol<<-numCol+addN
    }
    ## Then fill columns for the requested IDs
    newCol <- usedCol+(1:n)
    distMat[,newCol]<<- calcD(IDs)
    distIDs[newCol]<<-IDs
    usedCol<<-usedCol+n
    return(newCol)
  }# endOf addIDs

  list( # List holding all information to be access from outside
       
       #####################
       ## get
       ## origins: Vector of IDs used as origins 
       #####################
       get=function(origins){
         found<-match(origins,distIDs)
         miss<-is.na(found)
         if (any(miss)){
           found[miss]<-addIDs(origins[miss])
         }
         return(distMat[,found,drop=FALSE])
       }# endOf get
       ,
       getNrow=function(){
         return(numRow)
       },
       getNcol=function(){
         return(numCol)
       }
       ,
       wipe=function(){
         if (verbose) print(c("wiping ",usedCol))
         distMat<<-matrix(nrow=numRow,ncol=initSize)
         distIDs<<-rep(NA,initSize) 
         usedCol<<-0
         numCol<<-initSize
       }
       )# End of list

}## endOf constructDist


#################################
## Comments:
##
## - Consider when to remove origin locations (columns) from the matrix
##   -> When culled? 
##   -> When no longer infectious?
##   -> Or ...?
