#### MISC library 
####
#### Includes the following functions:
####
#### changelog
#### RandContacts
#### ReedFrost1
#### ReedFrost2
#### ReedFrost3
#### rpert
#### rgBeta
#### rtriang
####

### calculate the distance based on the Eucludian using UTM coordinates
calDist <<- function(Index){           
             distance <- sqrt((aHerd$north-aHerd$north[Index])^2 +(aHerd$east-aHerd$east[Index])^2)/1000
             return(distance)
           }

## Summary statistics
bsumm <- function(sim.vals,graph="none",rnd=0,alpha=0.05,...) {
  pct <- c(alpha/2,1-alpha/2)
  out <- c(mean(sim.vals),median(sim.vals),sqrt(var(sim.vals)),
           quantile(sim.vals,probs=pct))
  names(out) <- c("mean","median","sd",paste(pct[1]*100,"%",sep=""),
                  paste(pct[2]*100,"%",sep=""))
  if(graph=="density") plot(density(sim.vals,from=min(sim.vals)),...)
  else
    if(graph=="hist")
      hist(sim.vals,freq=F,xlab=paste("N =",length(sim.vals)),...)
  return(round(out,digits=rnd))
}

## randomize contact.vectors for fractions
RandContacts <- function(contact.vector) {
  floor(contact.vector)+rbinom(length(contact.vector),1,contact.vector -
                               floor(contact.vector))
}

## Tom's Reed-Frost function (recalculate entire intraherd epidemic curve each
## day)
ReedFrost1 <- function(T,iLag,susc,pr) {
  TotSus <- susc
  inf <- 0
  T <- max(1,T)
  if(iLag<1)
    iLag<-0
  m <- c(1,rep(0,T+iLag))
  for (x in 1:T) {
    inf <- inf+m[x]
    m[x+iLag] <- susc*(1-(1-pr)^inf)
    susc <- susc-m[x+iLag]
  }
  return(min(1,sum(m[1:T])/TotSus))
}

## Tom's version, corrected to change 0 iLag to 1; may not matter
ReedFrost2 <- function(T,iLag,susc,pr) {
  TotSus <- susc
  iLag <- iLag+(iLag==0)
  inf <- 0
  T <- max(1,T)
  m <- c(1,rep(0,50))
  for (x in 1:T) {
    inf <- inf+m[x]
    m[x+iLag] <- susc*(1-(1-pr)^inf)
    susc <- susc-m[x+iLag]
  }
  return(min(1,sum(m[1:T])/TotSus))
}

## 0/1 infectious status
ReedFrost3 <- function(T,iLag,susc,pr){
  (T > iLag)
}


## Function to generate n PERT beta random variates with min a, mode l, and max
## b.
rpert <- function(n,a,l,b) {
  mu <- (a+4*l+b)/6
  if (mu==l) v <- w <- 3 else {
    v <- (mu-a)*(2*l-a-b)/(l-mu)/(b-a)
    w <- v*(b-mu)/(mu-a)
    }
  a+(b-a)*rbeta(n,v,w)
}

## Function to generate n Generalized BETA random variates, alpha a, beta b,
## MIN min, MAX max.
rgBeta <-function(n, a, b, min, max) {
	min+((max-min)*(rbeta(n,a,b)))
}

## Function to generate n triangular distributed variates
rtriang <- function(n, min=0, mode=NULL, max=1) {
	if(is.null(mode))
          mode<-(max+min)/2
  ##Error checking
	if(min>max) {
          warning("rtriang: min>max, swx values");temp<-min;min<-max;max<-temp
        }
	if((mode<min)|(mode>max))
          stop("Mode(",mode,") outside range (",min,", ",max,")")
	if(max==min)
          return(rep(mode,n))
	##distrib creation
	xvals<-runif(n); hivals<-(xvals>(mode/(max-min)))
	xvals[!hivals]<-(min+(sqrt((xvals[!hivals])*((mode-min)*(max-min)))))
	xvals[hivals]<-(max-(sqrt((1-xvals[hivals])*((max-mode)*(max-min)))))
	return(xvals)
}


###########################
## lec.pmin
##
## A reduced version of 'pmin' taking two arguments.
## It does not checking for NA's
###########################
## lec.pmin<-function (el1,el2) {
##   change <- el1 > el2
##   el1[change]<-el2[change]
##   mostattributes(el1) <- attributes(el2)
##   el1
## }

###########################
## maxUnity
##
## Getting around pmin(1,...)
## and not checking for NA's
###########################
maxUnity<-function (el) {
  change <- el > 1
  el[el>1]<-1
  el
}


###########################################################################
## randInfoFile
##
## Generate random herd info file using:
## locations selected within a circle with approximately uniform pdf
## specified herd type frequencies (defaults from Bates 3-county data)
##
## Input:
## lat       latitude of circle center (degrees)
## long      longitude of circle center (degrees)
## radius    radius of circle (kilometers) -- length on earth surface
## dens      herd density (herds/km^2) -- only used if radius not specified
## n         number of herds to generate
## filename  output file name -- file only written if name provided
## freqs     vector of herd type frequencies
##
## Output:   nx5 data frame with herd IDs, latitude, longitude, herd type,
##           and herd status; optionally written to a file.
## 
###########################################################################
randInfoFile <- function(lat=36.4,long=-119.4,radius=NULL,dens=0.043,n=1000,
    filename=NULL,freqs=c(0.260,0.036,0.173,0.057,0.014,0.010,0.003,
    0.034,0.001,0.028,0.031,0.351,0.002)) {
  earthRad = 6371.3   # kilometers
  if (is.null(radius)) 
    arc =  acos(1- n/dens/(2*pi*earthRad^2)) else
      arc = radius/earthRad
  coords = rSphereSurf(n=n,arc=arc,center=c(long,lat),units="degrees")
  ID = 1:n
  herdType = sample(1:length(freqs),n,TRUE,freqs)    
  status = rep(1,n)
  out = data.frame(ID,lat=coords$lat,long=coords$long,herdType,status)
  if (!is.null(filename)) write.table(out,file=filename,row.names=F,sep=",")
  out
}


#######################################################################
## rSphereSurf
##
## Generate random locations within a circle on the surface of a sphere
## 
## Input:
## n        number of coordinate pairs    
## arc      arc length in radians from center to circle edge   
## center   coordinate pair (lon, lat) for center of circle
## units    units for both "center" and output
##
## Output:  nx2 data frame with longtiudes and latitudes 
##
#######################################################################
rSphereSurf <- function(n,arc,center=c(0,0),units="degrees") {
  degPerRad = 180 / pi   # degrees per radian
  if (units=="degrees") {
    center = center / degPerRad
  } else
    if (units!="radians") stop("Units not recognized.")   
  if (arc > pi) stop("Circle is too large--exceeds sphere circumference!")
  theta = runif(n,0,2*pi)  # random direction from center of circle
  r = acos(1-runif(n,0,1)*(1-cos(arc))) # random arc length
  lon1 = center[1]
  lat1 = center[2]
  lat = asin(sin(lat1)*cos(r)+cos(lat1)*sin(r)*cos(theta))
  if (cos(lat1) > .Machine$double.eps) 
    dlon = atan2(sin(theta)*sin(r)*cos(lat1),cos(r)-sin(lat1)*sin(lat)) else
      dlon = theta
  long = (lon1-dlon+pi)%%(2*pi) - pi
  coords = data.frame(long,lat)
  if (units=="degrees") coords = coords * degPerRad
  coords
}
 
##############################################################
##
## PreData is a function to prepare the data for the DADS
## models. this function is specific for the DADS model 
## as used in the Danish project (2009-2011). THH
##
###############################################################


PreData <- function(k){ #### is the name of the data set, it should be first imported into R ####

  k$herdsize<-k$Kvaeg1+k$Kvaeg2+k$Kvaeg4+k$Kvaeg5+k$Kvaeg6+k$Sows+k$Finishers+k$Goat_sheep ### a variable with the number of animals per herd

  k$xco<-substr(k$X,1,6)
  k$yco<-substr(k$X,8,14)
  k<-k[,!names(k)=="X"]

  ## CATTLE ##

  t<-k$milk=="NULL"
  tt<-k$Area=="NULL"
  k$type[t]<-0
  t<-k$milk=="milking"
  k$type[t&tt]<-1      ## non Bornholm Milking
  t<-k$milk=="not_m"
  k$type[t&tt]<-2      ## non Bornholm non-Milking
  tt<-k$Area=="B"
  t<-k$milk=="milking"
  k$type[t&tt]<-3      ## Bornholm Milking
  t<-k$milk=="not_m"
  k$type[t&tt]<-4      ## Bornholm non-Milking


  ## SWINE ##

  tt<-k$Area=="NULL"
  ttt<-k$Swine=="0"
  t<-k$Swine2=="AO"
  k$type[t&tt&ttt]<-5  ## non Bornholm conventional AO
  t<-k$Swine2=="Orne"
  k$type[t&tt&ttt]<-6  ## non Bornholm conventional Orne
  t<-k$Swine2=="OSow"
  k$type[t&tt&ttt]<-7  ## non Bornholm conventional OSow
  t<-k$Swine2=="PSow"
  k$type[t&tt&ttt]<-8  ## non Bornholm conventional PSow
  t<-k$Swine2=="OFar"
  k$type[t&tt&ttt]<-9  ## non Bornholm conventional OFar
  t<-k$Swine2=="PFar"
  k$type[t&tt&ttt]<-10 ## non Bornholm conventional PFar
  t<-k$Swine2=="OFin"
  k$type[t&tt&ttt]<-11 ## non Bornholm conventional OFin
  t<-k$Swine2=="PFin"
  k$type[t&tt&ttt]<-12 ## non Bornholm conventional PFin
  t<-k$Swine2=="OHob"
  k$type[t&tt&ttt]<-13 ## non Bornholm conventional OHob
  t<-k$Swine2=="PHob"
  k$type[t&tt&ttt]<-14 ## non Bornholm conventional PHob
  ttt<-k$Swine=="1"
  t<-k$Swine2=="AO"
  k$type[t&tt&ttt]<-15 ## non Bornholm SPF AO
  t<-k$Swine2=="Orne"
  k$type[t&tt&ttt]<-16 ## non Bornholm SPF Orne
  t<-k$Swine2=="OSow"
  k$type[t&tt&ttt]<-17 ## non Bornholm SPF OSow
  t<-k$Swine2=="PSow"
  k$type[t&tt&ttt]<-18 ## non Bornholm SPF PSow
  t<-k$Swine2=="OFar"
  k$type[t&tt&ttt]<-19 ## non Bornholm SPF OFar
  t<-k$Swine2=="PFar"
  k$type[t&tt&ttt]<-20 ## non Bornholm SPF PFar
  t<-k$Swine2=="OFin"
  k$type[t&tt&ttt]<-21 ## non Bornholm SPF OFin
  t<-k$Swine2=="PFin"
  k$type[t&tt&ttt]<-22 ## non Bornholm SPF PFin
  t<-k$Swine2=="OHob"
  k$type[t&tt&ttt]<-23 ## non Bornholm SPF OHob
  t<-k$Swine2=="PHob"
  k$type[t&tt&ttt]<-24 ## non Bornholm SPF PHob
  tt<-k$Area=="B"
  ttt<-k$Swine=="0"
  t<-k$Swine2=="AO"
  k$type[t&tt&ttt]<-25 ## Bornholm conventional AO
  t<-k$Swine2=="Orne"
  k$type[t&tt&ttt]<-26 ## Bornholm conventional Orne
  t<-k$Swine2=="OSow"
  k$type[t&tt&ttt]<-27 ## Bornholm conventional OSow
  t<-k$Swine2=="PSow"
  k$type[t&tt&ttt]<-28 ## Bornholm conventional PSow
  t<-k$Swine2=="OFar"
  k$type[t&tt&ttt]<-29 ## Bornholm conventional OFar
  t<-k$Swine2=="PFar"
  k$type[t&tt&ttt]<-30 ## Bornholm conventional PFar
  t<-k$Swine2=="OFin"
  k$type[t&tt&ttt]<-31 ## Bornholm conventional OFin
  t<-k$Swine2=="PFin"
  k$type[t&tt&ttt]<-32 ## Bornholm conventional PFin
  t<-k$Swine2=="OHob"
  k$type[t&tt&ttt]<-33 ## Bornholm conventional OHob
  t<-k$Swine2=="PHob"
  k$type[t&tt&ttt]<-34 ## Bornholm conventional PHob
  ttt<-k$Swine=="1"
  t<-k$Swine2=="AO"
  k$type[t&tt&ttt]<-35 ## Bornholm SPF AO
  t<-k$Swine2=="Orne"
  k$type[t&tt&ttt]<-36 ## Bornholm SPF Orne
  t<-k$Swine2=="OSow"
  k$type[t&tt&ttt]<-37 ## Bornholm SPF OSow
  t<-k$Swine2=="PSow"
  k$type[t&tt&ttt]<-38 ## Bornholm SPF PSow
  t<-k$Swine2=="OFar"
  k$type[t&tt&ttt]<-39 ## Bornholm SPF OFar
  t<-k$Swine2=="PFar"
  k$type[t&tt&ttt]<-40 ## Bornholm SPF PFar
  t<-k$Swine2=="OFin"
  k$type[t&tt&ttt]<-41 ## Bornholm SPF OFin
  t<-k$Swine2=="PFin"
  k$type[t&tt&ttt]<-42 ## Bornholm SPF PFin
  t<-k$Swine2=="OHob"
  k$type[t&tt&ttt]<-43 ## Bornholm SPF OHob
  t<-k$Swine2=="PHob"
  k$type[t&tt&ttt]<-44 ## Bornholm SPF PHob

  ## Sheep ##

  tt<-k$Area=="NULL"
  t<-k$sheep.goat=="Hobb"
  k$type[t&tt]<-45    ## non Bornholm Hobby
  t<-k$sheep.goat=="Comm"
  k$type[t&tt]<-46    ## non Bornholm Commercial

  tt<-k$Area=="B"
  t<-k$sheep.goat=="Hobb"
  k$type[t&tt]<-47    ## Bornholm Hobby
  t<-k$sheep.goat=="Comm"
  k$type[t&tt]<-48    ## Bornholm Commercial

  return(k)
  } # end PreData