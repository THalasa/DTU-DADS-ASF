#### GRAPHICS library version 0.05
#### Includes the following functions:
####
#### updateCharts
#### cbands
#### points1
#### points2
#### epimap
#### colScale
#### filled3
#### fmdplot
####

updateCharts <- function(step) {
  if (!hideMap && step=="time") {
    epimap(aHerd$lat,aHerd$long, aHerd$status, aHerd$vaccinated, vaccToday,
           showControls=outbreakDetected,
           main=paste("iteration",iteration,"\nDay",gTime), gTime=gTime,
           sub="", Lock, aHerd$timeInfected, aHerd$herdType)
    mtext(text=c("latent","Subclinical","Clinical","Diag/slaug"),
          side=1,line=2,adj=c(0.2,0.43,0.67,0.95))
    mtext(text=c("Herds:",sum(aHerd$status==2),sum(aHerd$status==3),
            sum(aHerd$status==4),sum(aHerd$status[aHerd$herdType!=13]==5)+
            sum(aHerd$status[aHerd$herdType!=13]==6)+
            sum(aHerd$aHerd$Diagnosed)),side=1,line=3,adj=(0:4)/4.5)
    mtext(text=c("Animals:",sum(aHerd$herdSize[aHerd$status==2]),
            sum(aHerd$herdSize[aHerd$status==3]),
            sum(aHerd$herdSize[aHerd$status==4]),
            sum(aHerd$herdSize[aHerd$status==5|aHerd$status==6])),
          side=1,line=4,adj=(0:4)/4.5)
    Sys.sleep(pause)    # wait for amount of time requested
  }

  if(is.null(initHideMap) && step=="iteration")
    hideMap<-TRUE
  if (!hideMap & step=="iteration")
    dev.set(dev.next())
  if (step=="iteration" && !hidePlots && iteration>1 && (iteration%%itupdate==0
        || iteration==n)) { 
    cbands(cumns[1:iteration,1:max(gEpiDur[1:iteration])],last=(iteration<n),
           xlab="Time (days)",ylab="Cumulative Incidence (herds)",
           main=paste("Epidemic Curves for",iteration,"iterations"))
    if (!is.na(gDaysUntilBaseline)){
      abline(v=gDaysUntilBaseline,lty=3)
      text(gDaysUntilBaseline,0,"Baseline controls established",
           cex=.7,adj=c(0,1))
    }
    title(sub=paste("Simulation:",runID))
  }
  if (!hidePlots && !hideMap && step=="iteration")
    dev.set(dev.next())
}


## graphing confidence bands for cumulative distributions
cbands <- function(matx,pv=c(0.025,0.05,0.25,0.75,0.95,0.975),
                   lty=c(4,3,2,2,3,4), col=c("turquoise","blue","purple",
                                         "purple","blue","turquoise"),
                   last=TRUE,txt=TRUE,...) {
  mdim <- dim(matx)
  xstps <- 1:mdim[2]
  matq <- matrix(0,length(pv),length(xstps))
  means <- rep(0,mdim[2])
  for(j in xstps) {
    matq[,j] <- quantile(matx[,j],p=pv)
    means[j] <- mean(matx[,j])
  }
  plot(xstps,means,ylim=c(0,matq[length(pv),mdim[2]]),type="l",lwd=2,...)
  for(i in length(pv):1) lines(xstps,matq[i,],type="l",lty=lty[i],col=col[i])
  if(txt)
    legend(0,matq[length(pv),mdim[2]],c("95% prediction interval","90% prediction interval","50% prediction interval","mean"),lty=c(lty[1:3],1),col=c(col[1:3],"black"),cex=.7)
  if (last) {
    lines(xstps,matx[mdim[1],],col="red",type="s")
    if(txt) text(mdim[2]*.8,matx[mdim[1],mdim[2]]+2,"current",col="red")
  }
}

## Point Mass plot of infected over multiple iterations (sales yards omitted)
points1 <- function(aHerd,gInfByEnd,n,...) {
  crude <- gInfByEnd/n
  mc <- max(crude)
  l.col <- heat.colors(5)
  plot(aHerd$long[aHerd$herdCategory!=6 ],aHerd$lat[aHerd$herdCategory!=6],
       xlab="latitude", ylab="longitude",cex=.4)
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = "springgreen4", border = "black")
  points(aHerd$long[crude < .1*mc & aHerd$herdCategory!=6],
         aHerd$lat[crude < .1*mc & aHerd$herdCategory!=6], col=l.col[5],pch=18)
  points(aHerd$long[crude >= .1*mc & crude < .2*mc & aHerd$herdCategory!=6],
         aHerd$lat[crude >= .1*mc & crude < .2*mc & aHerd$herdCategory!=6],
         col=l.col[4],pch=18)
  points(aHerd$long[crude>=(.2*mc) & crude < (.3*mc) & aHerd$herdCategory!=6],
         aHerd$lat[crude>=(.2*mc) & crude < (.3*mc) & aHerd$herdCategory!=6],
         col=l.col[3],pch=18)
  points(aHerd$long[crude>=(.3*mc) & crude < (.4*mc) & aHerd$herdCategory!=6],
         aHerd$lat[crude >= (.3*mc) & crude < (.4*mc) & aHerd$herdCategory!=6],
         col=l.col[2],pch=18)
  points(aHerd$long[crude >= (.4*mc) & aHerd$herdCategory!=6],
         aHerd$lat[crude>=(.4*mc) & aHerd$herdCategory!=6],col=l.col[1],pch=18)
  title(main="Estimated risk of premises becoming infected",
        sub="Sales Yards Omitted")
  legend(80.9,37.15,c(paste("prob >",.4*mc),paste(.3*mc,"< p <",.4*mc),
                      paste(.2*mc,"< p <",.3*mc), paste(.1*mc,"< p <",.2*mc),
                      paste("prob <",.1*mc)),fill=l.col, cex=.65, bg="white")
}

## Point Mass plot of infected over multiple iterations (sales yards included)
points2 <- function(aHerd,gInfByEnd,n,...) {
  crude <- gInfByEnd/n
  mc <- max(crude)
  l.col <- heat.colors(5)
  plot(aHerd$long,aHerd$lat, xlab="latitude", ylab="longitude",cex=.4)
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = "springgreen4", border = "black")
  points(aHerd$long[crude<.1*mc], aHerd$lat[crude<.1*mc], col=l.col[5],pch=18)
  points(aHerd$long[crude >= .1*mc & crude < .2*mc],
         aHerd$lat[crude >= .1*mc & crude < .2*mc], col=l.col[4],pch=18)
  points(aHerd$long[crude >= (.2*mc) & crude < (.3*mc)],
         aHerd$lat[crude >= (.2*mc) & crude < (.3*mc)], col=l.col[3],pch=18)
  points(aHerd$long[crude >= (.3*mc) & crude < (.4*mc)],
         aHerd$lat[crude >= (.3*mc) & crude < (.4*mc)], col=l.col[2],pch=18)
  points(aHerd$long[crude >= (.4*mc)],aHerd$lat[crude >= (.4*mc)],
         col=l.col[1],pch=18)
  title(main="Estimated risk of premises becoming infected",
        sub="Sales Yards Included")
  legend(80.9,37.15,c(paste("prob >",.4*mc),paste(.3*mc,"< p <",.4*mc),
                      paste(.2*mc,"< p <",.3*mc),paste(.1*mc,"< p <",.2*mc),
                      paste("prob <",.1*mc)),fill=l.col, cex=.65, bg="white")
  ##legend(79.2,36.6,"Scan Cluster, p < 0.001",fill="blue",cex=0.65,bg="white")
}


## Mapping the epidemic
epimap <- function(lat, long, status, vaccinated, vaccToday, showControls,
                   gTime, Lock, InfxnTime, type, ...) {
  plot(long[status==1],lat[status==1],xlab="",ylab="latitude",bty="o",type="n",
       cex=.3,col="darkgrey",...)
  infected <- status==3 | status==4
  theta <- seq(0,2*pi,length=20)
  xcirc <- cos(theta)*0.008993288
  ycirc <- sin(theta)*0.008993288
  if (showControls) {
    tmpCol<-colScale(length(controlMethods))
    for (i in 1:length(controlMethods)){
      tmpID<-controlMethods[[i]]$getIn()
      points(long[tmpID],lat[tmpID],pch=19,col=tmpCol[i])
    }
  }
  points(long[status==1],lat[status==1],col="darkgrey",pch=".")
  points(long[vaccinated],lat[vaccinated],col="dark green",pch=20)
  points(long[(status==2)&(InfxnTime<gTime)],
         lat[(status==2)&(InfxnTime<gTime)],col="red",pch=1)

  points(long[status==5],lat[status==5],col="black",pch=3)
  points(long[type==13],lat[type==13], col="black", pch=0)
  points(long[type==13&infected],lat[type==13&infected],col="red", pch=21)
  points(long[type==13&status==5],lat[type==13&status==5], col="black", pch=12)
  points(long[infected&(InfxnTime<gTime)],lat[infected&(InfxnTime<gTime) ],
         col="maroon",pch=16)
}

colScale<-function(n,half){
  if(missing(half))
    rgb(red=rep(255,n),
        green=as.integer(seq(255,0,length=n)),
        blue=as.integer(seq(127,0,length=n)),
        alpha=rep(255,n),
        maxColorValue=255)
  else
    rgb(red=rep(255,n),
        green=as.integer(c(rep(255,floor(n*half)),
          seq(255,0,length=n-floor(n*half)))),
        blue=c(as.integer(seq(210,0,length=floor(n*half))),rep(0,n-floor(n*half))),
        alpha=rep(255,n),
        maxColorValue=255)
}


## Filled contour plot of infected over multiple iterations
## This form of mapping was not satisfactory because it predicts infection
## where there is no herd but I left it in in case Tim wated to use it anyway.
filled3 <- function(aHerd,gInfByEnd,n,...) {
  crude <- gInfByEnd/n
  gInf.loess <- loess(crude~long*lat,data=aHerd, span=0.1)
  gInf.mar<-list(long=seq(min(aHerd$long),max(aHerd$long),
                   length=gMaxHerds^0.5),
                 lat=seq(min(aHerd$lat),max(aHerd$lat),length=1+gMaxHerds^0.5))
  gInf.z <- matrix(predict(gInf.loess,expand.grid(gInf.mar)),
                   length(gInf.mar$long),length(gInf.mar$lat))
  gInf.z[gInf.z<0] <- 0
  filled.contour(gInf.mar$long, gInf.mar$lat, gInf.z, zlim=c(0,max(gInf.z)),
                 color = terrain.colors,
                 plot.title = title(main = "Estimated risk of premises becoming infected",
                   xlab = "longitude", ylab = "latitude"),
                 plot.axes={axis(1);
                            axis(2);
                            points(aHerd$long,aHerd$lat, cex=.3)},
                 key.title = title(main="Prob"))
  mtext(paste("filled.contour(.) from", R.version.string),
        side = 1, line = 4, adj = 1, cex = .65,)
}

## plotting means and printing summary statistics for multiple runs
fmdplot <- function(robj.names="fmdout",out.name="",...) {
  x11(height=5)
  n <- length(robj.names)
  means <- matrix(0,n,180)
  meanDur <- rep(0,n)
  summ <- NULL
  xstps <- 1:180
  for (i in 1:n) {
    load(paste(robj.names[i],".robj",sep=""))
    for(j in xstps) {
      means[i,j] <- mean(cumns[,j])
    }
    meanDur[i] <- mean(gEpiDur)
  }
  for (i in 1:n) {
    if (i==1)
      plot(xstps,means[i,],type="l",lty=1,lwd=2,col=1,
           ylab="Cumulative Incidence (herds)", xlab="Time (days)",
           xlim=c(1,max(meanDur)),ylim=c(1,max(means)),
           main="FMD Epidemic Curves",...)
    else
      lines(xstps,means[i,],lty=i,lwd=2,col=i,...)
  }
  legend(max(meanDur),0,robj.names,lty=1:n,xjust=1,yjust=0,col=1:n,
         bty="n",cex=.8)
  dev.copy(file=out.name,device=jpeg);dev.off()
}


