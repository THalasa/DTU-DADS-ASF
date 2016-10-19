
 k2 <- read.table("DataDADSASF.csv",header=T,sep=";")


plot(k2$east[k2$east<8e+05],k2$north[k2$east<8e+05],ylim=c(6e+06,65e+05))

Bornholm <- which(k2$east>8e+05)

k3 <- k2[-Bornholm,]
k4 <- k3[-which(is.na(k3$chr)),]
summary(k4$herdSize)

k4$herdSize[k4$herdSize<2]<- 2

k4$LamAll[is.na(k4$LamAll)] <- 0
k4$LambdaWeaners[is.na(k4$LambdaWeaners)] <- 0
k4$LamAb <- as.numeric(as.character(k4$LamAb))
k4$LamAb[is.na(k4$LamAb)] <- 0

k4$herdSizeCat <- 1
k4$herdSizeCat[k4$herdSize>5&k4$herdSize<=300] <- 2
k4$herdSizeCat[k4$herdSize>300&k4$herdSize<=1200] <- 3
k4$herdSizeCat[k4$herdSize>1200&k4$herdSize<=2250] <- 4
k4$herdSizeCat[k4$herdSize>2250] <- 5

table(k4$herdSizeCat)

write.table(k4,"DataDADSASF2.csv",sep=";")
