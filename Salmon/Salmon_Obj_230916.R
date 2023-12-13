# Martin fish allergy data. 
# Version 220916
library(fitdistrplus)
library(readxl)
library(flexsurv)
# Reads data from clipboard
#Plotting raw data. Salmon objective 

setwd("C:/Users/ayoung/Desktop/Thesis")
getwd()

data <- read_excel("Martin_23.09.16.xlsx")

# data<-read.table("clipboard",dec=",")

salmonobj<-data[,2:3]
salmonobj<-na.omit(salmonobj)
fish1<-(salmonobj) 
colnames(fish1)<- c("left","right") #change colume name
fish1<-fish1[order(fish1$left,fish1$right),]

fish1<-data.frame(
  left=fish1$left, right=fish1$right
)


fish.norm<-fitdistcens(fish1,"norm")
pdf(file = "result/salmon_distribution.pdf", width = 7, height = 7*sqrt(2))

plotdistcens(fish1,Turnbull=FALSE)

cdfcompcens(list(fish.norm),main="Salmon. Objective reactions",xlab="Mg protein",ylab="Cumulative proportion of responses")

dev.off()
# data<-read.table("clipboard",dec=",")
# salmonobj<-data[,2:3]

salmonobj<-na.omit(salmonobj)
salmonobjlog<-log10(salmonobj)

colnames(salmonobjlog)<-c("left","right")
fish2<-salmonobjlog
fish2<-fish2[order(fish2$left,fish2$right),]

fish2 <- lapply(fish2, function(x) ifelse(is.infinite(x) & x < 0, fish2$right, x)) # change -inf to right value

fish2<-data.frame(
  left=fish2$left,
  right=fish2$right
)



fish2.norm<-fitdistcens(fish2,"norm")
pdf(file = "result/salmon_log_distribution.pdf", width = 7, height = 7*sqrt(2))
plotdistcens(fish2,Turnbull=FALSE)
#Plotting log-normal data. Salmon objective

cdfcompcens(list(fish.norm),main="Salmon. Objective reactions",xlab="Mg(log10) protein",ylab="Cumulative proportion of responses")

dev.off()

#Bootstrap. Salmon objective log-normal
salmonobjb1<-bootdistcens(fish.norm,niter=1001)
#salmonobjb1
#summary(salmonobjb1)
#plot(salmonobjb1)
#quantile(salmonobjb1)

CIcdfplot(salmonobjb1,main="Salmon. Objective reactions",xlab="Mg(log10) protein",ylab="Cumulative proportion of responses",CI.output="probability",CI.type="two.sided",CI.level=0.95,CI.col="black",CI.lty=2,CI.fill=NULL,CI.only=FALSE)

#CIcdfplot(salmonobjb1,main="Salmon. Subjective reactions",xlab="Mg protein",ylab="Cumulative proportion of responses", xaxt="n" ,  CI.output="probability",CI.type="two.sided",CI.level=0.95,CI.col="black",CI.lty=2,CI.fill=NULL,CI.only=FALSE)

#axis(1, at=0:4, labels=c("1","10","100","1000","10000"))



#Collecting 95% CI ED10 boundaries

plotdistcens(fish2,ylim=c(0,1),las=1)
points(ci950s[,1],pnorm(ci950s[,1],fish.norm$estimate[1],fish.norm$estimate[2]),
       xlab="Log Salmon Objective (mg)",ylab="",ylim=c(0,1),type="l",col="red",lwd=2,las=1)
points(ci950s[,1],ci950s[,2],type="l",lwd=1)
points(ci950s[,1],ci950s[,3],type="l",lwd=1)


sample01<-qnorm(0.1,fish.norm$estimate[1],fish.norm$estimate[2])
low01<-which(ci950s[,2]>0.0975 & ci950s[,2]<0.1005)
upp01<-which(ci950s[,3]>0.0975 & ci950s[,3]<0.1005)
l01<-min(ci950s[low01,1])
u01<-min(ci950s[upp01,1])
sample01

#Plotting raw data.Salmon subjective
# data<-read.table("clipboard",dec=",")

salmonsub<-data[,4:5]
salmonsub<-na.omit(salmonsub)
fish3<-(salmonsub) 

colnames(fish3)<- c("left","right")

fish3<-data.frame(
  left=fish3$left,
  right=fish3$right
)
fish3<-fish3[order(fish3$left,fish3$right),]

fish.norm<-fitdistcens(fish3,"norm")

pdf(file = "result/salmon_sub_distribution.pdf", width = 7, height = 7*sqrt(2))

plotdistcens(fish3,Turnbull=FALSE)
cdfcompcens(list(fish.norm),main="Salmon. Subjective reactions",xlab="Mg protein",ylab="Cumulative proportion of responses")


#Plotting log-normal data. Salmon subjective
# data<-read.table("clipboard",dec=",")
# salmonsub<-data[,3:4]
salmonsub<-na.omit(salmonsub)
salmonsublog<-log10(salmonsub)

colnames(salmonsublog)<-c("left","right")

fish4<-salmonsublog
fish4<-fish4[order(fish4$left,fish4$right),]
fish4 <- lapply(fish4, function(x) ifelse(is.infinite(x) & x < 0, fish4$right, x)) # change -inf to right value

fish4<-data.frame(
  left=fish4$left,
  right=fish4$right
)


fish.norm<-fitdistcens(fish4,"norm")
plotdistcens(fish4,Turnbull=FALSE)
cdfcompcens(list(fish.norm),main="Salmon. Subjective reactions",xlab="Mg(log10) protein",ylab="Cumulative proportion of responses")
dev.off()

#Bootstrap. Salmon subjective log-normal
salmonsubb1<-bootdistcens(fish.norm,niter=1001)
salmonsubb1
summary(salmonsubb1)
plot(salmonsubb1)
quantile(salmonsubb1)
CIcdfplot(salmonsubb1,main="Salmon. Subjective reactions",xlab="Mg(log10) protein",ylab="Cumulative proportion of responses",CI.output="probability",CI.type="two.sided",CI.level=0.95,CI.col="black",CI.lty=2,CI.fill=NULL,CI.only=FALSE)

#Plotting raw data. Mackerel objective 
# data<-read.table("clipboard",dec=",")

mackerelobj<-data[,5:6]
mackerelobj<-na.omit(mackerelobj)
fish5<-(mackerelobj) 
colnames(fish5)<- c("left","right")
fish5<-fish5[order(fish5$left,fish5$right),]
fish.norm<-fitdistcens(fish5,"norm")
plotdistcens(fish5,Turnbull=FALSE)
cdfcompcens(list(fish.norm),main="Mackerel. Objective reactions",xlab="Mg protein",ylab="Cumulative proportion of responses")

#Plotting log-normal data. Mackerel objective
# data<-read.table("clipboard",dec=",")
mackerelobj<-data[,5:6]
mackerelobj<-na.omit(mackerelobj)
mackerelobjlog<-log10(mackerelobj)
colnames(mackerelobjlog)<-c("left","right")
fish6<-mackerelobjlog
fish6<-fish6[order(fish6$left,fish6$right),]
fish.norm<-fitdistcens(fish6,"norm")
plotdistcens(fish6,Turnbull=FALSE)
cdfcompcens(list(fish.norm),main="Mackerel. Objective reactions",xlab="Mg(log10) protein",ylab="Cumulative proportion of responses")

#Bootstrap. Mackerel objective log-normal
mackerelobjb1<-bootdistcens(fish.norm,niter=1001)
mackerelobjb1
summary(mackerelobjb1)
plot(mackerelobjb1)
quantile(mackerelobjb1)
CIcdfplot(mackerelobjb1,main="Mackerel. Objective reactions",xlab="Mg(log10) protein",ylab="Cumulative proportion of responses",CI.output="probability",CI.type="two.sided",CI.level=0.95,CI.col="black",CI.lty=2,CI.fill=NULL,CI.only=FALSE)

#Plotting raw data. Mackerel subjective 
# data<-read.table("clipboard",dec=",")
mackerelsub<-data[,7:8]
mackerelsub<-na.omit(mackerelsub)
fish7<-(mackerelsub) 
colnames(fish7)<- c("left","right")
fish7<-fish7[order(fish7$left,fish7$right),]
fish.norm<-fitdistcens(fish7,"norm")
plotdistcens(fish7,Turnbull=FALSE)
cdfcompcens(list(fish.norm),main="Mackerel. Subjective reactions",xlab="Mg protein",ylab="Cumulative proportion of responses")

#Plotting log-normal data. Mackerel subjective
# data<-read.table("clipboard",dec=",")
mackerelsub<-data[,7:8]
mackerelsub<-na.omit(mackerelsub)
mackerelsublog<-log10(mackerelsub)
colnames(mackerelsublog)<-c("left","right")
fish8<-mackerelsublog
fish8<-fish8[order(fish8$left,fish8$right),]
fish.norm<-fitdistcens(fish8,"norm")
plotdistcens(fish8,Turnbull=FALSE)
cdfcompcens(list(fish.norm),main="Mackerel. Subjective reactions",xlab="Mg(log10) protein",ylab="Cumulative proportion of responses")

#Bootstrap. Mackerel subjective log-normal
mackeresubb1<-bootdistcens(fish.norm,niter=1001)
mackerelsubb1
summary(mackerelsubb1)
plot(mackeresubb1)
quantile(mackerelsubb1)
CIcdfplot(mackerelosubb1,main="Mackerel. Subjective reactions",xlab="Mg(log10) protein",ylab="Cumulative proportion of responses",CI.output="probability",CI.type="two.sided",CI.level=0.95,CI.col="black",CI.lty=2,CI.fill=NULL,CI.only=FALSE)

#Plotting raw data. Cod objective 
# data<-read.table("clipboard",dec=",")
codobj<-data[,9:10]
codobj<-na.omit(codobj)
fish9<-(codobj) 
colnames(fish9)<- c("left","right")
fish9<-fish9[order(fish9$left,fish9$right),]
fish.norm<-fitdistcens(fish9,"norm")
plotdistcens(fish1,Turnbull=FALSE)
cdfcompcens(list(fish.norm),main="Cod. Objective reactions",xlab="Mg protein",ylab="Cumulative proportion of responses")

#Plotting log-normal data. Cod objective
# data<-read.table("clipboard",dec=",")
codobj<-data[,9:10]
codobj<-na.omit(codobj)
codobjlog<-log10(codobj)
colnames(codobjlog)<-c("left","right")
fish10<-codobjlog
fish10<-fish10[order(fish10$left,fish10$right),]
#plot(fish11)
fish.norm<-fitdistcens(fish10,"norm")
plotdistcens(fish10,Turnbull=FALSE)
cdfcompcens(list(fish.norm),main="Cod. Objective reactions",xlab="Mg(log10) protein",ylab="Cumulative proportion of responses")

#Bootstrap. Cod objective log-normal
codobjb1<-bootdistcens(fish.norm,niter=1001)
codobjb1
summary(codobjb1)
plot(codobjb1)
quantile(codobjb1)
CIcdfplot(codobjb1,main="Cod. Objective reactions",xlab="Mg(log10) protein",ylab="Cumulative proportion of responses",CI.output="probability",CI.type="two.sided",CI.level=0.95,CI.col="black",CI.lty=2,CI.fill=NULL,CI.only=FALSE)

#Plotting raw data.Cod subjective
data<-read.table("clipboard",dec=",")
codsub<-data[,11:12]
codsub<-na.omit(codsub)
fish11<-(codsub) 
colnames(fish11)<- c("left","right")
fish11<-fish11[order(fish11$left,fish11$right),]
fish.norm<-fitdistcens(fish11,"norm")
plotdistcens(fish11,Turnbull=FALSE)
cdfcompcens(list(fish.norm),main="Cod. Subjective reactions",xlab="Mg protein",ylab="Cumulative proportion of responses")

#Plotting log-normal data. Cod subjective
data<-read.table("clipboard",dec=",")
codsub<-data[,11:12]
codsub<-na.omit(codsub)
codsublog<-log10(codsub)
colnames(codsublog)<-c("left","right")
fish12<-codsublog
fish12<-fish12[order(fish12$left,fish12$right),]
plot(fish12)
fish.norm<-fitdistcens(fish12,"norm",ylim=(0-1))
plotdistcens(fish12,Turnbull=FALSE)
cdfcompcens(list(fish.norm),main="Cod. Objective reactions",xlab="Mg(log10) protein",ylab="Cumulative proportion of responses")

#Bootstrap. Cod subjective log-normal
codsubb1<-bootdistcens(fish.norm,niter=1001)
codsubb1
summary(codcsubb1)
plot(codsubb1)
quantile(codsubb1)
CIcdfplot(codsubb1,main="Cod. Objective reactions",xlab="Mg(log10) protein",ylab="Cumulative proportion of responses",CI.output="probability",CI.type="two.sided",CI.level=0.95,CI.col="black",CI.lty=2,CI.fill=NULL,CI.only=FALSE)
