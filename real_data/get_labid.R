library(survival)
library(survminer)
library(dplyr)
library(fitdistrplus)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(gmodels)
library(gdata)
library(openxlsx)
source("C:/Users/ayoung/Desktop/Thesis/real_data/plots.R")
setwd("C:/Users/ayoung/Desktop/Thesis/real_data")
getwd()

file <- "LRE_070723_Samlet.xlsx"
file2 <- "Fasit_kres_eu.xlsx"

Fasit_data<-read.xlsx(xlsxFile=file2)
Lre_data <- read.xlsx(xlsxFile=file)


names(Fasit_data) <- unname(Fasit_data[2, ])

Fasit_data <- Fasit_data[-1,]
Fasit_data <- Fasit_data[-1,]
colnames(Fasit_data) <- c("Strain_no", "LIN_mm_zone", "LIN_mm_zone_1","LIN_MIC_MTS","LIN_MIC_Etest","LIN_MIC_Etest_1")


Lre_data$Species <- as.factor(Lre_data$Species)
Lre_data$Strain_no <- as.factor(Lre_data$Strain_no)

Gradient_data <-  Lre_data[,c("lab_id", "Strain_no","Species", names(Lre_data)[grepl("Gradient", names(Lre_data))])]
Gradient_data <-  merge(Gradient_data, Fasit_data, by = "Strain_no", all.x = TRUE)

Gradient_data <- na.omit(Gradient_data)
Fasit_data <- na.omit(Fasit_data)
# write.xlsx(Gradient_data,sheetName="sheet1",file="Gradient_data_without_NA.xlsx")

Gradient_data <- within(Gradient_data, {
  MIC <- ifelse(grepl("E-test/bioMerieux", rownames(Gradient_data)), Fasit_data$`LIN_MIC_Etest`, Fasit_data$`LIN_MIC_MTS`)
  
  # event type: 3 = interval, 2 = left, 0 = right
  event <- rep(3, nrow(Gradient_data))
  event[union(grep("<", MIC), grep("<=", MIC))] <- 2
  event[union(grep(">", MIC), grep(">=", MIC))] <- 0
  # as numeric
  MIC.num <- as.numeric(sub("<", "", sub("<=", "" ,sub(">", "", sub(c(">="), "", ifelse(MIC==">256", 512,MIC))))))
  # log MIC
  log.MIC <- log2(MIC.num)
  # log MIC for naive mean
  log.MIC.naive <- ifelse(event!=0, log.MIC-0.5, log.MIC+0.5)
  # as Surv
  log.MIC.surv <- Surv(time = ifelse(event==3, log.MIC-1, log.MIC), time2 = log.MIC, event = event, type = "interval")
  # lower and upper for WinBUGS
  lower <- ifelse(event==3, log.MIC-1, ifelse(event==2, -100, log.MIC))
  upper <- ifelse(event==3, log.MIC, ifelse(event==2, log.MIC, 100))
})

Gradient_Etest<- data.frame()
Gradient_MTS <- data.frame()

Gradient_Etest <- Gradient_data[Gradient_data$Gradient_test == "E-test/bioMerieux", ]
Gradient_MTS <- Gradient_data[Gradient_data$Gradient_test != "E-test/bioMerieux", ]

Gradient_Etest <- subset(Gradient_Etest)


data$Strain_no <- as.factor(data$Strain_no)
data$Species <- as.factor(data$Species)

data$Strain_no<-drop.levels(data$Strain_no)
head(data$Strain_no)

data$Species<-drop.levels(data$Species)
head(data$Species)

Fasit_data$Strain_no<-drop.levels(Fasit_data$Strain_no)
head(Fasit_data$Strain_no)

Fasit_data$LIN_MIC_MTS <-as.numeric(sub("<", "", sub("<=", "" ,sub(">", "", sub(c(">="), "", ifelse(Fasit_data$LIN_MIC_MTS==">256", 512,Fasit_data$LIN_MIC_MTS))))))
# target_strain <- 2 #1,3,5,6,20 
# arr <- c(1,4,5,6,20);
# 
# for(i in 1:length(arr)){
  target_strain <- 8

  
  data_strain <- subset(data, Strain_no == target_strain)
  
  data_strain.newdata <- with(data_strain, expand.grid(lab_id = levels(lab_id)))
  
  data_strain.sub <- within(data_strain.newdata, {lower.log.MIC.ref <- NA; upper.log.MIC.ref <- NA;})
  data_strain.newdata <- within(data_strain.newdata, {mode.log.MIC<-NA; E.log.MIC.naive <- NA; se.log.MIC.naive <- NA;})
  
  data_strain$Strain_no<-drop.levels(data_strain$Strain_no)
  head(data_strain$Strain_no)
  
  
  data_strain$lab_id<-drop.levels(data_strain$lab_id)
  head(data_strain$lab_id)
  
  c_value <- subset(Fasit_data, Strain_no == target_strain)
  
  file_path <-(paste0("figure_Strain_",target_strain,".pdf"))
  pdf(file=file_path)
  MIC_frequency <- table(data_strain$MIC.num)
  barplot(MIC_frequency, 
          main = paste0("MIC Frequency Strain_",target_strain,"/ Correct Answer = " , c_value$LIN_MIC_Etest), 
          xlab = "MIC", 
          ylab = "Frequency",
          col = "skyblue"  
  )
  dev.off();
  system(paste("open",file_path))
# }
# 
# k<-1
# up<-512
# 
# for (i in 1:nlevels(data_strain$lab_id)) {
#   set.seed(i)
#   lab_id_v <- levels(data_strain$lab_id)[i]
# 
#   
#   data_strain.data.sub <- subset(data_strain, lab_id == lab_id_v)
# 
#   mod <- with(data_strain.data.sub, lm(log.MIC.naive ~ 1))
#   
#   Y_act <- rnorm(n=1,mean=c_value$LIN_MIC_Etest,sd=c_value$LIN_MIC_Etest*2)
#   if(Y_act <0) Y_act <- Y_act*-1
#   
#   
#   data_strain.sub[k,"lower.log.MIC.ref"]<-2^floor(log2(Y_act))
#   data_strain.sub[k,"upper.log.MIC.ref"]<- 2^ceiling(log2(max(Y_act,up)))
# 
#   data_strain.newdata[k, "mode.log.MIC"] <- log2(data_strain.data.sub$MIC.num)
#   data_strain.newdata[k, c("E.log.MIC.naive", "se.log.MIC.naive")] <- c(coef(mod), sqrt(vcov(mod)))
#  #sample variance =0 or difference =0
#    k <- k+1
# }
# summary(mod)

# 
# data_strain.newdata <- within(data_strain.newdata, {
#   lower.log.MIC.naive <- E.log.MIC.naive - 1.96 * se.log.MIC.naive
#   upper.log.MIC.naive <- E.log.MIC.naive + 1.96 * se.log.MIC.naive
# })
# 
# X.strain <- model.matrix(~(lab_id)^2, data = data_strain.newdata)
# 
# mu_strain <- t(X.strain%*%t(beta))
# 
# data_strain.newdata <- within(data_strain.newdata, {
#   E.log.MIC <- colMeans(mu_strain)
#   lower.log.MIC <- apply(mu_strain, 2, quantile, 0.025)
#   upper.log.MIC <- apply(mu_strain, 2, quantile, 0.975)
# })
# 
# data_strain.newdata <- merge(data_strain.newdata, data_strain.sub, sort = F)

# n <- ncol(X.strain)
# d <- 0.25
#

