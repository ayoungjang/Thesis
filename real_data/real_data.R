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
source("C:/Users/ayoung/Desktop/Thesis/real_data/winbugs.R")
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



execute_winbugs(Gradient_Etest,"Etest")

execute_winbugs(Gradient_MTS,"MTS")

