library(survival)
library(R2WinBUGS)
library(survminer)
library(dplyr)
library(fitdistrplus)
library(openxlsx)
library(ggplot2)
library(reshape2)

setwd("C:/Users/ayoung/Desktop/Thesis/real_data")
getwd()

file <- "LRE_070723_Samlet.xlsx"

data <- read.xlsx(xlsxFile=file)

data[is.na(data)] <- 0 #NA -> 0

data$Species <- as.factor(data$Species)
data$Strain_no <- as.factor(data$Strain_no)

for(i in 1:nlevels(data$Species)){
  for(j in 1:nlevels(data$Strain_no)){
    species_value <- levels(data$Species[i])
    strain_value <- levels(data$Strain_no[j]);

    data.sub<-(subset(data, Species == species_value & Strain_no == strain_value))
      
  
  }
}