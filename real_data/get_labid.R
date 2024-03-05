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

# working dir
wbwd <- file.path(getwd(), "WinBUGS")

# bugs model
cat("model {
  for (i in 1:n) {
     y[i] ~ dnorm(mu[i], tau)I(lower[i], upper[i])
     mu[i] <- inprod2(X[i, ], beta[]) + b.lab[lab[i]]
  }
  for (j in 1:n.lab) {
    b.lab[j] ~ dnorm(0.0, tau.lab)
  }
  for (k in 1:n.beta) {
    beta[k] ~ dnorm(0.0, 1.0E-4)
  }
  tau <- pow(sigma, -2)
  sigma ~ dunif(0.01, 100)
  tau.lab <- pow(sigma.lab, -2)
  sigma.lab ~ dunif(0.01, 100)
  }
}", file = file.path(wbwd, "model.txt"))


# bugs data
X <- model.matrix(~ (lab_id)^2, data = data)
n <- nrow(X);
n.beta <- ncol(X); 

data$lab_id <- as.factor(data$lab_id)
n.lab <- nlevels(data$lab_id)

bugs.data <- with(data, list(n = n, lower = lower, upper = upper, X = X, n.beta = n.beta, n.lab = n.lab, lab = as.numeric(lab_id)))


# bugs inits
bugs.inits <- function()
  list(y = with(data, runif(n, lower, upper)), b.lab = rnorm(n.lab, 0, 0.1), beta = rnorm(n.beta), sigma = runif(1), sigma.lab = runif(1))

# bugs fit
bugs.fit <- bugs(model.file = "model.txt", data = bugs.data, inits = bugs.inits,
                 parameters.to.save = c("beta", "b.lab", "sigma", "sigma.lab"),
                 n.chains = 2, n.iter = 5100, n.burnin = 100, n.thin = 10, debug = F, DIC = F, working.directory = wbwd)


read.bugsfit <- T
if (read.bugsfit) {
  old.wd <- getwd()
  setwd(wbwd)
  bugs.fit <<- R2WinBUGS:::bugs.sims(n.chains = 2, n.iter = 5100, n.burnin = 100, n.thin = 10, parameters.to.save = c("beta", "b.lab", "sigma", "sigma.lab"), DIC = F)
  class(bugs.fit) <- "bugs"
  setwd(old.wd)
}
attach.bugs(bugs.fit)

# labs compared to mean
lab.data <- data.frame(
  lab = levels(data$lab_id),
  diff.log.MIC = colMeans(b.lab),
  lower.diff.log.MIC = apply(b.lab, 2, quantile, 0.025),
  upper.diff.log.MIC = apply(b.lab, 2, quantile, 0.975))

Fasit_data$LIN_MIC_MTS <-as.numeric(sub("<", "", sub("<=", "" ,sub(">", "", sub(c(">="), "", ifelse(Fasit_data$LIN_MIC_MTS==">256", 512,Fasit_data$LIN_MIC_MTS))))))

target_strain <- 2 #1,3,5,6,20 

data_strain <- subset(data, Strain_no == target_strain)

data_strain.newdata <- with(data_strain, expand.grid(lab_id = levels(lab_id)))

data_strain.sub <- within(data_strain.newdata, {lower.log.MIC.ref <- NA; upper.log.MIC.ref <- NA;})
data_strain.newdata <- within(data_strain.newdata, {mode.log.MIC<-NA; E.log.MIC.naive <- NA; se.log.MIC.naive <- NA;})

data_strain$Strain_no<-drop.levels(data_strain$Strain_no)
head(data_strain$Strain_no)


data_strain$lab_id<-drop.levels(data_strain$lab_id)
head(data_strain$lab_id)

c_value <- subset(Fasit_data, Strain_no == target_strain)

k<-1
up<-512

for (i in 1:nlevels(data_strain$lab_id)) {
  set.seed(i)
  lab_id_v <- levels(data_strain$lab_id)[i]

  
  data_strain.data.sub <- subset(data_strain, lab_id == lab_id_v)

  mod <- with(data_strain.data.sub, lm(log.MIC.naive ~ 1))
  
  Y_act <- rnorm(n=1,mean=c_value$LIN_MIC_MTS,sd=c_value$LIN_MIC_MTS*2)
  if(Y_act <0) Y_act <- Y_act*-1
  
  
  data_strain.sub[k,"lower.log.MIC.ref"]<-2^floor(log2(Y_act))
  data_strain.sub[k,"upper.log.MIC.ref"]<- 2^ceiling(log2(max(Y_act,up)))

  data_strain.newdata[k, "mode.log.MIC"] <- log2(data_strain.data.sub$MIC.num)
  data_strain.newdata[k, c("E.log.MIC.naive", "se.log.MIC.naive")] <- c(coef(mod), sqrt(vcov(mod)))
 #sample variance =0 or difference =0
   k <- k+1
}
summary(mod)


data_strain.newdata <- within(data_strain.newdata, {
  lower.log.MIC.naive <- E.log.MIC.naive - 1.96 * se.log.MIC.naive
  upper.log.MIC.naive <- E.log.MIC.naive + 1.96 * se.log.MIC.naive
})

X.strain <- model.matrix(~(lab_id)^2, data = data_strain.newdata)

mu_strain <- t(X.strain%*%t(beta))

data_strain.newdata <- within(data_strain.newdata, {
  E.log.MIC <- colMeans(mu_strain)
  lower.log.MIC <- apply(mu_strain, 2, quantile, 0.025)
  upper.log.MIC <- apply(mu_strain, 2, quantile, 0.975)
})

data_strain.newdata <- merge(data_strain.newdata, data_strain.sub, sort = F)


n <- ncol(X.strain)
d <- 0.25

file_path <-(paste0("figure_Strain_",target_strain,".pdf"))
pdf(file = file_path, width = 7, height = 7 * sqrt(2))
par(mar = c(8,5.5,2, 5), yaxs = "i");
plot.new();
plot.window(xlim = c(-10, 8), c(n+0.5, 0.5))

abline(h = seq(0.5, n+0.5, 1), col = 8);
abline(h = c(10.5, 20.5), lwd = 2); box();
axis(1, at = seq(-9, 7, 2), labels = signif(2^seq(-9, 7, 2), 3), cex.axis = 0.7)
axis(2, at = 1:n, labels = with(data_strain.newdata, levels(interaction(lab_id, sep = " "))), las = 1, cex.axis = 0.7)
title(xlab = "MIC")

with(data_strain.newdata, {
  points(mode.log.MIC, 1:n, pch = 0, cex = 0.7) #open squares - mod MICs
  points(E.log.MIC, 1:n-d, pch = 15, cex = 0.7) #solid square - Predicted mean MICs and
  segments(lower.log.MIC, 1:n-d, upper.log.MIC, 1:n-d)
  points(lower.log.MIC.ref, 1:n+d, col = 1, cex = 0.7) #open circles - lower boundary accounting for interval censoring
  # points(upper.log.MIC.ref, 1:n+d, col = 1, pch = 16, cex = 0.7) # solid circles - reference MICs
  # segments(lower.log.MIC.ref, 1:n+d, upper.log.MIC.ref, 1:n+d, col = 1)
})
par(xpd=TRUE)
par(new=T)
par(fig=c(0, 1, 0, 0.2), mar=c(2,2,2,2))
#lower boundary accounting for 
legend("bottomright", legend = c("Mode MICs", "mean MICs", "interval censoring", " reference MICs"),
       pch = c(0,
               15, 1, 16, 1), col = c("black", "black", "black", "black"),
       cex = 0.7, bty = "n")
dev.off()
system(paste("open",file_path))
