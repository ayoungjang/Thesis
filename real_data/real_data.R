library(survival)
library(R2WinBUGS)
library(survminer)
library(dplyr)
library(fitdistrplus)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(gmodels)

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
attach(Gradient_data)
CrossTable(lab_id,Gradient_test)

print(nlevels(Gradient_data$lab_id))

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

Gradient_Etest.sub <- subset(Gradient_Etest, select = c(Strain_no, Species,lab_id, upper, lower, log.MIC))

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
X <- model.matrix(~ (Strain_no+Species)^2, data = Gradient_Etest)
n <- nrow(X);
n.beta <- ncol(X); 

Gradient_Etest$lab_id <- as.factor(Gradient_Etest$lab_id)
n.lab <- nlevels(Gradient_Etest$lab_id)

bugs.data <- with(Gradient_Etest, list(n = n, lower = lower, upper = upper, X = X, n.beta = n.beta, n.lab = n.lab, lab = as.numeric(lab_id)))


# bugs inits
bugs.inits <- function()
  list(y = with(Gradient_Etest, runif(n, lower, upper)), b.lab = rnorm(n.lab, 0, 0.1), beta = rnorm(n.beta), sigma = runif(1), sigma.lab = runif(1))

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
  lab = levels(Gradient_Etest$lab_id),
  diff.log.MIC = colMeans(b.lab),
  lower.diff.log.MIC = apply(b.lab, 2, quantile, 0.025),
  upper.diff.log.MIC = apply(b.lab, 2, quantile, 0.975))


Gradient_Etest.newdata <- with(Gradient_Etest, expand.grid(Strain_no = levels(Strain_no), Species = levels(Species))) 


# naive mean and mode
Gradient_Etest.newdata <- within(Gradient_Etest.newdata, {mode.log.MIC <- NA; E.log.MIC.naive <- NA; se.log.MIC.naive <- NA;})


k<-1
Gradient_Etest$Species <- as.factor(Gradient_Etest$Species)
Gradient_Etest$Strain_no <- as.factor(Gradient_Etest$Strain_no)


  for (i in 1:nlevels(Gradient_Etest$Strain_no)) {
    for (j in 1:nlevels(Gradient_Etest$Species)) {
      
      strain_value <- levels(Gradient_Etest$Strain_no)[i]
      species_value <- levels(Gradient_Etest$Species)[j]
   
      Gradient_Etest.data.sub <- subset(Gradient_Etest, Strain_no == strain_value & Species == species_value)
      if (nrow(Gradient_Etest.data.sub) == 0) {
        next
      }
      
      mod <- with(Gradient_Etest.data.sub, lm(log.MIC.naive ~ 1))
    
      Gradient_Etest.newdata[k, "mode.log.MIC"] <- log2(as.numeric(sub("<", "", sub("<=", "" ,sub(">", "", sub(c(">="), "", names(which.max(table(Gradient_Etest.data.sub$MIC)))))))))
      Gradient_Etest.newdata[k, c("E.log.MIC.naive", "se.log.MIC.naive")] <- c(coef(mod), sqrt(vcov(mod)))

      k <- k+1
    }
  }


Gradient_Etest.newdata <- na.omit(Gradient_Etest.newdata)

Gradient_Etest.newdata <- within(Gradient_Etest.newdata, {
  lower.log.MIC.naive <- E.log.MIC.naive - 1.96 * se.log.MIC.naive
  upper.log.MIC.naive <- E.log.MIC.naive + 1.96 * se.log.MIC.naive
})

X.samplepred <- model.matrix(~(Species+Strain_no)^2, data = Gradient_Etest.newdata)

mu_sample <- t(X.samplepred%*%t(beta))

Gradient_Etest.newdata <- within(Gradient_Etest.newdata, {
  E.log.MIC <- colMeans(mu_sample)
  lower.log.MIC <- apply(mu_sample, 2, quantile, 0.025)
  upper.log.MIC <- apply(mu_sample, 2, quantile, 0.975)
})

Gradient_Etest.newdata <- merge(Gradient_Etest.newdata, Gradient_Etest.sub, sort = F)

# labs compared to mean
pdf(file = "figure_3.pdf", width = 7, height = 7)
par(mar = c(4.5, 4.5, 0.5, 0.5), yaxs = "i")
plot.new()
plot.window(xlim = c(-2, 2), c(n.lab + 1, 0))
abline(h = 1:n.lab, col = 8, lty = 3)
box()
axis(1, at = seq(-2, 2, 1), labels = 2^seq(-2, 2, 1))
axis(2, at = 1:n.lab, labels = lab.data$lab)
abline(v = 0, lty = 2)
title(xlab = "Fold difference in MIC compared to the expected mean", ylab = "Laboratory number")
with(lab.data, {
  points(diff.log.MIC, 1:n.lab, pch = 15)
  segments(lower.diff.log.MIC, 1:n.lab, upper.diff.log.MIC, 1:n.lab)
})
dev.off()
system(paste("open", "figure_3.pdf"))


# Gradient_Etest.newdata$Species <- factor(Gradient_Etest.newdata$Strain_no, levels = unique(Gradient_Etest.newdata$Strain_no[order(as.numeric(Gradient_Etest.newdata$Species))]))
# Gradient_Etest.newdata$Strain_no <- factor(Gradient_Etest.newdata$Strain_no, levels = unique(Gradient_Etest.newdata$Strain_no[order(as.numeric(Gradient_Etest.newdata$Strain_no))]))
n <- ncol(X.samplepred)
d <- 0.25
print(n)



pdf(file = "figure_4.pdf", width = 7, height = 7 * sqrt(2))
# Gradient_Etest.newdata$Strain_no <- factor(Gradient_Etest.newdata$Strain_no, levels = unique(as.character(1:n)))

n <- ncol(Gradient_Etest.newdata)
d <- 0.25

par(mar = c(8,5.5,2, 5), yaxs = "i");
plot.new();

# print(levels(interaction(Gradient_Etest.newdata$Strain_no,Gradient_Etest.newdata$Species, sep = "")))

plot.window(xlim = c(-10, 8), c(n+0.5, 0.5))

abline(h = seq(0.5, n+0.5, 1), col = 8);
abline(h = c(10.5, 20.5), lwd = 2); box();
axis(1, at = seq(-9, 7, 2), labels = signif(2^seq(-9, 7, 2), 3), cex.axis = 0.7)
axis(2, at = 1:n, labels = with(Gradient_Etest.newdata, levels(interaction(Strain_no, sep = " "))), las = 1, cex.axis = 0.7)
# axis(2, at = 1:n, labels = with(Gradient_Etest.newdata, levels(interaction(Strain_no,Species, sep = " "))), las = 1, cex.axis = 0.7)

title(xlab = "MIC")

with(Gradient_Etest.newdata, {
  points(mode.log.MIC, 1:n, pch = 0, cex = 0.7) # open squares - mode MICs
  points(E.log.MIC, 1:n-d, pch = 15, cex = 0.7) #solid square - Predicted mean MICs and 95% confidence intervals
  segments(lower.log.MIC, 1:n-d, upper.log.MIC, 1:n-d)
  points(lower.log.MIC.ref, 1:n+d, col = 1, cex = 0.7)#open circles - lower boundary accounting for interval censoring
  points(upper.log.MIC.ref, 1:n+d, col = 1, pch = 16, cex = 0.7) #solid circles - reference MICs
  segments(lower.log.MIC.ref, 1:n+d, upper.log.MIC.ref, 1:n+d, col = 1)
})
par(xpd=TRUE)
par(new=T)
par(fig=c(0, 1, 0, 0.2), mar=c(2,2,2,2))
#lower boundary accounting for 
legend("bottomright", legend = c("Mode MICs", "mean MICs", "interval censoring", " reference MICs"),
       pch = c(0, 15, 1, 16, 1), col = c("black", "black", "black", "black"),
       cex = 0.7, bty = "n")
dev.off()
system(paste("open", "figure_4.pdf"))

