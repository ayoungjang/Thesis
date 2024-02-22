# Adjusted code for Han de Neeling data
library(survival)
library(R2WinBUGS)
library(fitdistrplus)
library(ggplot2)
library(dplyr)
library(tidyr)
setwd("C:/Users/ayoung/Desktop/Thesis/RCode")
getwd()

# read data
gono.data <- read.table("gono_in.txt", header = T, sep = "\t") # gonokokken data
refmic.data <- read.table("reference MICs_U.txt", header = T, sep = "\t") # reference MICs

# reshape gono.data into long format
gono.data <- reshape(data = gono.data,
  varying = names(gono.data)[!is.element(names(gono.data), c("lab", "ant", "agar"))],
  v.names = "MIC", timevar = "stam",
  times = names(gono.data)[!is.element(names(gono.data), c("lab", "ant", "agar"))],
  direction = "long")

# create censored data objects
gono.data <- within(gono.data, {
  # as character
  MIC <- as.character(MIC)
  # event type: 3 = interval, 2 = left, 0 = right
  event <- rep(3, nrow(gono.data))
  event[union(grep("<", MIC), grep("<=", MIC))] <- 2
  event[union(grep(">", MIC), grep(">=", MIC))] <- 0
  # as numeric
  MIC.num <- as.numeric(sub("<", "", sub("<=", "" ,sub(">", "", sub(c(">="), "", MIC)))))
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

refmic.data <- within(refmic.data, {
  # as character
  MIC <- as.character(MIC)
  # event type: 3 = interval, 2 = left, 0 = right
  event <- rep(3, nrow(refmic.data))
  event[union(grep("<", MIC), grep("<=", MIC))] <- 2
  event[union(grep(">", MIC), grep(">=", MIC))] <- 0
  # as numeric
  MIC.num <- as.numeric(sub("<", "", sub("<=", "" ,sub(">", "", sub(c(">="), "", MIC)))))
  # log MIC
  log.MIC <- ceiling(log2(MIC.num))
  # as Surv
  log.MIC.surv <- Surv(time = ifelse(event==3, log.MIC-1, log.MIC), time2 = log.MIC, event = event, type = "interval")
  # lower and upper for WinBUGS
  lower <- ifelse(event==3, log.MIC-1, ifelse(event==2, -100, log.MIC))
  upper <- ifelse(event==3, log.MIC, ifelse(event==2, log.MIC, 100))
})

# "stam" -> "strain" & lab, stam as factor
gono.data <- within(gono.data, {
  stam <- gsub("stam", "strain", stam)
  lab <- factor(lab)
  stam <- factor(stam)
})

refmic.data <- within(refmic.data, {
  stam <- gsub("stam", "strain", stam) #replace
  stam <- factor(stam)
})

# remove cef data
gono.data <- droplevels(subset(gono.data, ant!="cef"))
refmic.data <- droplevels(subset(refmic.data, ant!="cef"))

# selection of refmic.data
refmic.data.sub <- subset(refmic.data, select = c(stam, ant, upper, lower, log.MIC))
names(refmic.data.sub) <- c("stam", "ant", "upper.log.MIC.ref", "lower.log.MIC.ref", "log.MIC.ref")

#
# model with WinBUGS
#

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
X <- model.matrix(~ (stam+ant)^2, data = gono.data)
n <- nrow(X); n.beta <- ncol(X); n.lab <- 17
bugs.data <- with(gono.data, list(n = n, lower = lower, upper = upper, X = X, n.beta = n.beta, n.lab = n.lab, lab = as.numeric(lab)))

# bugs inits
bugs.inits <- function()
  list(y = with(gono.data, runif(n, lower, upper)), b.lab = rnorm(n.lab, 0, 0.1), beta = rnorm(n.beta), sigma = runif(1), sigma.lab = runif(1))

# bugs fit
bugs.fit <- bugs(model.file = "model.txt", data = bugs.data, inits = bugs.inits,
  parameters.to.save = c("beta", "b.lab", "sigma", "sigma.lab"),
  n.chains = 2, n.iter = 5100, n.burnin = 100, n.thin = 10, debug = F, DIC = F, working.directory = wbwd)

#
# post process
# first run "fit with WinBUGS" part except bugs.fit
#

# read simulations
read.bugsfit <- T
if (read.bugsfit) {
  old.wd <- getwd()
  setwd(wbwd)
  bugs.fit <<- R2WinBUGS:::bugs.sims(n.chains = 2, n.iter = 5100, n.burnin = 100, n.thin = 10,
    parameters.to.save = c("beta", "b.lab", "sigma", "sigma.lab"), DIC = F)
  class(bugs.fit) <- "bugs"
  setwd(old.wd)
}
attach.bugs(bugs.fit)

# labs compared to mean
lab.data <- data.frame(
  lab = levels(gono.data$lab),
  diff.log.MIC = colMeans(b.lab),
  lower.diff.log.MIC = apply(b.lab, 2, quantile, 0.025),
  upper.diff.log.MIC = apply(b.lab, 2, quantile, 0.975))

# make refmic.newdata
refmic.newdata <- with(gono.data, expand.grid(stam = levels(stam), ant = levels(ant))) 


# naive mean and mode
refmic.newdata <- within(refmic.newdata, {mode.log.MIC <- NA; E.log.MIC.naive <- NA; se.log.MIC.naive <- NA;})


k<-1
gono.data$ant <- as.factor(gono.data$ant)

for (i in 1:nlevels(gono.data$ant)) { #nlevels -> delete duplicated factor level count return (=distinct value)

  for (j in 1:nlevels(gono.data$stam)) {

    ant_value <- levels(gono.data$ant)[i]
    stam_value <- levels(gono.data$stam)[j]
  
    gono.data.sub <- subset(gono.data, ant == ant_value & stam == stam_value)

    mod <- with(gono.data.sub, lm(log.MIC.naive ~ 1))
   
    refmic.newdata[k, "mode.log.MIC"] <- log2(as.numeric(sub("<", "", sub("<=", "" ,sub(">", "", sub(c(">="), "", names(which.max(table(gono.data.sub$MIC)))))))))
    refmic.newdata[k, c("E.log.MIC.naive", "se.log.MIC.naive")] <- c(coef(mod), sqrt(vcov(mod)))
   
    k <- k+1
  }
}



refmic.newdata <- within(refmic.newdata, {
  lower.log.MIC.naive <- E.log.MIC.naive-1.96*se.log.MIC.naive
  upper.log.MIC.naive <- E.log.MIC.naive+1.96*se.log.MIC.naive
})

# winbugs
X.pred <- model.matrix(~ (stam+ant)^2, data = refmic.newdata)

mu <- t(X.pred%*%t(beta))
refmic.newdata <- within(refmic.newdata, {
  E.log.MIC <- colMeans(mu)
  lower.log.MIC <- apply(mu, 2, quantile, 0.025)
  upper.log.MIC <- apply(mu, 2, quantile, 0.975)
})
refmic.newdata <- merge(refmic.newdata, refmic.data.sub, sort = F)

#
# plot
#

# labs compared to mean
pdf(file = "resultaten/figure 3.pdf", width = 7, height = 7)
par(mar = c(4.5, 4.5, 0.5, 0.5), yaxs = "i");
plot.new(); 
plot.window(xlim = c(-2, 2), c(n.lab+1, 0))
abline(h = 1:n.lab, col = 8, lty = 3)
box();
axis(1, at = seq(-2, 2, 1), labels = 2^seq(-2, 2, 1)); 
axis(2, at = 1:n.lab, labels = lab.data$lab);
abline(v = 0, lty = 2)
title(xlab = "Fold difference in MIC compared to the expected mean", ylab = "Laboratory number")
with(lab.data, {
  points(diff.log.MIC, 1:n.lab, pch = 15)
  segments(lower.diff.log.MIC, 1:n.lab, upper.diff.log.MIC, 1:n.lab)
})
dev.off()


# refmics
pdf(file = "resultaten/figure 4.pdf", width = 7, height = 7*sqrt(2))
n <- ncol(X.pred); d <- 0.25
#windows(7, 10)
par(mar = c(4.5, 5.5, 0.5, 0.5), yaxs = "i"); 
plot.new(); 
plot.window(xlim = c(-10, 6), c(n+0.5, 0.5))

abline(h = seq(0.5, n+0.5, 1), col = 8); 
abline(h = c(10.5, 20.5), lwd = 2); 
box()
axis(1, at = seq(-10, 6, 2), labels = signif(2^seq(-10, 6, 2), 3), cex.axis = 0.7)
axis(2, at = 1:n, labels = with(refmic.newdata, levels(interaction(stam, ant, sep = " "))), las = 1, cex.axis = 0.7)
title(xlab = "MIC")


with(refmic.newdata, {
  points(mode.log.MIC, 1:n, pch = 0, cex = 0.7)
  points(E.log.MIC, 1:in-d, pch = 15, cex = 0.7)
  segments(lower.log.MIC, 1:n-d, upper.log.MIC, 1:n-d)
  points(lower.log.MIC.ref, 1:n+d, col = 1, cex = 0.7)
  points(upper.log.MIC.ref, 1:n+d, col = 1, pch = 16, cex = 0.7)
  segments(lower.log.MIC.ref, 1:n+d, upper.log.MIC.ref, 1:n+d, col = 1)
})

dev.off()

####################################
# refmics


# grouped
data_grouped <- data_long %>% filter(!is.na(ant))


for( i in 1:nlevels(data_grouped$ant)){
  
  ant_value <- levels(data_grouped$ant)[i]
  
  pdf(file = "resultaten/figure 4_"+ant_value+"pdf", width = 7, height = 7*sqrt(2))
  # 
  # 
  # n <- ncol(X.pred); d <- 0.25
  # #windows(7, 10)
  # par(mar = c(4.5, 5.5, 0.5, 0.5), yaxs = "i"); plot.new(); plot.window(xlim = c(-10, 6), c(n+0.5, 0.5))
  # abline(h = seq(0.5, n+0.5, 1), col = 8); abline(h = c(10.5, 20.5), lwd = 2); box()
  # axis(1, at = seq(-10, 6, 2), labels = signif(2^seq(-10, 6, 2), 3), cex.axis = 0.7)
  # axis(2, at = 1:i, labels = with(refmic.newdata, levels(interaction(stam, ant, sep = " "))), las = 1, cex.axis = 0.7)
  # title(xlab = "MIC")
  # 
  # 
  # with(refmic.newdata, {
  #   points(mode.log.MIC, 1:n, pch = 0, cex = 0.7)
  #   points(E.log.MIC, 1:n-d, pch = 15, cex = 0.7)
  #   segments(lower.log.MIC, 1:n-d, upper.log.MIC, 1:n-d)
  #   points(lower.log.MIC.ref, 1:n+d, col = 1, cex = 0.7)
  #   points(upper.log.MIC.ref, 1:n+d, col = 1, pch = 16, cex = 0.7)
  #   segments(lower.log.MIC.ref, 1:n+d, upper.log.MIC.ref, 1:n+d, col = 1)
  # })
  # 
  # 
  # dev.off()
  
  
  
}


####################################



pdf(file = "resultaten/figure 4_survival.pdf", width = 7, height = 7*sqrt(2))

glimpse(refmic.newdata)
surv_object <- Surv(time = refmic.newdata$mode.log.MIC, event = refmic.newdata$upper.log.MIC)

f1 <- fitdist(surv_object ~ stam ,"norm")
summary(gofstat(f1))
summary(f1)
plot(f1)



dev.off()


# plotdist(as.numeric(refmic.newdata$mode.log.MIC), histo = TRUE, demp = TRUE) #The descdist function provides classical descriptive statistics (minimum, maximum, median, mean, standard deviation)
# descdist(refmic.newdata$MIC, boot=1000)


# #
# # write tables
# #
# 
# tmp <- within(refmic.newdata, {
#   MIC.unemo <- MIC
#   MIC.unemo.lower <- 2^lower
#   MIC.unemo.upper <- 2^upper
#   MIC.hanjan <- 2^E.log.MIC
#   MIC.hanjan.lower <- 2^lower.log.MIC
#   MIC.hanjan.upper <- 2^upper.log.MIC
# })
# 
# ix <- with(tmp, order(volgstam, volgant))
# 
# 
# tmp <- tmp[ix, c("stam", "ant", "volgstam", "volgant", "MIC.unemo", "MIC.unemo.lower", "MIC.unemo.upper",
#                  "MIC.hanjan", "MIC.hanjan.lower", "MIC.hanjan.upper")]
# write.table(tmp, file = "resultaten/tabel refmic unemo hanjan.txt", quote = F, sep = "\t", row.names = F)
# 
# #
# # calculate mean refmic difference per ant
# #
# 
# calc.reldiff <- function(ix) {
#   x1 <- mu[, ix] # make subset van mu
#   x2 <- NULL; for (i in 1:length(ix)) x2 <- cbind(x2, runif(1000, refmic.newdata[ix[i], "lower"], refmic.newdata[ix[i], "upper"]))
#   dx <- as.vector(x1-x2)
#   cat(ix, "\n")
#   cat(round(2^c(mean(dx), quantile(dx, c(0.025, 0.975))), 2), "\n")
# }
# with(refmic.newdata, calc.reldiff(which(upper!=100 & ant=="cip"))) # niet right cens, cip
# with(refmic.newdata, calc.reldiff(which(upper!=100 & ant=="pen"))) # niet right cens, pen
# with(refmic.newdata, calc.reldiff(which(upper!=100 & ant=="pen" & stam!="strain0909"))) # niet right cens, pen, -strain0909
# with(refmic.newdata, calc.reldiff(which(upper!=100 & ant=="tet"))) # niet right cens, tet
# with(refmic.newdata, calc.reldiff(which(upper!=100))) # niet right cens, alles
# with(refmic.newdata, calc.reldiff(which(upper!=100 & !(stam=="strain0909" & ant=="pen")))) # niet right cens, alles, -strain0909 & -pen
