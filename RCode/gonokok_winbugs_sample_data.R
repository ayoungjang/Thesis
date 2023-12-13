# nolint: commented_code_linter.
# Adjusted code for Han de Neeling data
library(survival)
library(R2WinBUGS)
library(survminer)
library(dplyr)
library(fitdistrplus)
library(ggplot2)
library(truncnorm)

setwd("C:/Users/ayoung/Desktop/Thesis/RCode")
getwd()


min_v <- 1
max_v <- 32
s <- 20 # strain numbers

n_o <- 10 # number of observations
# 10 lab 20 strains

# random number within boundaries
lb <- 0.1 # lower boundary
ub <- 32 # upper boundary

left <- 0.03125
right <- 16
# 1< MIC < 32 increasing MIC
mic_arr <- c( 0.25, 1, 2, 4, 8, 16, 32)
sample.data <- data.frame()



# Function to generate random numbers with specified mean and sd within a range
generate_random_numbers <- function(n, mean_value, sd_value, min_range, max_range,mic_value) {
  # Generate random numbers from a normal distribution with specified mean and sd
  random_numbers <- rnorm(n, mean = mean_value, sd = sd_value)
  
  # Allow a certain proportion of values to go beyond the specified range
  # You can adjust the proportion based on your requirements
  proportion_out_of_range <- 0.1
  num_out_of_range <- round(n * proportion_out_of_range)
  
  # Randomly select indices to go beyond the range
  out_of_range_indices <- sample(1:n, num_out_of_range)
  
  # Generate values beyond the range
  # random_numbers[out_of_range_indices] <- runif(num_out_of_range, min_range, max_range)
  random_numbers[out_of_range_indices] <-rtruncnorm(n = num_out_of_range, a = mic_arr[1], b = mic_arr[length(mic_arr)], mean = mic_value, sd = sd_value)
  
  return(random_numbers)
}




for (x in 1:s) {
  set.seed(x) # reset
  # with chosen mean and SD

  lab <- 1:n_o
  stam <- as.character(x)
  Y_act <- numeric(n_o)
  Y_obs <- numeric(n_o) # can divide log2
  lower <- numeric(n_o)
  upper <- numeric(n_o)
  MIC <- numeric(n_o)

  mic_idx <- sample(1:length(mic_arr), 1)
  mic_value <- mic_arr[mic_idx]
  m <- round(runif(1, min = 1, max = 32), 1)
  sd <- mic_value * 2
  min_idx <- ifelse(mic_idx == 1, 1, mic_idx - 1)
  max_idx <- ifelse(mic_idx == length(mic_arr), length(mic_arr), mic_idx + 1)

  min_value <- mic_arr[min_idx]
  max_value <- mic_arr[max_idx]

  Y_act <- generate_random_numbers(n_o, mic_value, sd, min_value, max_value,mic_value)
  # Y_act <- rtruncnorm(n = n_o, a = min_value, b = max_value, mean = mic_value, sd = sd)
  
  for (i in 1:n_o) {
    Y_act[i] <- max(lb, Y_act[i])

    Y_obs[i] <- 2^ceiling(log2(Y_act[i]))

    lower[i] <- 2^floor(log2(Y_act[i]))

    upper[i] <- 2^ceiling(log2(Y_act[i]))
    MIC[i] <- mic_value
  }
  # mean should be close to 4
  sample.data <- rbind(sample.data, data.frame(lab, stam, MIC, Y_act, Y_obs, lower, upper))
}



####
####


# make MIC column

#
# # reshape into long format
# sample.data <- reshape(data = sample.data,
#                      varying = names(sample.data)[!is.element(names(sample.data), c("Strain"))],
#                      v.names = "MIC", timevar = "stam",
#                      times = names(sample.data)[!is.element(names(sample.data), c("Strain"))],
#                      direction = "long")


# create censored data objects
sample.data <- within(sample.data, {
  # as character
  # MIC <- as.character(MIC)
  # event type: 3 = interval, 2 = left, 0 = right
  event <- rep(3, nrow(sample.data))
  event[Y_act < left] <- 2

  # Set event to 2 for Y_act greater than right threshold
  event[Y_act > right] <- 0

  # log MIC
  log.MIC <- log2(MIC)
  # log MIC for naive mean
  log.MIC.naive <- ifelse(event != 0, log.MIC - 0.5, log.MIC + 0.5)
  # as Surv
  mode.log.MIC <- log2(as.numeric(names(which.max(table(sample.data$MIC)))))
  log.MIC.surv <- Surv(time = ifelse(event == 3, log.MIC - 1, log.MIC), time2 = log.MIC, event = event, type = "interval")
  # lower and upper for WinBUGS
  # lower <- ifelse(event == 3, log.MIC - 1, ifelse(event == 2, -100, log.MIC))
  # upper <- ifelse(event == 3, log.MIC, ifelse(event == 2, log.MIC, 100))
})


sample.data.sub <- subset(sample.data, select = c(stam, upper, lower, log.MIC))
names(sample.data.sub) <- c("stam", "upper.log.MIC.ref", "lower.log.MIC.ref", "log.MIC.ref")


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
X <- model.matrix(~ (stam)^2, data = sample.data)
n <- nrow(X)
n.beta <- ncol(X)
n.lab <- 10
bugs.data <- with(sample.data, list(n = n, lower = lower, upper = upper, X = X, n.beta = n.beta, n.lab = n.lab, lab = as.numeric(lab)))

# bugs inits
bugs.inits <- function() {
  list(y = with(sample.data, runif(n, lower, upper)), b.lab = rnorm(n.lab, 0, 0.1), beta = rnorm(n.beta), sigma = runif(1), sigma.lab = runif(1))
}


# bugs fit
bugs.fit <- bugs(
  model.file = "model.txt", data = bugs.data, inits = bugs.inits,
  parameters.to.save = c("beta", "b.lab", "sigma", "sigma.lab"),
  n.chains = 2, n.iter = 5100, n.burnin = 100, n.thin = 10, debug = F, DIC = F, working.directory = wbwd
)

read.bugsfit <- T

if (read.bugsfit) {
  old.wd <- getwd()
  setwd(wbwd)
  bugs.fit <<- R2WinBUGS:::bugs.sims(
    n.chains = 2, n.iter = 5100, n.burnin = 100, n.thin = 10,
    parameters.to.save = c("beta", "b.lab", "sigma", "sigma.lab"), DIC = F
  )
  class(bugs.fit) <- "bugs"
  setwd(old.wd)
}

attach.bugs(bugs.fit)

colSD <- function(x) { # sd
  sqrt(colSums((x - colMeans(x))^2) / (nrow(x) - 1))
}

sample.data$lab <- as.factor(sample.data$lab)
# labs compared to mean
# confidence value

lab.data <- data.frame(
  lab = levels(sample.data$lab),
  diff.log.MIC = colMeans(b.lab),
  lower.diff.log.MIC = apply(b.lab, 2, quantile, 0.025),
  upper.diff.log.MIC = apply(b.lab, 2, quantile, 0.975)
)

# lab.data<-data.frame(
#   lab=levels(sample.data$lab),
#   diff.log.MIC = colMeans(b.lab),
#   lower.diff.log.MIC = colMeans(b.lab) - 1.96 *  colSD(b.lab)/ sqrt(s),
#   upper.diff.log.MIC= colMeans(b.lab) + 1.96 * colSD(b.lab) / sqrt(s)
# )

sample.data$stam <- as.factor(sample.data$stam)

# make sample new data
sample.newdata <- with(sample.data, expand.grid(stam = levels(stam)))
sample.newdata <- within(sample.newdata, {
  mode.log.MIC <- NA
  E.log.MIC.naive <- NA
  se.log.MIC.naive <- NA
})

for (j in 1:nlevels(sample.data$stam)) {
  stam_value <- levels(sample.data$stam)[j]

  sample.data.new_sub <- subset(sample.data, stam == stam_value)
  mod <- with(sample.data.new_sub, lm(log.MIC.naive ~ 1))

  sample.newdata$mode.log.MIC[j] <- log2(as.numeric(names(which.max(table(sample.data.new_sub$MIC)))))
  sample.newdata[j, c("E.log.MIC.naive", "se.log.MIC.naive")] <- c(coef(mod), sqrt(vcov(mod)))
}

sample.newdata <- within(sample.newdata, {
  lower.log.MIC.naive <- E.log.MIC.naive - 1.96 * se.log.MIC.naive
  upper.log.MIC.naive <- E.log.MIC.naive + 1.96 * se.log.MIC.naive
})

X.samplepred <- model.matrix(~ (stam)^2, data = sample.newdata)

mu <- X.samplepred %*% t(beta)

sample.newdata <- within(sample.newdata, {
  E.log.MIC <- colMeans(mu)
  lower.log.MIC <- apply(mu, 2, quantile, 0.025)
  upper.log.MIC <- apply(mu, 2, quantile, 0.975)
})

sample.newdata <- merge(sample.newdata, sample.data.sub, sort = F)

# labs compared to mean
pdf(file = "resultaten/figure_sample_3.pdf", width = 7, height = 7)
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

# pdf(file = "resultaten/figure_sample_4.pdf", width = 7, height = 7*sqrt(2))
# n <- ncol(X.samplepred); d <- 0.25
#
# par(mar = c(4.5, 5.5, 0.5, 0.5), yaxs = "i");
# plot.new();
# plot.window(xlim = c(-10, 6), c(n+0.5, 0.5))
# abline(h = seq(0.5, n+0.5, 1), col = 8);
# abline(h = c(10.5, 20.5), lwd = 2); box()
# axis(1, at = seq(-10, 6, 2), labels = signif(2^seq(-10, 6, 2), 3), cex.axis = 0.7)
# axis(2, at = 1:n, labels = with(sample.newdata, levels(interaction(stam, sep = " "))), las = 1, cex.axis = 0.7)
#
#
# title(xlab = "MIC")
# with(sample.newdata, {
#   points(mode.log.MIC, 1:n, pch = 0, cex = 0.7)
#   points(E.log.MIC, 1:n-d, pch = 15, cex = 0.7)
#   segments(lower.log.MIC, 1:n-d, upper.log.MIC, 1:n-d)
#   points(lower.log.MIC.ref, 1:n+d, col = 1, cex = 0.7)
#   points(upper.log.MIC.ref, 1:n+d, col = 1, pch = 16, cex = 0.7)
#   segments(lower.log.MIC.ref, 1:n+d, upper.log.MIC.ref, 1:n+d, col = 1)
# })
#

# sample.newdata$stam <- factor(sample.newdata$stam, levels = unique(sample.newdata$stam[order(as.numeric( sample.newdata$stam))]))
# pdf_path = "resultaten/figure_sample_4.pdf"
# pdf(file = pdf_path,width = 7, height = 7*sqrt(2))
# n <- ncol(X.samplepred);
# d <- 0.25
#
# par(mar = c(4.5, 5.5, 0.5, 0.5), yaxs = "i");
# plot.new();
# plot.window(xlim = c(-10, 6), c(n+0.5, 0.5))
# abline(h = seq(0.5, n+0.5, 1), col = 8);
# abline(h = c(10.5, 20.5), lwd = 2); box()
# axis(1, at = seq(-10, 6, 2), labels = signif(2^seq(-10, 6, 2), 3), cex.axis = 0.7)
# axis(2, at = 1:n, labels = with(sample.newdata, levels(interaction(stam, sep = " "))), las = 1, cex.axis = 0.7)
# title(xlab = "MIC")
#
# with(sample.newdata, {
#   points(mode.log.MIC, 1:length(mode.log.MIC), pch = 0, cex = 0.7)
#   points(E.log.MIC, 1:length(E.log.MIC)-d, pch = 15, cex = 0.7)
#   segments(lower.log.MIC, 1:length(lower.log.MIC)-d, upper.log.MIC, 1:length(upper.log.MIC)-d)
#   points(lower.log.MIC.ref, 1:length(lower.log.MIC.ref)+d, col = 1, cex = 0.7)
#   points(upper.log.MIC.ref, 1:length(upper.log.MIC.ref)+d, col = 1, pch = 16, cex = 0.7)
#   segments(lower.log.MIC.ref, 1:length(lower.log.MIC.ref)+d, upper.log.MIC.ref, 1:length(upper.log.MIC.ref)+d, col = 1)
# })
#
# dev.off()


sample.newdata$stam <- factor(sample.newdata$stam, levels = unique(sample.newdata$stam[order(as.numeric(sample.newdata$stam))]))
pdf_path <- "resultaten/figure_sample_4_plot.pdf"
pdf(file = pdf_path, width = 7, height = 7 * sqrt(2))
n <- ncol(X.samplepred)
d <- 0.25

par(mar = c(4.5, 5.5, 0.5, 0.5), yaxs = "i")
plot.new()
plot.window(xlim = c(-10, 6), c(n + 0.5, 0.5))
abline(h = seq(0.5, n + 0.5, 1), col = 8)
abline(h = c(10.5, 20.5), lwd = 2)
box()
axis(1, at = seq(-10, 6, 2), labels = signif(2^seq(-10, 6, 2), 3), cex.axis = 0.7)
axis(2, at = 1:n, labels = with(sample.newdata, levels(interaction(stam, sep = " "))), las = 1, cex.axis = 0.7)
title(xlab = "MIC")

with(sample.newdata, {
  points(mode.log.MIC, 1:length(mode.log.MIC), pch = 0, cex = 0.7)
  points(E.log.MIC, 1:length(E.log.MIC) - d, pch = 15, cex = 0.7)
  segments(lower.log.MIC, 1:length(lower.log.MIC) - d, upper.log.MIC, 1:length(upper.log.MIC) - d)
  points(lower.log.MIC.ref, 1:length(lower.log.MIC.ref) + d, col = 1, cex = 0.7)
  points(upper.log.MIC.ref, 1:length(upper.log.MIC.ref) + d, col = 1, pch = 16, cex = 0.7)
  segments(lower.log.MIC.ref, 1:length(lower.log.MIC.ref) + d, upper.log.MIC.ref, 1:length(upper.log.MIC.ref) + d, col = 1)
})

dev.off()
sample.newdata$stam <- factor(sample.newdata$stam, levels = unique(sample.newdata$stam[order(as.numeric(sample.newdata$stam))]))

d <- 0.25

pdf_path <- "resultaten/figure_sample_4.pdf"

pdf(file = pdf_path, width = 7, height = 7 * sqrt(2))
# Set up the ggplot object
ggplot(sample.newdata, aes(x = mode.log.MIC, y = as.numeric(as.character(stam)))) +
  geom_point(pch = 0, size = 1.2) +
  geom_point(aes(x = E.log.MIC, y = as.numeric(as.character(stam)) + 0.2), col = "black", pch = 15, size = 2) +
  geom_segment(aes(x = lower.log.MIC, xend = upper.log.MIC, y = as.numeric(as.character(stam)), yend = as.numeric(as.character(stam))), col = "blue") +
  geom_point(aes(x = lower.log.MIC.ref, y = as.numeric(as.character(stam)) + 0.3), col = "black", size = 0.7) +
  geom_point(aes(x = upper.log.MIC.ref, y = as.numeric(as.character(stam))), col = "black", pch = 16, size = 0.7) +
  geom_segment(aes(x = lower.log.MIC.ref, xend = upper.log.MIC.ref, y = as.numeric(as.character(stam)), yend = as.numeric(as.character(stam))), col = "black") +
  scale_x_continuous(breaks = seq(-10, 6, 2), labels = signif(2^seq(-10, 6, 2), 3), expand = c(0, 0), limits = c(-10, 6)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "MIC") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    panel.grid.major.y = element_line(color = "gray", linetype = "dashed"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "none",
    plot.margin = margin(1, 1, 1, 1, "cm") # Adjust the margin values as needed
  )
dev.off()
system(paste("open", pdf_path))
