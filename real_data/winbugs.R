
library(R2WinBUGS)



execute_winbugs <- function(data,name){
  data$Strain_no <- as.factor(data$Strain_no)
  data$Species <- as.factor(data$Species)
  
  data$Strain_no<-drop.levels(data$Strain_no)
  head(data$Strain_no)
  
  
  data$Species<-drop.levels(data$Species)
  head(data$Species)
  
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
  X <- model.matrix(~ (Strain_no)^2, data = data)
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
  
  # naive mean and mode
  # data.newdata <- with(data, expand.grid(Strain_no = levels(Strain_no)))
  
  data.newdata <- with(data, expand.grid(Strain_no = levels(Strain_no)))
  
  data.sub <- within(within(data.newdata, {mode.log.MIC <- NA; upper.log.MIC.ref <- NA;}))
  data.newdata <- within(data.newdata, {lower.log.MIC.ref <- NA; E.log.MIC.naive <- NA; se.log.MIC.naive <- NA;})
  
  k<-1
  for (i in 1:nlevels(data$Strain_no)) {
    strain_value <- levels(data$Strain_no)[i]
  
    data.data.sub <- subset(data, Strain_no == strain_value)
    
    
    mod <- with(data.data.sub, lm(log.MIC.naive ~ 1))
    
    max_indices <- which.max(table(data.data.sub$MIC.num))
    max_value <- as.numeric(names(max_indices))
    
    subset_rows <- data.data.sub$MIC.num == max_value
    
    lower_value <- min(data.data.sub$lower.log.MIC.ref[subset_rows])
    upper_value <- max(data.data.sub$upper.log.MIC.ref[subset_rows])
    
    
    # max_indices <- which(table(data.data.sub$MIC.num) == max(table(data.data.sub$MIC.num)))
    # max_values <- as.numeric(names(max_indices))
    # 
    # max_value <- max(max_values)
    # print(max_value)
    
    data.sub[k,"lower.log.MIC.ref"]<-lower_value
    data.sub[k,"upper.log.MIC.ref"]<-upper_value
      
    data.newdata[k, "mode.log.MIC"] <- log2(max_values[1])
    data.newdata[k, c("E.log.MIC.naive", "se.log.MIC.naive")] <- c(coef(mod), sqrt(vcov(mod)))
    
    k <- k+1
    # }
  }
  # }
  # 
  # debug(test)
  # test();
  
  data.newdata <- within(data.newdata, {
    lower.log.MIC.naive <- E.log.MIC.naive - 1.96 * se.log.MIC.naive
    upper.log.MIC.naive <- E.log.MIC.naive + 1.96 * se.log.MIC.naive
  })
  
  X.samplepred <- model.matrix(~(Strain_no)^2, data = data.newdata)
  
  
  mu_sample <- t(X.samplepred%*%t(beta))
  
  
  
  
  data.newdata <- within(data.newdata, {
    E.log.MIC <- colMeans(mu_sample)
    lower.log.MIC <- apply(mu_sample, 2, quantile, 0.025)
    upper.log.MIC <- apply(mu_sample, 2, quantile, 0.975)
  })
  
  data.newdata <- merge(data.newdata, data.sub, sort = F)
  draw_plot(data, name)
}
