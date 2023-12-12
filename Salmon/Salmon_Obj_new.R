# Martin fish allergy data. 

library(fitdistrplus)
library(readxl)
library(flexsurv)

#Plotting raw data. Salmon objective 

setwd("C:/Users/ayoung/Desktop/Thesis/Salmon")
getwd()

data <- read_excel("Martin_23.09.16.xlsx") #read data from the file


  for(x in 1:6){
    
    st <- 2 * x
    en <- 2 * x + 1
    
    fish_temp <- data[, st:en]
    fish_temp <- na.omit(fish_temp)
    
    #get_fish_name
    col_names <- colnames(fish_temp)
    filename <- sapply(col_names[1][1],function(name){
      letter <- unlist(strsplit(name,"_"))[1]
      letter2 <-unlist(strsplit(name,"_"))[2]
      return (paste(letter,"_",letter2))
    }) 
    
    fish<-data.frame()
    
    # Check if fish_temp has any rows before adding it to fish
    if (nrow(fish_temp) > 0) {
      colnames(fish_temp) <- c("left", "right")
      fish <- rbind(fish, fish_temp)
    }
    
    
    
    fish <-  data.frame(left = fish[, 1], right = fish[, 2]) #make data frame
    
    fish.norm <-fitdistcens(fish,"norm")
    
    #######save file#######
    pdf(file = paste("result/",filename,"_distribution.pdf"), width = 7, height = 7*sqrt(2))
    plotdistcens(fish,Turnbull=FALSE)
    cdfcompcens(list(fish.norm),main=paste(filename," Objective reactions"),xlab="Mg protein",ylab="Cumulative proportion of responses")
    dev.off()
    #######save file#######
    fish_obj_log<-data.frame()
    
    fish_obj_log<-log10(fish)
    fish_obj_log <-  data.frame(left = fish_obj_log[, 1], right = fish_obj_log[, 2]) #make data frame
    
    for (col in colnames(fish_obj_log)) {
      fish_obj_log[fish_obj_log == -Inf] <- fish_obj_log[fish_obj_log == -Inf & col == "left", "right"]
    }
    
    fish_obj_log.norm<-fitdistcens(fish_obj_log,"norm")
    
    
    #######save file#######
    pdf(file = paste("result/",filename,"_log_distribution.pdf"), width = 7, height = 7*sqrt(2))
    plotdistcens(fish_obj_log,Turnbull=FALSE)
    
    #Plotting log-normal data. Salmon objective
    cdfcompcens(list(fish.norm),main=paste(filename," Objective reactions"),xlab="Mg(log10) protein",ylab="Cumulative proportion of responses")
    dev.off()
    #######save file#######
    
    
    pdf(file = paste("result/",filename,"_subjective log-normal.pdf"), width = 7, height = 7*sqrt(2))
    #Bootstrap. subjective log-normal
    subb1<-bootdistcens(fish.norm,niter=1001)
    
    summary(subb1)
    plot(subb1)
    quantile(subb1)
    CIcdfplot(subb1,main=paste(filename,"Subjective reactions"),xlab="Mg(log10) protein",ylab="Cumulative proportion of responses",CI.output="probability",CI.type="two.sided",CI.level=0.95,CI.col="black",CI.lty=2,CI.fill=NULL,CI.only=FALSE)
    dev.off()
  }


  
 
