setwd('c:/users/elad/Documents')
LabData <- read.csv('LabData.csv')

# define parameters
rHz2PARi = 0.05311
rHz2PARo = 0.01488
BgLi = 0.099
BgLo = 0.167
rLi_0 = 14576 - BgLi
rLo_0 = 3028 - BgLo
RedPARTransCalib = (rLo_0*rHz2PARo) / (rLi_0 * rHz2PARi)
rEpsa = 0.07523
rEpsb = 0.07815
rEpsc = 0.0
width = 4.0

# load data with additional parameters 
LabData<-data.frame(LabData, RedPARTrans=(LabData$rLoref*rHz2PARo)/(LabData$rLiref*rHz2PARi), Chlb_=LabData$Chlb/LabData$Chla, Car_=LabData$Car/LabData$Chla, DW_=LabData$DW/LabData$Chla)

# remove outliers
LabData <- LabData[-c(5,11),]

# plot
# install.packages("ggplot2")
library(ggplot2)
# ggplot(LabData, aes(x=RedPARTrans, y=Chla)) + geom_point()

Chlb_avg = mean(LabData$Chlb_)
Car_avg = mean(LabData$Car_)
DW_avg = mean(LabData$DW_)

SAC_red = rEpsa+rEpsb*Chlb_avg+rEpsc*Car_avg*0.0006*DW_avg
SAC_red

# create smaller dataset - all constants are represented by Ka and kb
dataSet <- data.frame(y  = LabData$Chla,
                      x  = LabData$RedPARTrans,
                      ka = 1/RedPARTransCalib,
                      kb = (rEpsa+rEpsb*Chlb_avg+rEpsc*Car_avg)*DW_avg*width)

dataSet  <- dataSet[order(dataSet$x),]
# # reindex rows
row.names(dataSet) <- 1:nrow(dataSet)
# # splot according to heuristic
dataSet1 <- dataSet[c(15,14,13,11,10,9,8,6,4),]
dataSet2 <- dataSet[-c(15,14,13,11,10,9,8,6,4),]
# # plot 2 datasets
ggplot(dataSet1, mapping=aes(x=x, y=y,colour = "accepted observations")) + geom_point(shape=3) +
  geom_point(data = dataSet2, mapping=aes(x=x, y=y, colour = "outliers"),shape=4) +
  scale_colour_manual("", breaks = c("accepted observations", "outliers"), values = c("orange", "black"),guide = guide_legend(override.aes = list(
    linetype = c("blank", "blank"), shape=c(3,4)))) +
  labs(x="Red PAR transmittance",y="Chloropyll a [mg/L]")

# take bigger dataset
dataSet <- dataSet1
# 
# # devide to 2 datasets
# split_datasets = TRUE
# if (split_datasets){
#   1
#   # dataSet
#   # # notice that the data quality sucks! - fix it!
#   
# }



# dataSet <- dataSet[-c(1,2,4,5,7,12),]
ggplot(dataSet, aes(x=x, y=y)) + geom_point()

#################### optimize k(lambda) #################### 

# function to calculate SSE used by the optimization function
sse = function(beta,y,x,ka,kb) {
  y_pred = (-log10(x*ka)/(kb*beta))^0.5
  errors = y - y_pred
  squared_errors = errors^2
  return(sum(squared_errors))
} 

# run optimization to minimize SSE using beta (Kd parameter)
optimization <- optimise(function(x) sse(x, dataSet$y,dataSet$x,dataSet$ka,dataSet$kb),interval = c(0.0000001,0.1), tol = 0.0000001,maximum=FALSE)
optimization
Kd_opt = optimization$minimum

# calculate results
SST <- sum((dataSet$y - mean(dataSet$y))^2)
SSE <- sse(optimization$minimum, dataSet$y,dataSet$x,dataSet$ka,dataSet$kb)
R_2 <- 1-(SSE/SST)
R_2

# manually check and plot optimum (redundant)
# betas      <- seq(0,0.01,0.00001)
# SSEs       <- sapply(betas, function(x) sse(x, dataSet$y,dataSet$x,dataSet$ka,dataSet$kb))
# betas_SSEs <- data.frame(beta = betas, sse = SSEs)
# library('ggplot2')
# ggplot(data = betas_SSEs, aes(x=beta,y=sse)) + geom_point()

# function that predicts y and return dataframe with y, y_pred and sse   
predict_ <- function(beta, x, y, ka, kb) {
  y_pred = (-log10(x*ka)/(kb*beta))^0.5
  errors = y - y_pred
  return (data.frame(x = x, y = y, y_pred = y_pred, error = errors))
}
prediction_data <- predict_(Kd_opt, dataSet$x, dataSet$y, dataSet$ka, dataSet$kb )

# plot results
ggplot(prediction_data, mapping=aes(x=x, y=y,colour = "y")) + geom_point(shape=3) +
  geom_line(data = prediction_data, mapping=aes(x=x, y=y_pred, colour = "y_pred"),size = 1) + 
  scale_colour_manual("", breaks = c("y", "y_pred"), values = c("red", "green"),guide = guide_legend(override.aes = list(
    linetype = c("blank", "solid"), shape=c(3,NA)))) + 
  # ggtitle("Model 1: Y vs. Y prediction") +
  labs(x="Red PAR transmittance",y="Chloropyll a [mg/L]") 


#################### optimize both pow and k(lambda) ####################

sse1 = function(beta,y,x,ka,kb,pow) {
  y_pred = (-log10(x*ka)/(kb*beta))^pow
  errors = y - y_pred
  squared_errors = errors^2
  return(sum(squared_errors))
} 

sse2 = function(pow,y,x,ka,kb) {
  optimization <- optimise(function(x) sse1(x, dataSet$y,dataSet$x,dataSet$ka,dataSet$kb, pow),interval = c(0.0000001,0.1), tol = 0.0000001,maximum=FALSE)
  beta         <- optimization$minimum
  y_pred = (-log10(x*ka)/(kb*beta))^pow
  errors = y - y_pred
  squared_errors = errors^2
  return(sum(squared_errors))
}

### try optimise using NLS ###

# optimize
optimization2 <- optimise(function(x) sse2(x, dataSet$y,dataSet$x,dataSet$ka,dataSet$kb),interval = c(-2,2), tol = 0.1,maximum=FALSE)
optimization2
pow_opt <- optimization2$minimum
# find matching kd_opt
optimization1 <- optimise(function(x) sse1(x, dataSet$y,dataSet$x,dataSet$ka,dataSet$kb, pow_opt),interval = c(0.0000001,0.1), tol = 0.0000001,maximum=FALSE)
optimization1
kd_opt2 = optimization1$minimum

SST <- sum((dataSet$y - mean(dataSet$y))^2)
SSE <- sse1(kd_opt2, dataSet$y,dataSet$x,dataSet$ka,dataSet$kb, pow_opt)
R_2 <- 1-(SSE/SST)
R_2

# function that predicts y and return dataframe with y, y_pred and sse   
predict2_ <- function(beta, x, y, ka, kb,pow) {
  y_pred = (-log10(x*ka)/(kb*beta))^pow
  errors = y - y_pred
  return (data.frame(x = x, y = y, y_pred = y_pred, error = errors))
}
prediction_data <- predict2_(kd_opt2, dataSet$x, dataSet$y, dataSet$ka, dataSet$kb, pow_opt)

# plot results
ggplot(prediction_data, mapping=aes(x=x, y=y,colour = "y")) + geom_point(shape=3) +
  geom_line(data = prediction_data, mapping=aes(x=x, y=y_pred, colour = "y_pred"),size=1) + 
  scale_colour_manual("", breaks = c("y", "y_pred"), values = c("red", "green"),guide = guide_legend(override.aes = list(
    linetype = c("blank", "solid"), shape=c(3,NA)))) + 
  # ggtitle("Model 1: Y vs. Y prediction - optimize K and pow (partial data set)") +
  labs(x="Red PAR transmittance",y="Chloropyll a [mg/L]") 


################## predict using a model of the form y_hat = 10^-y ####################

model <- lm(-log10(y) ~ x ,data=dataSet)
summary(model)
prediction_data <- data.frame(x=dataSet$x, y=dataSet$y, y_pred=10^-predict(model, data.frame(x=dataSet$x)), error = dataSet$y-10^-predict(model, data.frame(x=dataSet$x)))

SST <- sum((dataSet$y - mean(dataSet$y))^2)
SSE <- sum((prediction_data$y-prediction_data$y_pred)^2)
R_2 <- 1-(SSE/SST)
R_2

# plot results
ggplot(prediction_data, mapping=aes(x=x, y=y,colour = "y")) + geom_point(shape=3) +
  geom_line(data = prediction_data, mapping=aes(x=x, y=y_pred, colour = "y_pred"),size=1) + 
  scale_colour_manual("", breaks = c("y", "y_pred"), values = c("red", "green"),guide = guide_legend(override.aes = list(
    linetype = c("blank", "solid"), shape=c(3,NA)))) + 
  # ggtitle("Model 1: Y vs. Y prediction - optimize K and pow") +
  labs(x="Red PAR transmittance",y="Chloropyll a [mg/L]") 



# plot original data
# plot: add linear regression line
# LabData$predicted <- model$fitted.values
# LabData <- data.frame(LabData[1:12], lapply(LabData[13], function(x) 10^(-x)))
# plot(-log10(Chla) ~ RedPARTrans , data = LabData)
# points(-log10(LabData$predicted) ~ LabData$RedPARTrans, col='red')
