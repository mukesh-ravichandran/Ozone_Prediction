rm(list = ls())
library(haven)
library(forecast)
library(fma)
library(tseries)
library(expsmooth)
library(lmtest)
library(zoo)
library(ggplot2)
library(lubridate)
library(gtable)
library(dplyr)
library(MLmetrics)
library(stinepack)
library(data.table)
library(imputeTS)
library(gridExtra)

setwd("/Users/mukeshravichandran/OneDrive - North Carolina State University/MSA-Mukesh’s MacBook Pro/Courses/Fall/AA502/Fall 1/Time Series/Project/Code/Ozone_Prediction/Data")

#Theme for graphs
my_theme = theme_minimal()+theme(plot.title = element_text(hjust = 0.5), 
                      panel.background = element_rect(fill='white')) + theme(panel.grid.major = element_blank(),
                                                                          panel.grid.minor = element_blank())

#Reading all the project files in
temp = list.files(pattern="*.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))


#Creating all the functions necessary

#white noise plot
white.noise = function(df) {
  White.LB <- rep(NA, 10)
  for(i in 1:10){
    White.LB[i] <- Box.test(df, lag = i, type = "Lj", fitdf = 1)$p.value
  }
  
  White.LBred <- pmin(White.LB, 0.2)
  barplot(White.LBred, main = "Ljung-Box Test P-values", ylab = "Probabilities", xlab = "Lags", ylim = c(0, 0.2))
  abline(h = 0.01, lty = "dashed", col = "black")
  abline(h = 0.05, lty = "dashed", col = "black")
  return(White.LB)
}

#adftest for checking stationarity

adf= function(series) {

ADF.Pvalues <- rep(NA, 3)
for(i in 0:2){
  ADF.Pvalues[i+1] <- adf.test(series, alternative = "stationary", k = i)$p.value
}

ADF.Pvalues

}

ozone = select(Ozone_Raleigh2.csv,c("Date", "Daily.Max.8.hour.Ozone.Concentration"))
ozone$Date = mdy(ozone$Date)

#Converting in to Date time and imputing missing values
Date <- seq(min(ozone$Date), max(ozone$Date), by = 1) 
Date= as.data.frame(Date)

ozone= full_join(Date,ozone,by="Date", suffix=c("",""))
sum(is.na(ozone$Daily.Max.8.hour.Ozone.Concentration)) #60 values missing

ozone$Daily.Max.8.hour.Ozone.Concentration= na_interpolation(ozone$Daily.Max.8.hour.Ozone.Concentration,option = "spline")
sum(is.na(ozone$Daily.Max.8.hour.Ozone.Concentration)) #no values missing

#converting in to time series
ozone_ts= ts(ozone$Daily.Max.8.hour.Ozone.Concentration,start= 2014,frequency = 365.25)
plot1= autoplot(ozone_ts, color= "blue")+my_theme


training= subset(ozone_ts,end=length(ozone_ts)-(28+14))
validation= subset(ozone_ts,start=length(ozone_ts)-(28+14-1), end=length(ozone_ts)-14)
testing= subset(ozone_ts,start=length(ozone_ts)-(13), end=length(ozone_ts))

autoplot(training)+my_theme

#checking if the split is done properly
length(training)
length(validation)
length(testing)
length(ozone_ts)
sum(ozone_ts == c(training,validation,testing)) == length(ozone_ts)

#####ARIMA Modelling#####

#creating regression variables
# one for easonaluty and one for trend
xreg1<-cbind(fourier(training,K=6),seq(1,length(training)))
colnames(xreg1)<-c('s1','c1','s2','c2','s3','c3','s4','c4','s5','c5','s6','c6','time')


Trial_1_Seasonal<-Arima(training,order=c(0,0,0),seasonal=c(0,0,0),xreg=xreg1)


autoplot(Trial_1_Seasonal$residuals,color="red")+ geom_hline(yintercept=0,linetype="dashed")

adf(Trial_1_Seasonal$residuals) # seems like the residuals are stationary

#lets check for white noise

white.noise(Trial_1_Seasonal$residuals) #--> No white noise, so there is some correlation

Acf(Trial_1_Seasonal$residuals,lag.max = 400)$acf
Pacf(Trial_1_Seasonal$residuals,lag.max = 400)$acf

#no idea what they mean

#lets see what the computer is telling us

auto.arima(Trial_1_Seasonal$residuals)

#it says an AR2 model, lets check that

Trial_2_Seasonal<-Arima(training,order=c(3,0,0),seasonal=c(0,0,0),xreg=xreg1)
autoplot(Trial_2_Seasonal$residuals,color="red")+ geom_hline(yintercept=0,linetype="dashed") #definitely better than the first one

adf(Trial_2_Seasonal$residuals)

white.noise(Trial_2_Seasonal$residuals) # yaay white noise
summary(Trial_2_Seasonal)



autoplot(forecast(Trial_2_Seasonal,xreg = xreg1[109.5:137.5,])) + my_theme + labs(title = "ARIMA Model", y= "Daily 8 Hour Maximum Ozone Concentration", X= "Date")

#Mape = 17.74585
#MAE = 0.005807725
#on the training data

arima.forecast= forecast(Trial_2_Seasonal,xreg = xreg1[109.5:137.5,])[["mean"]]


MAE(arima.forecast,validation) 
#MAPE and MAE of ARIMAX model on the validation == 11.7 and 0.0045


##### ARIMAX Modelling #####


#imputing CO missing values
CO = select(CO_Raleigh.csv,c("Date", "Daily.Max.8.hour.CO.Concentration"))
CO$Date = mdy(CO$Date)

Date <- seq(min(CO$Date), max(CO$Date), by = 1) 
Date= as.data.frame(Date)

CO= full_join(Date,CO,by="Date", suffix=c("",""))
sum(is.na(CO$Daily.Max.8.hour.CO.Concentration))

CO$Daily.Max.8.hour.CO.Concentration= na_interpolation(CO$Daily.Max.8.hour.CO.Concentration,option = "spline")

CO_ts= ts(CO$Daily.Max.8.hour.CO.Concentration, start = 2014, frequency = 365.25)

plot2= autoplot(CO_ts, color= "black")+my_theme


#imputing SO2 missing values
SO2 = select(SO2_Raleigh.csv,c("Date", "Daily.Max.1.hour.SO2.Concentration"))
SO2$Date = mdy(SO2$Date)

Date <- seq(min(SO2$Date), max(SO2$Date), by = 1) 
Date= as.data.frame(Date)

SO2= full_join(Date,SO2,by="Date", suffix=c("",""))
sum(is.na(SO2$Daily.Max.1.hour.SO2.Concentration))

SO2$Daily.Max.1.hour.SO2.Concentration= na_interpolation(SO2$Daily.Max.1.hour.SO2.Concentration,option = "spline")

SO2_ts= ts(SO2$Daily.Max.1.hour.SO2.Concentration, start = 2014, frequency = 365.25)

plot3= autoplot(SO2_ts, color= "green")+my_theme



#imputing NO2 missing values

NO2 = select(NO_Raleigh.csv,c("Date", "Daily.Max.1.hour.NO2.Concentration"))
NO2$Date = mdy(NO2$Date)

Date <- seq(min(NO2$Date), max(NO2$Date), by = 1) 
Date= as.data.frame(Date)

NO2= full_join(Date,NO2,by="Date", suffix=c("",""))
sum(is.na(NO2$Daily.Max.1.hour.NO2.Concentration))

NO2$Daily.Max.1.hour.NO2.Concentration= na_interpolation(NO2$Daily.Max.1.hour.NO2.Concentration,option = "spline")

NO2_ts= ts(NO2$Daily.Max.1.hour.NO2.Concentration, start = 2014, frequency = 365.25)

plot4= autoplot(NO2_ts, color= "red")+my_theme


#imputing TMAX missing values

tmax = select(Raleigh_weather.csv,c("DATE", "TMAX"))
tmax$DATE = mdy(tmax$DATE)

DATE <- seq(min(tmax$DATE), max(tmax$DATE), by = 1) 
DATE= as.data.frame(DATE)

tmax= full_join(DATE,tmax,by="DATE", suffix=c("",""))
sum(is.na(tmax$TMAX))

tmax$TMAX= na_interpolation(tmax$TMAX,option = "spline")

tmax_ts= ts(tmax$TMAX, start = 2014, frequency = 365.25)

plot5= autoplot(tmax_ts, color= "skyblue")+my_theme






grid.arrange(plot1+
               theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank()), plot2+
               theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank()),plot3+
               theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank()),plot4+
               theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank()),plot5, ncol=1)



####Arima to model NO2####
training_NO2= subset(NO2_ts,end=length(NO2_ts)-(28+14))
validation_NO2= subset(NO2_ts,start=length(NO2_ts)-(28+14-1), end=length(NO2_ts)-14)
testing_NO2= subset(ozone_ts,start=length(NO2_ts)-(13), end=length(NO2_ts))

xreg1_NO2<-cbind(fourier(training_NO2,K=6),seq(1,length(training_NO2)))
colnames(xreg1_NO2)<-c('s1','c1','s2','c2','s3','c3','s4','c4','s5','c5','s6','c6','time')


NO2_Seasonal<-Arima(training_NO2,order=c(0,0,0),seasonal=c(0,0,0),xreg=xreg1_NO2)

autoplot(NO2_Seasonal$residuals)

adf(NO2_Seasonal$residuals)
white.noise(NO2_Seasonal$residuals)

Acf(NO2_Seasonal$residuals)$acf
Pacf(NO2_Seasonal$residuals,lag.max = 400)$acf


auto.arima(NO2_Seasonal$residuals)

NO2_Seasonal_2<-Arima(training_NO2,order=c(2,0,0),seasonal=c(0,0,0),xreg=xreg1_NO2)

adf(NO2_Seasonal_2$residuals)
white.noise(NO2_Seasonal_2$residuals)    #yaay white noise





ccf(NO2_Seasonal_2$residuals,Trial_2_Seasonal$residuals,lag.max = 10) 
#there is a spike at lag 0 -- so it is the immediate effect


####Arima to model TMAX####

training_tmax= subset(tmax_ts,end=2301)
validation_tmax= subset(tmax_ts,start=2302, end=2329)
testing_tmax= subset(tmax_ts,start=2330, end=2343)

xreg1_tmax<-cbind(fourier(training_tmax,K=6),seq(1,length(training_tmax)))
colnames(xreg1_tmax)<-c('s1','c1','s2','c2','s3','c3','s4','c4','s5','c5','s6','c6','time')






tmax_seasonal<-Arima(training_tmax,order=c(0,0,0),seasonal=c(0,0,0),xreg=xreg1_tmax)

autoplot(tmax_seasonal$residuals)

adf(tmax_seasonal$residuals)
white.noise(tmax_seasonal$residuals)

Acf(tmax_seasonal$residuals)$acf
Pacf(tmax_seasonal$residuals,lag.max = 400)$acf


auto.arima(tmax_seasonal$residuals)

tmax_seasonal_2<-Arima(training_tmax,order=c(3,0,0),seasonal=c(0,0,0),xreg=xreg1_tmax)

adf(tmax_seasonal_2$residuals)
white.noise(tmax_seasonal_2$residuals)    #yaay white noise





ccf(tmax_seasonal_2$residuals,Trial_2_Seasonal$residuals,lag.max = 10) 






arimax= Arima(training,order=c(0,0,0),seasonal=c(0,0,0),xreg=cbind(xreg1,training_NO2,training_tmax))

adf(arimax$residuals)
white.noise(arimax$residuals)

auto.arima(arimax$residuals)

arimax_2= Arima(training,order=c(1,0,0),seasonal=c(0,0,0),xreg=cbind(xreg1,training_NO2,training_tmax))

summary(arimax_2)

#MAPE 0.005375724
#MAE 16.52266

arimax.forecast= forecast(arimax_2,xreg = cbind(xreg1[109.5:136.5,],validation_NO2,validation_tmax))[["mean"]]


MAPE(arimax.forecast,validation) 


#MAPE 0.086
#MAE 0.00348

#####NEURAL NETWORKS#####

NN.Model <- nnetar(training, p = 3, size = 2)
NN.Forecast <- forecast(NN.Model, h = 200)
plot(NN.Forecast)


Model.four<-Arima(training,order=c(0,0,0),xreg=xreg1)
NN.Model2<-nnetar(Model.four$residuals,p=3,size=2)
NN.Forecast2<-forecast(NN.Model2,h=500)
plot(NN.Forecast2)

Pass.Forecast <- rep(NA, 28)
for(i in 1:365){
  Pass.Forecast[i] <- training[length(training) - 365 + i] + forecast(NN.Model2, h = 28)$mean[i]
}

Pass.Forecast <- ts(Pass.Forecast, start = 2020.29979466119, frequency = 365.24)

plot(training, xlim = c(2014, 2021))
lines(Pass.Forecast, col = "blue")
abline(v = 2020.29979466119, col = "red", lty = "dashed")

NN.forecast= Pass.Forecast[1:28]

MAE(NN.forecast,validation)

#MAPe = 23.12412
#MAe = 0.009

##Can definitely do better with NN - something is off


#ENSEMBLE MODELLING


  
Ensemble1= (arima.forecast + arimax.forecast + NN.forecast)/3

MAPE(Ensemble1, validation)

Ensemble2= (arima.forecast + arimax.forecast)/2
  
MAPE(Ensemble2, validation)


