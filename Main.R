#In this code snippet, we are plotting the in-sample and out-of-sample returns with the measure (e.g., VIX, RQ10, etc.) as a line. 
#We also estimate the parameters for GARCHx and CAViaRx models and forecast the in-sample and out-of-sample returns.
#For GARCHx set-up we forecast VaR using the estimated parameters.

#installing all the libraries
library(DEoptim)
library(ggplot2)
library(tseries)
library(knitr)
library(quantmod)
library(tidyverse)
library(grf) 
library(lmtest)

#input data and transformation
IBM_ret <- read.table("IBM_ret.csv", header = TRUE, sep = ",")
IBM_returns<-as.numeric(IBM_ret$IBM_Returns)
vix<-as.numeric(IBM_ret$VIX)
RV<-as.numeric(IBM_ret$RV)
RQ_05<-as.numeric(IBM_ret$RQ_05)
RQ_10<-as.numeric(IBM_ret$RQ_10)

#initial parameters for base models and x-specification model
init_par <- c( 0.1, 0.8, 0.1) #base model
init_params <- c( 0.1, 0.8, 0.1, 0.1) #x-specification model

split_index <- floor(0.75 * length(IBM_returns))
in_sample_returns <- IBM_returns[1:split_index]
out_of_sample_returns <- IBM_returns[(split_index + 1):length(IBM_returns)]

in_sample_vix<-vix[1:split_index]
out_of_sample_vix<-vix[(split_index + 1):length(vix)]

in_sample_rv<-RV[1:split_index]
out_of_sample_rv<-RV[(split_index + 1):length(RV)]

in_sample_RQ_10<-RQ_10[1:split_index]
out_of_sample_RQ_10<-RQ_10[(split_index + 1):length(RQ_10)]

in_sample_RQ_05<-RQ_05[1:split_index]
out_of_sample_RQ_05<-RQ_05[(split_index + 1):length(RQ_05)]

source("Model.R")
#Engel's CaViaR estimation: Symmetric absolute value model
#ft(β) = β1 + β2 ft−1(β)+ β3|yt−1|
alpha=0.01 #significance level which can be chosen manually.
#In optim function, we are using L-BFGS-B method for optimization with lower and upper bounds for the parameters. Since it is a constrained optimization problem, we need to provide specific lower and upper bounds for the parameters.
#Base model
cav<- optim(fn=parestim_base, returns = in_sample_returns, alpha= alpha, par = init_par, method = "L-BFGS-B",lower = c(-Inf, 0, -Inf),
             upper = c(0, 1, 0))
par_base<-cav$par
names(par_base)<-c("beta1","beta2","beta3")
par_base<-c(cav$par,NA)
par_base
#Caviarx model  with exogeneous variables
cavxvix<- optim(par=init_params, fn= parestim_cavx, returns = in_sample_returns, specx = in_sample_vix, alpha = alpha, method = "L-BFGS-B",lower = c(-Inf, 0, -Inf, -Inf),
             upper = c(0, 1, 0, 0) )
names(cavxvix$par)<-c("beta1","beta2","beta3","beta4")
par_vix<-cavxvix$par

cavxRV<- optim(par=init_params, fn= parestim_cavx, returns = in_sample_returns, specx = in_sample_rv, alpha = alpha, method = "L-BFGS-B",lower = c(-Inf, 0, -Inf, -Inf),
             upper = c(0, 1, 0, 0))
names(cavxRV$par)<-c("beta1","beta2","beta3","beta4")
par_rv<-cavxRV$par

cavxrq10<- optim(par=init_params, fn= parestim_cavx, returns = in_sample_returns, specx = in_sample_RQ_10, alpha = alpha, method = "L-BFGS-B",lower = c(-Inf, 0, -Inf, -Inf),
             upper = c(0, 1, 0, 0) )
names(cavxrq10$par)<-c("beta1","beta2","beta3","beta4")
par_rq10<-cavxrq10$par


cavxrq05<- optim(par=init_params, fn= parestim_cavx, returns = in_sample_returns, specx = in_sample_RQ_05, alpha = alpha, method = "L-BFGS-B",lower = c(-Inf, 0, -Inf, -Inf),
                 upper = c(0, 1, 0, 0) )
names(cavxrq05$par)<-c("beta1","beta2","beta3","beta4")
par_rq05<-cavxrq05$par

#table for all cAViaRx parameters. 
Caviarparam_table<-data.frame("Base"=par_base, "VIX"=par_vix, "RV"=par_rv, "RQ_10"=par_rq10, "RQ_05"=par_rq05)
rownames(Caviarparam_table)<-c("beta1","beta2","beta3","beta4")
Caviarparam_table

#next we forecast the in-sample and out-of-sample returns using the estimated parameters.
#In-sample and out-of-sample forecasts
insamcav<-Forecastcav_base(cav$par, in_sample_returns, alpha)

outsamcav<-OutForecastcav_base(cav$par, out_of_sample_returns, insam = insamcav, alpha)

insamcavvix <- Forecastcav_x(cavxvix$par, in_sample_returns, specx = in_sample_vix, alpha)

outsamcavvix <- OutForecastcav_x(cavxvix$par, out_of_sample_returns, insam= insamcavvix, specx= out_of_sample_vix, alpha)

insamcavRV <- Forecastcav_x(cavxRV$par, in_sample_returns, specx= in_sample_rv, alpha)

outsamcavRV <- OutForecastcav_x(cavxRV$par, out_of_sample_returns, insam= insamcavRV, specx= out_of_sample_rv, alpha)

insamcavRQ_10 <-Forecastcav_x(cavxrq10$par, in_sample_returns, specx= in_sample_RQ_10, alpha)

outsamcavRQ_10 <- OutForecastcav_x(cavxrq10$par, out_of_sample_returns, insam= insamcavRQ_10, specx= out_of_sample_RQ_10, alpha)

insamcavRQ_05 <- Forecastcav_x(cavxrq05$par, in_sample_returns, specx= in_sample_RQ_05, alpha)

outsamcavRQ_05 <- OutForecastcav_x(cavxrq05$par, out_of_sample_returns, insam= insamcavRQ_05, specx= out_of_sample_RQ_05, alpha)


#########################################################################################################################################3

#GARCH
garch <- optim(par = init_par, mu=0, fn = garch_loglik, returns = in_sample_returns, method = "BFGS")
parbase <- garch$par
parbase<-c(parbase,NA)
parbase

cat("Log-Likelihood Value: ", -garch$value, "\n")
#GARCHX
init_params <- c(0.1, 0.8, 0.1, 0.1)
garxvix<-optim(par = init_params, mu=0,
               fn = garchx_loglik,
               gr = NULL,
               returns = in_sample_returns,
               specx = in_sample_vix,
               method = "L-BFGS-B",
               lower = c(0, 0,0,0),
               upper = c(Inf, 1, 1, Inf),
               hessian=TRUE)

parvix<-garxvix$par
cat("Log-Likelihood Value: ", -garxvix$value, "\n")
garxrq10<-optim(par = init_params, mu=0,
                fn = garchx_loglik,
                gr = NULL,
                returns = in_sample_returns,
                specx = in_sample_RQ_10,
                method = "L-BFGS-B",
                lower = c(0,0,0,0),
                upper = c(Inf, 1, 1, Inf),
                hessian=TRUE)

names(garxrq10$par)<-c("omega","alpha","beta","gamma")
parrq10<-garxrq10$par
cat("Log-Likelihood Value: ", -garxrq10$value, "\n")
garxrq05<-optim(par = init_params, mu=0,
                fn = garchx_loglik,
                gr = NULL,
                returns = in_sample_returns,
                spec = in_sample_RQ_05,
                method = "L-BFGS-B",
                lower = c(0, 0, 0, 0),
                upper = c(Inf, 1, 1, Inf), hessian=TRUE)

names(garxrq05$par)<-c("omega","alpha","beta","gamma")
parrq05<-garxrq05$par
cat("Log-Likelihood Value: ", -garxrq05$value, "\n")

garxRV<-optim(par = init_params, mu=0,
              fn = garchx_loglik,
              gr = NULL,
              returns = in_sample_returns,
              specx = in_sample_rv,
              method = "L-BFGS-B",
              lower = c(0, 0, 0, 0),
              upper = c(Inf, 1, 1, Inf), hessian=TRUE)

names(garxRV$par)<-c("omega","alpha","beta","gamma")
parrv<-garxRV$par
cat("Log-Likelihood Value: ", -garxRV$value, "\n")
#table 
GARCHparam_table<-data.frame("Base"=parbase, "VIX"=parvix, "RV"=parrv, "RQ_10"=parrq10, "RQ_05"=parrq05)
GARCHparam_table

#loglikihood values
loglik<-data.frame("Base"=-garch$value, "VIX"=-garxvix$value, "RV"=-garxRV$value, "RQ_10"=-garxrq10$value, "RQ_05"=-garxrq05$value)
loglik

#forecasting
#insample
in_sample_variances <- ForecastGarch(garch$par, in_sample_returns)

in_sample_variancesx <- ForecastGarchx(garxvix$par, in_sample_returns, in_sample_vix)
in_sample_variancesrq10 <- ForecastGarchx(garxrq10$par, in_sample_returns, in_sample_RQ_10)
in_sample_variancesrq05 <- ForecastGarchx(garxrq05$par, in_sample_returns, in_sample_RQ_05)
in_sample_variancesRV <- ForecastGarchx(garxRV$par, in_sample_returns, in_sample_rv)

#out of sample
n_ahead <- 10 #forecast horizon which is commonly practiced in literature
outsamvariances <- OutForecastGarch(garch$par, out_of_sample_returns, n_ahead) 

outsamvariancesx <- OutForecastGarchx(garxvix$par, out_of_sample_returns, out_of_sample_vix, n_ahead)
outsamvariancesrq10 <- OutForecastGarchx(garxrq10$par, out_of_sample_returns, out_of_sample_RQ_10, n_ahead)
outsamvariancesrq05 <- OutForecastGarchx(garxrq05$par, out_of_sample_returns, out_of_sample_RQ_05, n_ahead)
outsamvariancesRV <- OutForecastGarchx(garxRV$par, out_of_sample_returns, out_of_sample_rv, n_ahead)

#next we want to calculate VaR with GARCH
#we find volatility by taking square root of the variance
volgarch<-calculate_volatility(outsamvariances)
volgarx<-calculate_volatility(outsamvariancesx)
volgarxrq10<-calculate_volatility(outsamvariancesrq10)
volgarxrq05<-calculate_volatility(outsamvariancesrq05)
volgarxRV<-calculate_volatility(outsamvariancesRV)

#VaR calculation
var_garch <- calculate_var(volgarch, alpha)
var_garx <- calculate_var(volgarx, alpha)
var_garxrq10 <- calculate_var(volgarxrq10, alpha)
var_garxrq05 <- calculate_var(volgarxrq05, alpha)
var_garxRV <- calculate_var(volgarxRV, alpha)
#remove last values for plotting and comparison
var_garch<-var_garch[1:(length(var_garch)-10)]
var_garx<-var_garx[1:(length(var_garx)-10)]
var_garxrq10<-var_garxrq10[1:(length(var_garxrq10)-10)]
var_garxrq05<-var_garxrq05[1:(length(var_garxrq05)-10)]
var_garxRV<-var_garxRV[1:(length(var_garxRV)-10)]


#Next we plot the in-sample and out-of-sample returns with the measure in Plot.R

