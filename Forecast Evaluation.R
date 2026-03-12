#Forecast Evaluation for CAViaRx models

#Loading the required libraries
library(grf) 
library(lmtest)
library(segMGarch)
library(dplyr)
library(sandwich)
#Dynamic Quantile test by Engel and Manganelli (2004)
#This function is used to test the null hypothesis that the forecast is equal to the actual value
#We use it manually for different quantiles
VaR_level <- alpha #manually set the quantile level as alpha
DQtest(out_of_sample_returns, outsamcav,VaR_level)
DQtest(out_of_sample_returns, outsamcavRQ_05,VaR_level)
DQtest(out_of_sample_returns, outsamcavRQ_10,VaR_level)
DQtest(out_of_sample_returns, outsamcavRV,VaR_level)
DQtest(out_of_sample_returns, outsamcavvix,VaR_level)

##The Tick loss function for the out-of-sample returns
tick_loss <- function(y, var, alpha) {
  loss <- ifelse(y >= var, 
                 (1 - alpha) * (y - var), 
                 alpha * (var - y)
  )
  return(loss)
}

loss_values <- tick_loss(out_of_sample_returns, outsamcav,alpha)
mean_tick_loss <- mean(loss_values)
mean_tick_loss
loss_valuesvix <- tick_loss(out_of_sample_returns, outsamcavvix, alpha)
mean_tick_lossvix <- mean(loss_valuesvix)

loss_valuesRV <- tick_loss(out_of_sample_returns, outsamcavRV, alpha)
mean_tick_lossRV <- mean(loss_valuesRV)

loss_valuesRQ05 <- tick_loss(out_of_sample_returns, outsamcavRQ_05, alpha)
mean_tick_lossRQ05 <- mean(loss_valuesRQ05)


loss_valuesRQ10 <- tick_loss(out_of_sample_returns, outsamcavRQ_10, alpha)
mean_tick_lossRQ10 <- mean(loss_valuesRQ10)

#table for tick loss
tick_loss_list <- c(mean_tick_loss, mean_tick_lossvix, mean_tick_lossRV, mean_tick_lossRQ10, mean_tick_lossRQ05)
names(tick_loss_list) <- c("Base-Caviar", "Caviarx (VIX)", "Caviarx (RV)", "Caviarx (RQ10)", "Caviarx (RQ05)")
kable(tick_loss_list, caption = "Mean out-of-sample tick loss")

tickloss <- data.frame(Time = out_of_sample_returns, TickLoss = loss_values, TickLossVIX = loss_valuesvix, TickLossRV = loss_valuesRV,  TickLossRQ05 = loss_valuesRQ05)

# Define the Diebold-Mariano test function (taken from the lecture notes from Prof. Dimitriades, "Predictive Modelling")
DM_test <- function(q1, q2, y, tau) {
  
  # Define the asymmetric piecewise linear (APL) loss function
  apl_score <- function(q, y, tau) {
    (1 - tau) * (q - y) * (q > y) + tau * (y - q) * (y >= q)
  }
  
  # Calculate the score difference between two sets of quantile forecasts
  score_diff <- apl_score(q1, y, tau) - apl_score(q2, y, tau)
  
  # Fit a regression model of the score difference on the intercept (constant)
  dm_model <- lm(score_diff ~ 1)
  
  # Apply Newey-West correction for heteroskedasticity and autocorrelation
  DM_quant <- coeftest(dm_model, vcov = NeweyWest(dm_model))
  
  return(DM_quant)  # Return the corrected test result
}

DM_test(outsamcav, outsamcavvix, out_of_sample_returns, alpha)
DM_test(outsamcavvix, outsamcavRV, out_of_sample_returns, alpha)
DM_test(outsamcavRQ_10, outsamcavRQ_05, out_of_sample_returns, alpha)
DM_test(outsamcav, outsamcavRQ_10, out_of_sample_returns, alpha)
DM_test(outsamcavRQ_05, outsamcavRV, out_of_sample_returns, alpha)
DM_test(outsamcavRQ_05, outsamcavvix, out_of_sample_returns, alpha)
DM_test(outsamcavRQ_10, outsamcavvix, out_of_sample_returns, alpha)
DM_test(outsamcavRQ_10, outsamcavRV, out_of_sample_returns, alpha)


#Mincer-Zarnowitz test (taken from the lecture notes from Prof. Dimitriades, "Predictive Modelling")
MZ.quant.test <- function(haty, y, tau=tau){
  MZ.QR.fit <- quantreg::rq(y~haty, tau=tau) 
  # se="boot" uses a bootstrap for the covariance (for simplicity)
  set.seed(1)
  Summary.MZ.QR.fit <- summary(MZ.QR.fit, covariance=TRUE, se="boot")
  
  W <- t(coef(MZ.QR.fit) - c(0,1)) %*% solve(Summary.MZ.QR.fit$cov) %*% (coef(MZ.QR.fit) - c(0,1))
  pval <- 1- pchisq(W,2) %>% round(3)
  
  # Add a respective plot
  p <- ggplot(data=data.frame(y=y, haty=haty), aes(x=haty, y=y)) +
    geom_abline(slope=1, intercept=0, col="black") +
    geom_point(color="grey50") +
    geom_quantile(quantiles=tau, color="red")
  
  return(list(estimate=Summary.MZ.QR.fit, WaldStat=as.numeric(W), Wald.pval=as.numeric(pval), plot=p))
}
tau=alpha
MZ.quant.test(outsamcav, out_of_sample_returns, tau)
MZ.quant.test(outsamcavvix, out_of_sample_returns, tau)
MZ.quant.test(outsamcavRV, out_of_sample_returns, tau)
MZ.quant.test(outsamcavRQ_10, out_of_sample_returns, tau)
MZ.quant.test(outsamcavRQ_05, out_of_sample_returns, tau)



#Forecast Evaluation for GARCH models

#Liklihood test for garch parameters
lr_test <- function(loglik_null, loglik_alt, df) { #df is the degrees of freedom
  LR_stat <- -2 * (sum(loglik_null) - sum(loglik_alt))
  p_value <- pchisq(LR_stat, df = df, lower.tail = FALSE)
  return(list(LR_stat = LR_stat, p_value = p_value))
}
#We take null hypothesis as the GARCHx-RQ model as the liklihood value is the highest for this model specification.
null<-(-garxrq05$value)

lr_test_vix<-lr_test(null,-garxvix$value,3)
lr_test_rq10<-lr_test(null,-garxrq10$value,3)
lr_test_rv<-lr_test(null,-garxRV$value,3)
lr_test_base<-lr_test(null,-garch$value,2)
#To make table for liklihood ratio statistic
lr_table<-data.frame("Base"=lr_test_base$LR_stat, "VIX"=lr_test_vix$LR_stat, "RV"=lr_test_rv$LR_stat, "RQ_10"=lr_test_rq10$LR_stat, "RQ05"=0)
lr_table

#Unconditional Coverage test for the GARCH models for alpha 5%
uc_test<- function(returns, VaR, alpha) {
  n <- length(returns)
  violations <- returns < VaR
  x <- sum(violations)
  obs_violation_rate <- x / n
  LR_uc <- -2 * (x * log(alpha) + (n - x) * log(1 - alpha) - 
                   (x * log(obs_violation_rate) + (n - x) * log(1 - obs_violation_rate)))
  p_value <- 1 - pchisq(LR_uc, df = 1)
  
  chi_squared_critical <- qchisq(alpha, df = 1) #against the confidence level
  test_result <- ifelse(LR_uc < alpha, "Reject Null", "Fail to Reject Null")
  
  return(list(
    LR_uc = LR_uc,
    p_value = p_value
  ))
}
uc_test_garch<-uc_test(out_of_sample_returns, var_garch, alpha)
uc_test_garx<-uc_test(out_of_sample_returns, var_garx, alpha)
uc_test_garxrq10<-uc_test(out_of_sample_returns, var_garxrq10, alpha)
uc_test_garxrq05<-uc_test(out_of_sample_returns, var_garxrq05, alpha)
uc_test_garxRV<-uc_test(out_of_sample_returns, var_garxRV, alpha)

uc_table<-data.frame("Base"=uc_test_garch$LR_uc, "VIX"=uc_test_garx$LR_uc, "RV"=uc_test_garxRV$LR_uc, "RQ_10"=uc_test_garxrq10$LR_uc, "RQ_05"=uc_test_garxrq05$LR_uc)
uc_table
p_table<-data.frame("Base"=uc_test_garch$p_value, "VIX"=uc_test_garx$p_value, "RV"=uc_test_garxRV$p_value, "RQ_10"=uc_test_garxrq10$p_value, "RQ_05"=uc_test_garxrq05$p_value)
p_table


#Tick loss function for GARCHx models
tick_loss <- function(y, var, alpha) {
  loss <- ifelse(y >= var, 
                 (1 - alpha) * (y - var), 
                 alpha * (var - y)
  )
  return(loss)
}
lossgar<-tick_loss(out_of_sample_returns, var_garch, alpha)
lossgarx<-tick_loss(out_of_sample_returns, var_garx, alpha)
lossgarxrq10<-tick_loss(out_of_sample_returns, var_garxrq10, alpha)
lossgarxrq05<-tick_loss(out_of_sample_returns, var_garxrq05, alpha)
lossgarxRV<-tick_loss(out_of_sample_returns, var_garxRV, alpha)


meanlossgar<-mean(lossgar)
meanlossgarx<-mean(lossgarx)
meanlossgarxrq10<-mean(lossgarxrq10)
meanlossgarxrq05<-mean(lossgarxrq05)
meanlossgarxRV<-mean(lossgarxRV)
meanloss_table<-data.frame("Base"=meanlossgar, "VIX"=meanlossgarx, "RV"=meanlossgarxRV, "RQ_10"=meanlossgarxrq10, "RQ_05"=meanlossgarxrq05)
meanloss_table

#Diebold-Mariano test for GARCH models
DM_test(var_garx, var_garch, out_of_sample_returns, alpha) #can be used for all the models
DM_test(var_garx, var_garxRV, out_of_sample_returns, alpha)
DM_test(var_garx, outsamcavvix, out_of_sample_returns, alpha)
DM_test(outsamcavRQ_10, var_garxrq10, out_of_sample_returns, alpha)
#Mincer-Zarnowitz test
MZ.quant.test(var_garch, out_of_sample_returns, tau)
MZ.quant.test(var_garx, out_of_sample_returns, tau)
MZ.quant.test(var_garxRV, out_of_sample_returns, tau)
MZ.quant.test(var_garxrq10, out_of_sample_returns, tau)
MZ.quant.test(var_garxrq05, out_of_sample_returns, tau)


