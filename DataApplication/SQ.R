setwd("/Users/shushuz/Dropbox (University of Michigan)/ES_HighD/Empirical_studies/DataApplication")
library(haven)
library(quantreg)
library(xtable)
# For data wrangling
#dplyr::select, mutate, select, recode
library(dplyr)
library(survey)
library(conquer)
library(glmnet)
library(janitor)

variance_estimation=function(d_x_1,d_x_2,y_1,y_2){
  # To find the columns that are constant
  constant_columns_1 = data.frame(d_x_1) %>%
    dplyr::select_at(setdiff(names(.), names(remove_constant(.)))) %>%
    unique()
  constant_columns_2 = data.frame(d_x_2) %>%
    dplyr::select_at(setdiff(names(.), names(remove_constant(.)))) %>%
    unique()
  constant_col = unique(c(colnames(constant_columns_1),colnames(constant_columns_2)))
  d_x_1 = as.matrix(dplyr::select(data.frame(d_x_1), -constant_col))
  d_x_2 = as.matrix(dplyr::select(data.frame(d_x_2), -constant_col))
  
  #######find variables on first half######## 
  #######three step procedure#########
  # step 1
  quan_model = conquer.cv.reg(d_x_1,y_1,
                              #lambda=0.03,
                              tau=alpha,
                              h=max(0.05,sqrt(alpha*(1-alpha))*(log(p)/n)^(1/4)))
  beta_hat = quan_model$coeff[-1]
  beta0 = quan_model$coeff[1]
  sq = sum(beta_hat!=0)
  # step 2
  Z_2step = apply(cbind(y_1,d_x_1),MARGIN = 1,FUN = Z_beta,beta=beta_hat,alpha=alpha,beta_0=beta0)
  ######### CV #########
  ##perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(d_x_1, Z_2step,standardize=FALSE)
  ##find optimal lambda value that minimizes test MSE
  lambda_theta <- cv_model$lambda.min
  theta_model = glmnet(d_x_1, Z_2step,lambda=lambda_theta,standardize=FALSE)
  theta_hat = theta_model$beta
  se = length(theta_hat@i)
  # step 3
  d_x_1_stand = sweep(d_x_1, 2, colMeans(d_x_1))
  cv_model <- cv.glmnet(d_x_1_stand[,-1], d_x_1_stand[,1],intercept=FALSE,standardize=FALSE)
  lambda_cv <- cv_model$lambda.1se
  dx1_model = glmnet(d_x_1[,-1], d_x_1[,1], lambda=lambda_cv,intercept=FALSE,standardize=FALSE)
  gamma_hat = dx1_model$beta
  nonzero_index = which(gamma_hat!=0)
  #print(nonzero_index)
  #nonzero_index = gamma_hat@i+2
  sm = length(nonzero_index)
  ##########fit on second half###################
  # step 1
  dx2_short = d_x_2[,beta_hat!=0]
  quan_model <- conquer(dx2_short,y_2,tau=alpha)
  beta_hat = quan_model$coeff[-1]
  beta0 = quan_model$coeff[1]
  epsilon_hat = y_2 - dx2_short%*%beta_hat - beta0
  epsilon_minus = ifelse(epsilon_hat<0,epsilon_hat,0)
  # step 2
  Z_2step = apply(cbind(y_2,dx2_short),MARGIN = 1,FUN = Z_beta,beta=beta_hat,alpha=alpha,beta_0=beta0)
  dx2_short_theta = d_x_2[,theta_hat@i+1]
  theta_model = lm(Z_2step~dx2_short_theta)
  theta_hat = theta_model$coefficients[-1]
  res3 = theta_model$residuals
  # step 3
  d_x_2_stand = sweep(d_x_2, 2, colMeans(d_x_2))
  if(sm==0){
    w_hat = d_x_2_stand[,1]
  } else {
    lm_model <- lm(d_x_2_stand[,1]~d_x_2_stand[,nonzero_index+1]+0)
    w_hat = lm_model$residuals}
  
  #results
  s2 = nrow(d_x_2)
  sigma_w_hat_square_1 = sum(w_hat^2)/(s2-sm-1)
  sigma_s_hat_1 = alpha^2*t(w_hat^2)%*%((res3)^2)/(s2-sm-sq-se-3)
  #######################The other way around#####################
  #######find variables on first half######## 
  #######three step procedure#########
  # step 1
  quan_model = conquer.cv.reg(d_x_2,y_2,
                              #lambda=0.03,
                              tau=alpha,
                              h=max(0.05,sqrt(alpha*(1-alpha))*(log(p)/n)^(1/4)))
  beta_hat = quan_model$coeff[-1]
  beta0 = quan_model$coeff[1]
  sq = sum(beta_hat!=0)
  # step 2
  Z_2step = apply(cbind(y_2,d_x_2),MARGIN = 1,FUN = Z_beta,beta=beta_hat,alpha=alpha,beta_0=beta0)
  ######### CV #########
  ##perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(d_x_2, Z_2step,standardize=FALSE)
  ##find optimal lambda value that minimizes test MSE
  lambda_theta <- cv_model$lambda.min
  theta_model = glmnet(d_x_2, Z_2step,lambda=lambda_theta,standardize=FALSE)
  theta_hat = theta_model$beta
  se = length(theta_hat@i)
  # step 3
  cv_model <- cv.glmnet(d_x_2_stand[,-1], d_x_2_stand[,1],intercept=FALSE,standardize=FALSE)
  lambda_cv <- cv_model$lambda.1se
  dx2_model = glmnet(d_x_2_stand[,-1], d_x_2_stand[,1], lambda=lambda_cv,intercept=FALSE,standardize=FALSE)
  gamma_hat = dx2_model$beta
  nonzero_index = which(gamma_hat!=0)
  #print(nonzero_index)
  #nonzero_index = gamma_hat@i+2
  sm = length(nonzero_index)
  ##########fit on second half###################
  # step 1
  dx1_short = d_x_1[,beta_hat!=0]
  quan_model <- conquer(dx1_short,y_1,tau=alpha)
  beta_hat = quan_model$coeff[-1]
  beta0 = quan_model$coeff[1]
  epsilon_hat = y_1 - dx1_short%*%beta_hat - beta0
  epsilon_minus = ifelse(epsilon_hat<0,epsilon_hat,0)
  # step 2
  Z_2step = apply(cbind(y_1,dx1_short),MARGIN = 1,FUN = Z_beta,beta=beta_hat,alpha=alpha,beta_0=beta0)
  dx1_short_theta = d_x_1[,theta_hat@i+1]
  theta_model = lm(Z_2step~dx1_short_theta)
  theta_hat = theta_model$coefficients[-1]
  res3 = theta_model$residuals
  # step 3
  if(sm==0){
    w_hat = d_x_1_stand[,1]
  } else {
    lm_model <- lm(d_x_1_stand[,1]~d_x_1_stand[,nonzero_index+1]+0)
    w_hat = lm_model$residuals}
  
  #results
  s1 = nrow(d_x_1)
  sigma_w_hat_square_2 = sum(w_hat^2)/(s1-sm-1)
  sigma_s_hat_2 = alpha^2*t(w_hat^2)%*%((res3)^2)/(s1-sm-sq-se-3)
  sigma_w_hat_square = (sigma_w_hat_square_1+sigma_w_hat_square_2)/2
  sigma_s_hat = (sigma_s_hat_1+sigma_s_hat_2)/2
  wid = sqrt(sigma_s_hat)*qnorm(0.975)/sqrt(n)/alpha/sigma_w_hat_square
  return(wid)
}

Z_beta = function(data,beta,alpha,beta_0){#data:each row of cbind(y,x)
  yi = data[1]
  xi = data[-1]
  res = yi-xi%*%beta-beta_0
  if (res<=0){
    return (res/alpha+xi%*%beta+beta_0)
  }
  else{
    return (xi%*%beta+beta_0)
  }
}

######### Two Step Regression ######
alpha = 0.20
design_matrix = read.csv("./design_matrix_new.csv")
y = -design_matrix[,"cotinine"]
x = subset(design_matrix,select = -c(cotinine))#,,overnight2,mental2
x = x[,-grep("HSQ*",colnames(x))] 
x = as.matrix(x)
p = ncol(x)
n = nrow(x)

cbind(rownames(theta_hat)[which(theta_hat!=0)],theta_hat[which(theta_hat!=0)])

############ES Estimation#######
# step 1: fit quantile regression for estimating beta
quan_model = conquer.cv.reg(x,y,
                            #lambda=0.03,
                            tau=alpha,
                            h=max(0.05,sqrt(alpha*(1-alpha))*(log(p)/n)^(1/4)))
beta_hat = quan_model$coeff[-1]
beta_0 = quan_model$coeff[1]
print(beta_hat[1:3])

####step 2###
Z_2step = apply(cbind(y,x),MARGIN = 1,FUN = Z_beta,beta=beta_hat,alpha=alpha,beta_0=beta_0)
lambda_max <- max(abs(colSums(x*Z_2step)))/n/1000
epsilon <- .0001
K <- 100
lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon),
                            length.out = K)), digits = 10)
cv_model <- cv.glmnet(x, Z_2step,lambda=lambdapath,standardize=FALSE)#,penalty.factor=c(rep(0,3),rep(1,ncol(x)-3)))
plot(cv_model)
lambda_theta <- cv_model$lambda.min
theta_model = glmnet(x, Z_2step,lambda=lambda_theta,standardize=FALSE)#,penalty.factor=c(rep(0,3),rep(1,ncol(x)-3)))
theta_hat = theta_model$beta
resid1 = Z_2step - predict(theta_model,newx=x)
print(theta_hat[1:3])
theta_refit_model = lm(Z_2step~x[,which(theta_hat!=0)])
theta_refit = theta_refit_model$coefficients[2:4]
print(theta_refit)

######## (1) Asian v.s. White ######
#x_stand = x - colMeans(x)
x_stand = sweep(x, 2, colMeans(x))
d = x_stand[,1]
x_tilde = x_stand[,-c(1,2,3)]
#### supply larger lambda sequence ###
lambda_max <- max(abs(colSums(x_tilde*d)))/n/100
epsilon <- .0001
K <- 100
lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon),
                            length.out = K)), digits = 10)
### CV ###
cv_model <- cv.glmnet(x_tilde, d,intercept=FALSE,standardize=FALSE)
plot(cv_model)
lambda_cv <- cv_model$lambda.1se
gamma_model = glmnet(x_tilde, d, lambda=lambda_cv,intercept=FALSE,standardize=FALSE)
gamma_hat = gamma_model$beta
nonzero_index = which(gamma_hat!=0)
print(nonzero_index)
resid2 = d - predict(gamma_model,newx=x_tilde)
w_hat = resid2
if(length(nonzero_index)==0){resid2_refit = resid2
} else {
  gamma_model_refit = lm(d~x_tilde[,nonzero_index]+0)
  resid2_refit = d - predict(gamma_model_refit,newx=x_tilde)}

theta_0_debias = theta_hat[1] + t(resid1)%*%resid2/x[,1]%*%resid2
theta_0_debias_refit = theta_hat[1] + t(resid1)%*%resid2_refit/x[,1]%*%resid2_refit

d = x[,1]
x_tilde = x[,-c(1,2,3)]
cv_model <- cv.glmnet(x_tilde, d)
plot(cv_model)
lambda_cv <- cv_model$lambda.1se
gamma_model = glmnet(x_tilde, d,lambda=lambda_cv)#,standardize=TRUE
gamma_hat = gamma_model$beta
resid2 = d - predict(gamma_model,newx=x_tilde)
theta_1_debias = theta_hat[1] + t(resid1)%*%resid2/d%*%resid2
print(theta_1_debias)

#Variance Estimation#
sample_index = sample(seq_len(nrow(x)),size = ceiling(n/2))
d_x_1 =x[sample_index,-c(2,3)]
d_x_2 =x[-sample_index,-c(2,3)]
y_1 = y[sample_index]
y_2 = y[-sample_index]
wid1 = variance_estimation(d_x_1,d_x_2,y_1,y_2)
LB1 = theta_1_debias-wid1
UB1 = theta_1_debias+wid1
print(LB1);print(UB1)

######## (2) Black v.s. White ########
x_stand = sweep(x, 2, colMeans(x))
d = x_stand[,2]
x_tilde = x_stand[,-c(1,2,3)]
### CV ###
cv_model <- cv.glmnet(x_tilde, d,standardize=FALSE,intercept=FALSE)
plot(cv_model)
lambda_cv <- cv_model$lambda.1se
gamma_model = glmnet(x_tilde, d, lambda=lambda_cv,intercept=FALSE,standardize=FALSE)
gamma_hat = gamma_model$beta
nonzero_index = which(gamma_hat!=0)
resid2 = d - predict(gamma_model,newx=x_tilde)
w_hat = resid2
if(length(nonzero_index)==0){resid2_refit = resid2
} else {
  gamma_model_refit = lm(d~x_tilde[,nonzero_index]+0)
  resid2_refit = d - predict(gamma_model_refit,newx=x_tilde)}

theta_2_debias = theta_hat[2] + t(resid1)%*%resid2/x[,2]%*%resid2
print(theta_2_debias)

#Variance Estimation#
sample_index = sample(seq_len(nrow(x)),size = ceiling(n/2))
d_x_1 =x[sample_index,-c(1,3)]
d_x_2 =x[-sample_index,-c(1,3)]
y_1 = y[sample_index]
y_2 = y[-sample_index]

wid2 = variance_estimation(d_x_1,d_x_2,y_1,y_2)
LB2 = theta_2_debias-wid2
UB2 = theta_2_debias+wid2
print(LB2);print(UB2)

######## (3) Latino v.s. White ########
x_stand = sweep(x, 2, colMeans(x))
d = x_stand[,3]
x_tilde = x_stand[,-c(1,2,3)]
### CV ###
cv_model <- cv.glmnet(x_tilde, d,standardize=FALSE,intercept=FALSE)
plot(cv_model)
lambda_cv <- cv_model$lambda.1se
gamma_model = glmnet(x_tilde, d, lambda=lambda_cv,intercept=FALSE,standardize=FALSE)
gamma_hat = gamma_model$beta
nonzero_index = which(gamma_hat!=0)
resid2 = d - predict(gamma_model,newx=x_tilde)
w_hat = resid2
if(length(nonzero_index)==0){resid2_refit = resid2
} else {
  gamma_model_refit = lm(d~x_tilde[,nonzero_index]+0)
  resid2_refit = d - predict(gamma_model_refit,newx=x_tilde)}

theta_3_debias = theta_hat[3] + t(resid1)%*%resid2/x[,3]%*%resid2
print(theta_3_debias)

#Variance Estimation#
sample_index = sample(seq_len(nrow(x)),size = ceiling(n/2))
d_x_1 =x[sample_index,-c(1,2)]
d_x_2 =x[-sample_index,-c(1,2)]
y_1 = y[sample_index]
y_2 = y[-sample_index]

wid3 = variance_estimation(d_x_1,d_x_2,y_1,y_2)
LB3 = theta_3_debias-wid3
UB3 = theta_3_debias+wid3
print(LB3);print(UB3)

#### compare with low-dim estimation ####
quan_model <- conquer(x,y,tau=alpha)
beta_hat = quan_model$coeff[-1]
beta0 = quan_model$coeff[1]
epsilon_hat = y - x%*%beta_hat - beta0
epsilon_minus = ifelse(epsilon_hat<0,epsilon_hat,0)
# step 2
Z_2step = apply(cbind(y,x),MARGIN = 1,FUN = Z_beta,beta=beta_hat,alpha=alpha,beta_0=beta0)
theta_model = lm(Z_2step~x)
theta_hat = theta_model$coefficients[-1]
theta0 = theta_model$coefficients[1][[1]]
res3 = theta_model$residuals
x_new = x[,!is.na(theta_hat)]
omega = epsilon_minus + alpha*(x_new%*%(beta_hat[!is.na(theta_hat)] - as.numeric(theta_hat[!is.na(theta_hat)]))+beta0 - theta0)
sigma_hat = t(x_new)%*%x_new/n
omega_hat = t(x_new)%*%diag(as.vector(omega^2))%*%x_new/n
sandwich = solve(sigma_hat)%*%omega_hat%*%solve(sigma_hat)
c(theta_hat[1] - 1.96/alpha/sqrt(n)*sqrt(sandwich[1,1]),
theta_hat[1] + 1.96/alpha/sqrt(n)*sqrt(sandwich[1,1])) #41.81836 -115.8880  199.5247 
c(theta_hat[2] - 1.96/alpha/sqrt(n)*sqrt(sandwich[2,2]),
  theta_hat[2] + 1.96/alpha/sqrt(n)*sqrt(sandwich[2,2])) #-65.02886 -195.07549 65.01777 
c(theta_hat[3] - 1.96/alpha/sqrt(n)*sqrt(sandwich[3,3]),
  theta_hat[3] + 1.96/alpha/sqrt(n)*sqrt(sandwich[3,3])) #24.85934 -109.6969  159.4156

####Summary of Results
print(theta_1_debias);print(LB1);print(UB1)
print(theta_2_debias);print(LB2);print(UB2)
print(theta_3_debias);print(LB3);print(UB3)
###### compare with LS lasso#####
### desparsified lasso
library(desla)
H = c(1,2,3)
despar = desla(x,y,H)#,alphas=0.01)
despar$intervals
