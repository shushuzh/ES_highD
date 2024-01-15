library(Matrix)
library(conquer)
library(glmnet)
library(MultiRNG)
library(dplyr)
library(janitor)

# main function: compute ES estimation
highdim_2step = function(x,y,alpha){
  # Function to calculate two-step high-dim ES Estimation
  #
  # Parameters:
  #   x: covariates (n*p matrix without intercept);
  #   y: observations;
  #   alpha: quantile level for ES;
  #  
  # Returns:
  #   theta0_hat: ES estimation (intercept)
  #   theta0_hat: ES estimation (slopes)
  #   resid: residuals from the second step (Y_i - X_i^T \hat \beta) 1(Y_i \le X_i^T \hat \beta)/alpha + X_i^T \hat \beta - X_i \hat \theta
  
  # step 1: fit quantile regression for estimating beta
  quan_model = conquer.cv.reg(x,y,
                              #lambda=0.03,
                              tau=alpha,
                              h=max(0.05,sqrt(alpha*(1-alpha))*(log(p)/n)^(1/4)))
  lambda_beta = quan_model$lambda
  beta_hat = quan_model$coeff.min[-1]
  beta0 = quan_model$coeff.min[1]
  # step 2: fit ES regression for estimating theta
  Z_2step = apply(cbind(y,x),MARGIN = 1,FUN = Z_beta,beta=beta_hat,alpha=alpha,beta_0=beta0)
  ## CV: perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(x, Z_2step,standardize=FALSE)
  ## find optimal lambda value that minimizes test MSE
  lambda_theta <- cv_model$lambda.min
  theta_model = glmnet(x, Z_2step,lambda=lambda_theta,standardize=FALSE)
  theta_hat = theta_model$beta
  theta0_hat = theta_hat[1]
  resid = Z_2step - predict(theta_model,newx=x)
  return(list(theta0_hat = theta0_hat, theta_hat = theta_hat,resid = resid, lambda_beta = lambda_beta, lambda_theta = lambda_theta))
}

# main function: compute ES inference
highdim_inf = function(x,y,alpha,col = 1,res_est=NULL){
  # Inference of high-dim ES on one parameter of interest using debiased idea
  #
  # Parameters:
  #   x: covariates (n*p matrix without intercept);
  #   y: observations;
  #   alpha: quantile level for ES;
  #   col: the index of the variable that we are interested in inferencing (scaler);
  #   res_est: results from two-step estimation, i.e., return from function "highdim_2step", 
  #            NULL if two-step estimation was not called;
  #  
  # Returns:
  #   theta_debias: debiased estimator for the ES estimator at the parameter of interest (i.e., col)
  #   Conf.int: confidence interval for the debiased ES estimator
  
  if (is.null(res_est)){
    res_est = highdim_2step(x,y,alpha)
  }
  theta0_hat = res_est$theta0_hat
  theta_hat = res_est$theta_hat
  resid1 = res_est$resid
  
  # standardize x
  x_stand = sweep(x, 2, colMeans(x))
  d = x_stand[,col]
  x_tilde = x_stand[,-col]
  ### CV ###
  cv_model <- cv.glmnet(x_tilde, d,intercept=FALSE,standardize=FALSE)
  lambda_cv <- cv_model$lambda.1se
  gamma_model = glmnet(x_tilde, d, lambda=lambda_cv,intercept=FALSE,standardize=FALSE)
  gamma_hat = gamma_model$beta
  nonzero_index = which(gamma_hat!=0)
  resid2 = d - predict(gamma_model,newx=x_tilde)
  w_hat = resid2
  theta_debias = theta_hat[col] + t(resid1)%*%resid2/x[,col]%*%resid2
  #########################Variance Estimation using RCV################################
  ### split sample into two
  sample_index = sample(seq_len(nrow(x)),size = ceiling(n/2))
  d_x_1 =x[sample_index,]
  d_x_2 =x[-sample_index,]
  y_1 = y[sample_index]
  y_2 = y[-sample_index]
  wid = variance_estimation(d_x_1,d_x_2,y_1,y_2)
  return(list(theta_debias = theta_debias, Conf.int = c(theta_debias-wid,theta_debias+wid)))
}

# supp function: compute RCV variance estimator
variance_estimation=function(d_x_1,d_x_2,y_1,y_2,conf.level = 0.95){
  # Remove the columns that are constant
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
  beta_hat = quan_model$coeff.min[-1]
  beta0 = quan_model$coeff.min[1]
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
  dx1_model = glmnet(d_x_1_stand[,-1], d_x_1_stand[,1], lambda=lambda_cv,intercept=FALSE,standardize=FALSE)
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
  beta_hat = quan_model$coeff.min[-1]
  beta0 = quan_model$coeff.min[1]
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
  CI_level = 1-(1-conf.level)/2
  wid = sqrt(sigma_s_hat)*qnorm(CI_level)/sqrt(n)/alpha/sigma_w_hat_square
  return(wid)
}

# supp function: compute adjusted response for the second step
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

