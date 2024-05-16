library(Matrix)
library(conquer)
library(glmnet)
library(MultiRNG)
library(dplyr)
library(janitor)

# main function: compute ES estimation
highdim_2step = function(x,y,alpha,standardize=T,tuning.method = "CV"){
  # Function to calculate two-step high-dim ES Estimation
  #
  # Parameters:
  #   x: covariates (n*p matrix without intercept);
  #   y: observations;
  #   alpha: quantile level for ES;
  #   standardize: (optional) Logical flag for x variable standardization, prior to fitting the model sequence. The coefficients are always returned on the original scale. Default is standardize=TRUE;
  #   tuning.method: "CV" (for cross-validation) or "BIC";
  
  # Returns:
  #   theta0_hat: ES estimation (intercept)
  #   theta0_hat: ES estimation (slopes)
  #   resid: residuals from the second step (Y_i - X_i^T \hat \beta) 1(Y_i \le X_i^T \hat \beta)/alpha + X_i^T \hat \beta - X_i \hat \theta
  
  if (tuning.method == "BIC"){
    quant_model_bic = quantile_bic(x,y,alpha,standardize=standardize)
    lambda_beta = quant_model_bic$lambda_select
    model_select = quant_model_bic$model_select
    beta_hat = model_select$coeff[-1]
    beta0 = model_select$coeff[1]
    # step 2: fit ES regression for estimating theta
    #Z_2step = beta0+x%*%beta_hat + (y-beta0-x%*%beta_hat)*ifelse(y<beta0+x%*%beta_hat,1,0)/alpha
    Z_2step = apply(cbind(y,x),MARGIN = 1,FUN = Z_beta,beta=beta_hat,alpha=alpha,beta_0=beta0)
    ES_model = LS_bic(x,Z_2step,standardize=standardize)
    lambda_theta = ES_model$lambda_select
    theta_model = ES_model$model_select
  } else {
    # step 1: fit quantile regression for estimating beta
    quan_model = conquer.cv.reg(x,y,
                                #lambda=0.03,
                                tau=alpha,
                                h=max(0.05,sqrt(alpha*(1-alpha))*(log(p)/n)^(1/4)))
    lambda_beta = quan_model$lambda
    beta_hat = quan_model$coeff.min[-1]
    beta0 = quan_model$coeff.min[1]
    # step 2: fit ES regression for estimating theta
    Z_2step = beta0+x%*%beta_hat + (y-beta0-x%*%beta_hat)*ifelse(y<beta0+x%*%beta_hat,1,0)/alpha
    #Z_2step = apply(cbind(y,x),MARGIN = 1,FUN = Z_beta,beta=beta_hat,alpha=alpha,beta_0=beta0)
    ## CV: perform k-fold cross-validation to find optimal lambda value
    cv_model <- cv.glmnet(x, Z_2step,standardize=standardize)
    ## find optimal lambda value that minimizes test MSE
    lambda_theta <- cv_model$lambda.min
    theta_model = glmnet(x, Z_2step,lambda=lambda_theta,standardize=standardize)
}
  theta_hat = theta_model$beta
  theta0_hat = theta_model$a0
  resid = Z_2step - predict(theta_model,newx=x)
  return(list(theta0_hat = theta0_hat, theta_hat = theta_hat,beta0_hat = beta0, beta_hat = beta_hat,resid = resid, 
              lambda_beta = lambda_beta, lambda_theta = lambda_theta))
}

# main function: compute ES inference
highdim_inf = function(x,y,alpha,
                       col = 1,res_est=NULL,
                       conf.level = 0.95, standardize=T,
                       variance.method = "RCV"){
  # Inference of high-dim ES on one parameter of interest using debiased idea
  #
  # Parameters:
  #   x: covariates (n*p matrix without intercept);
  #   y: observations;
  #   alpha: quantile level for ES;
  #   col: the index of the variable that we are interested in inferencing (scaler);
  #   res_est: results from two-step estimation, i.e., return from function "highdim_2step", 
  #            NULL if two-step estimation was not called;
  #   conf.level: confidence level;
  #   variance.method: Methods for estimating the asymptotic variance. Choose from "RCV" (refitted cross-validation), "refit" (valinilla refit), and "plug-in" 
  
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
  beta0 = res_est$beta0_hat
  beta_hat = res_est$beta_hat
  
  d = x[,col]
  x_tilde = x[,-col]
  ### CV ###
  cv_model <- cv.glmnet(x_tilde, d,standardize=standardize)
  lambda_cv <- cv_model$lambda.1se
  gamma_model = glmnet(x_tilde, d, lambda=lambda_cv,standardize=standardize)
  gamma_hat = gamma_model$beta
  nonzero_index = which(gamma_hat!=0)
  resid2 = d - predict(gamma_model,newx=x_tilde)
  w_hat = resid2
  theta_debias = theta_hat[col] + t(resid1)%*%resid2/x[,col]%*%resid2
  sq = sum(beta_hat!=0)
  se = length(theta_hat@i)
  sm = length(nonzero_index)
  CI_level = 1-(1-conf.level)/2
  ############ Variance estimation using simple plug-in ################
  if (variance.method == "plug-in"){
    sigma_w_hat_square = sum(w_hat^2)/(n-sm-1)
    sigma_s_hat = alpha^2*t(w_hat^2)%*%((resid1)^2)/(n-sm-sq-se-3)
    var_plugin = sigma_s_hat/sigma_w_hat_square^2
    wid = sqrt(var_plugin)*qnorm(CI_level)/sqrt(n)/alpha
  } 
  else if (variance.method == "refit"){
    #########################Variance Estimation using vanilla refit################################
    Z_2step = apply(cbind(y,x),MARGIN = 1,FUN = Z_beta,beta=beta_hat,alpha=alpha,beta_0=beta0)
    if(se==0){
      res3 = Z_2step - mean(Z_2step)
    } else {
      x_refit = x[,theta_hat@i+1]
      theta_model_refit = lm(Z_2step~x_refit)
      theta_hat_refit = theta_model_refit$coefficients[-1]
      res3 = theta_model_refit$residuals
    }
    if(sm==0){
      w_hat = d - mean(d)
    } else {
      lm_model <- lm(d~x_tilde[,nonzero_index])
      w_hat = lm_model$residuals}
    sigma_w_hat_square = sum(w_hat^2)/(n-sm-1)
    sigma_s_hat = alpha^2*t(w_hat^2)%*%((res3)^2)/(n-sm-sq-se-3)
    var_vanilla = sigma_s_hat/sigma_w_hat_square^2
    wid = sqrt(var_vanilla)*qnorm(CI_level)/sqrt(n)/alpha
  } else {
  #########################Variance Estimation using RCV################################
  ### split sample into two
  sample_index = sample(seq_len(nrow(x)),size = ceiling(n/2))
  d_x_1 =x[sample_index,]
  d_x_2 =x[-sample_index,]
  y_1 = y[sample_index]
  y_2 = y[-sample_index]
  var = variance_estimation(d_x_1,d_x_2,y_1,y_2,alpha,standardize=standardize)
  wid = sqrt(var)*qnorm(CI_level)/sqrt(n)/alpha
  }
  return(list(theta_debias = theta_debias, Conf.int = c(theta_debias-wid,theta_debias+wid)))
}

# supp function: compute RCV variance estimator
variance_estimation=function(d_x_1,d_x_2,y_1,y_2,alpha,standardize=T){
  # asymptotic variance estimation of the high-dim ES on one parameter of interest using refitted cross-validation
  #
  # Parameters:
  #   d_x_1: first part of the covariates;
  #   d_x_2: second part of the covariates;
  #   y_1: first part of the observations (corresponding with d_x_1);
  #   y_2: second part of the observations (corresponding with d_x_2);
  #   alpha: quantile level for ES;
  #  
  # Returns:
  #   wid: width of the confidence interval
  
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
  cv_model <- cv.glmnet(d_x_1, Z_2step,standardize=standardize)
  ##find optimal lambda value that minimizes test MSE
  lambda_theta <- cv_model$lambda.min
  theta_model = glmnet(d_x_1, Z_2step,lambda=lambda_theta,standardize=standardize)
  theta_hat = theta_model$beta
  se = length(theta_hat@i)
  # step 3
  #d_x_1_stand = sweep(d_x_1, 2, colMeans(d_x_1))
  d_x_1_stand = d_x_1
  cv_model <- cv.glmnet(d_x_1_stand[,-1], d_x_1_stand[,1],standardize=standardize)
  lambda_cv <- cv_model$lambda.1se
  dx1_model = glmnet(d_x_1_stand[,-1], d_x_1_stand[,1], lambda=lambda_cv,standardize=standardize)
  gamma_hat = dx1_model$beta
  nonzero_index = which(gamma_hat!=0)
  #print(nonzero_index)
  #nonzero_index = gamma_hat@i+2
  sm = length(nonzero_index)
  ##########fit on second half###################
  # step 1
  if(sq==0){
    epsilon_hat = y_2 - beta0
    Z_2step = beta0 + ifelse(epsilon_hat<0,epsilon_hat,0)/alpha
  } else {
    dx2_short = matrix(d_x_2[,beta_hat!=0],ncol=sq)
    quan_model <- conquer(dx2_short,y_2,tau=alpha)
    beta_hat = quan_model$coeff[-1]
    beta0 = quan_model$coeff[1]
    #epsilon_hat = y_2 - dx2_short%*%beta_hat - beta0
    #epsilon_minus = ifelse(epsilon_hat<0,epsilon_hat,0)
    Z_2step = apply(cbind(y_2,dx2_short),MARGIN = 1,FUN = Z_beta,beta=beta_hat,alpha=alpha,beta_0=beta0)
  }
  # step 2
  if(se==0){
    res3 = Z_2step - mean(Z_2step)
  } else {
    dx2_short_theta = matrix(d_x_2[,theta_hat@i+1],ncol=se)
    theta_model = lm(Z_2step~dx2_short_theta)
    theta_hat = theta_model$coefficients[-1]
    res3 = theta_model$residuals
  }
  # step 3
  #d_x_2_stand = sweep(d_x_2, 2, colMeans(d_x_2))
  d_x_2_stand = d_x_2
  if(sm==0){
    w_hat = d_x_2_stand[,1] - mean(d_x_2_stand[,1])
  } else {
    lm_model <- lm(d_x_2_stand[,1]~d_x_2_stand[,nonzero_index+1])
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
  cv_model <- cv.glmnet(d_x_2, Z_2step,standardize=standardize)
  ##find optimal lambda value that minimizes test MSE
  lambda_theta <- cv_model$lambda.min
  theta_model = glmnet(d_x_2, Z_2step,lambda=lambda_theta,standardize=standardize)
  theta_hat = theta_model$beta
  se = length(theta_hat@i)
  # step 3
  cv_model <- cv.glmnet(d_x_2_stand[,-1], d_x_2_stand[,1],standardize=standardize)
  lambda_cv <- cv_model$lambda.1se
  dx2_model = glmnet(d_x_2_stand[,-1], d_x_2_stand[,1], lambda=lambda_cv,standardize=standardize)
  gamma_hat = dx2_model$beta
  nonzero_index = which(gamma_hat!=0)
  #print(nonzero_index)
  #nonzero_index = gamma_hat@i+2
  sm = length(nonzero_index)
  ##########fit on second half###################
  # step 1
  if(sq==0){
    epsilon_hat = y_1 - beta0
    Z_2step = beta0 + ifelse(epsilon_hat<0,epsilon_hat,0)/alpha
  } else {
    dx1_short = matrix(d_x_1[,beta_hat!=0],ncol=sq)
    quan_model <- conquer(dx1_short,y_1,tau=alpha)
    beta_hat = quan_model$coeff[-1]
    beta0 = quan_model$coeff[1]
    #epsilon_hat = y_1 - dx1_short%*%beta_hat - beta0
    #epsilon_minus = ifelse(epsilon_hat<0,epsilon_hat,0)
    Z_2step = apply(cbind(y_1,dx1_short),MARGIN = 1,FUN = Z_beta,beta=beta_hat,alpha=alpha,beta_0=beta0)
  }
  # step 2
  if(se==0){
    res3 = Z_2step - mean(Z_2step)
  } else {
    dx1_short_theta = matrix(d_x_1[,theta_hat@i+1],ncol=se)
    theta_model = lm(Z_2step~dx1_short_theta)
    theta_hat = theta_model$coefficients[-1]
    res3 = theta_model$residuals
  }
  # step 3
  if(sm==0){
    w_hat = d_x_1_stand[,1] - mean(d_x_1_stand[,1])
  } else {
    lm_model <- lm(d_x_1_stand[,1]~d_x_1_stand[,nonzero_index+1])
    w_hat = lm_model$residuals}
  
  #results
  s1 = nrow(d_x_1)
  sigma_w_hat_square_2 = sum(w_hat^2)/(s1-sm-1)
  sigma_s_hat_2 = alpha^2*t(w_hat^2)%*%((res3)^2)/(s1-sm-sq-se-3)
  sigma_w_hat_square = (sigma_w_hat_square_1+sigma_w_hat_square_2)/2
  sigma_s_hat = (sigma_s_hat_1+sigma_s_hat_2)/2
  #CI_level = 1-(1-conf.level)/2
  var = sigma_s_hat/sigma_w_hat_square^2
  return(var)
}

# supp functions

### compute adjusted response for the second step
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


self_tuning <- function(x,tau, nlambda,standardize=T) {
  if (standardize) {
    x <- scale(x,center = T, scale = T)
  }
  
  lambda_sim <- numeric(nlambda)
  for (b in 1:nlambda) {
    u <- runif(nrow(x), 0, 1) <= tau
    lambda_sim[b] <- max(abs(crossprod(x, (tau - u))))
  }
  
  return(lambda_sim / n)
}

quantile_bic = function(x,y,tau, lambda_seq=numeric(), nlambda=100,
                        h=max(0.05,sqrt(alpha*(1-alpha))*(log(p)/n)^(1/4)),
                        Cn=max(2, log(log(n))),standardize=T){
  n = nrow(x)
  p = ncol(x)
  if (length(lambda_seq) == 0) {
    lam_max <- max(self_tuning(x, tau, nlambda, standardize=standardize))
    lambda_seq <- seq(0.25 * lam_max, lam_max, length.out=nlambda)
  } else {
    nlambda <- length(lambda_seq)
  }
  check_sum <- function(x, tau) {
    sum(ifelse(x >= 0, tau * x, (tau - 1) * x))
  }
  model_all = conquer.reg(x,y,lambda=lambda_seq,tau=tau,h=h)
  beta_hat = model_all$coeff
  BIC <- log(apply(matrix(rep(y,nlambda),ncol=nlambda)-matrix(rep(beta_hat[1,],each=n),ncol=nlambda)-x%*%beta_hat[-1,], 2, check_sum, tau=tau)) + apply(beta_hat,2,function(x){sum(x!=0)}) * log(p) * Cn / n
  # plot(lambda_seq,BIC)
  bic_select <- which.min(BIC)
  lambda_select <- lambda_seq[bic_select]
  model_select <- conquer.reg(x,y,lambda=lambda_select,tau=tau,h=h)
  return(list(lambda_select=lambda_select,model_select=model_select))
}

LS_bic = function(x,z,lambda_seq=numeric(), nlambda=100,
                  Cn=max(2, log(log(n))),epsilon=.0001,standardize=T){
  n = nrow(x)
  p = ncol(x)
  model_all = glmnet(x,z,nlambda = nlambda, standardize=standardize)
  theta_hat = model_all$beta
  theta0 = model_all$a0
  lambda_seq = model_all$lambda
  nlambda = length(lambda_seq)
  BIC <- log(apply((matrix(rep(z,nlambda),ncol=nlambda)-matrix(rep(theta0,each=n),ncol=nlambda)-x%*%theta_hat)^2, 2, sum)/2/n) + apply(theta_hat,2,function(x){sum(x!=0)}) * log(p) * Cn / n
  bic_select <- which.min(BIC)
  lambda_select <- lambda_seq[bic_select]
  model_select <- glmnet(x,z,lambda=lambda_select,standardize=standardize)
  return(list(lambda_select=lambda_select,model_select=model_select))
}


