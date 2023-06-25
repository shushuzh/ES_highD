# estimate theta
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

runsim_onepara = function(n,p,seed,alpha,s){
  set.seed(seed)
  
  ES_epsilon = -dnorm(qnorm(alpha))/alpha
  # generate sigma
 #sigma = matrix(c(1:p^2),nrow = p)
 # for (i in 1:p){
 #   for (j in 1:p){
 #     sigma[i,j] = (.6)^{abs(i-j)}
 #  }
 # }
  #generate x
 # x = abs(rmvnorm(n,rep(0,p),sigma))
  x = abs(rmvnorm(n,rep(0,p),diag(p)))
  #generate response
  gamma_star = c(rep(2,ceiling(s/2)),rep(1,s-ceiling(s/2)),rep(0,p-s))
  eta_star = c(rep(1/3,ceiling(s/2)),rep(0,p-ceiling(s/2)))
  beta_star = gamma_star + qnorm(alpha)*eta_star
  theta_star = gamma_star + ES_epsilon*eta_star
  epsilon = rnorm(n,0,1)
  y = x%*%gamma_star + diag(epsilon)%*%x%*%eta_star

################ compute oracle method: only fit nonzero entries (low-dim) #######
x_short = x[,1:s]
beta_star_short = beta_star[1:s]
theta_star_short = theta_star[1:s]
    # estimate theta
Z_short = apply(cbind(y,x_short),MARGIN = 1,FUN = Z_beta,beta=beta_star_short,alpha=alpha,beta_0=0)
lm_model_short = lm(Z_short~x_short)
theta_hat_short = lm_model_short$coefficients[-1]
theta0_hat_lowdim=theta_hat_short[1]
# outcome
est_error_short = norm(theta_hat_short-theta_star_short,type="2")/norm(theta_star,type="2")

####### true parameters ######
x_stand = sweep(x, 2, colMeans(x))
w_hat = x_stand[,1]
sigma_w_hat_true = sum(w_hat^2)/(n-1)
#method1
epsilon_hat = y - x%*%beta_star
epsilon_minus = ifelse(epsilon_hat<0,epsilon_hat,0)
sigma_s_hat_square_empirical = t(w_hat^2)%*%((epsilon_minus-alpha*x%*%(theta_star-beta_star))^2)/n
#sigma_s_hat_square_empirical = t(w_hat^2)%*%((epsilon_minus-mean(epsilon_minus))^2)/(n-1)
#method2
Z_2step = apply(cbind(y,x),MARGIN = 1,FUN = Z_beta,beta=beta_star,alpha=alpha,beta_0=0)
res4 = Z_2step - x%*%theta_star
sigma_s_hat_true = alpha^2*t(w_hat^2)%*%((res4)^2)/n
wid_true = sqrt(sigma_s_hat_true)*qnorm(0.975)/sqrt(n)/alpha/sigma_w_hat_true
wid_true_emp = sqrt(sigma_s_hat_square_empirical)*qnorm(0.975)/sqrt(n)/alpha/sigma_w_hat_true

############################high-dim Estimation####################
  #step 1: fit quantile regression for estimating beta
  quan_model = conquer.cv.reg(x,y,
                              #lambda=0.03,
                              tau=alpha,
                              h=max(0.05,sqrt(alpha*(1-alpha))*(log(p)/n)^(1/4)))
  lambda_beta = quan_model$lambda
  beta_hat = quan_model$coeff[-1]
  beta0 = quan_model$coeff[1]
  #step 2: fit ES regression for estimating theta
  Z_2step = apply(cbind(y,x),MARGIN = 1,FUN = Z_beta,beta=beta_hat,alpha=alpha,beta_0=beta0)
  ## CV: perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(x, Z_2step,standardize=FALSE)
  ##find optimal lambda value that minimizes test MSE
  lambda_theta <- cv_model$lambda.min
  theta_model = glmnet(x, Z_2step,lambda=lambda_theta,standardize=FALSE)
  theta_hat = theta_model$beta
  theta0_hat = theta_hat[1]
  resid1 = Z_2step - predict(theta_model,newx=x)
  nonzero_index = which(theta_hat!=0)
  ####### refit model ####################
  if(length(nonzero_index)==0){theta_hat_refit_long=theta_hat
  } else {refit_model = lm(Z_2step~x[,nonzero_index])
  theta_hat_refit = refit_model$coefficients[-1]
  theta_hat_refit_long = rep(0,p)
  theta_hat_refit_long[nonzero_index] = theta_hat_refit}
  ##outcome
  est_error = norm(theta_hat[1:s]-theta_star[1:s],type="2")/norm(theta_star,type="2")
  est_error_FP = norm(theta_hat[(s+1):p],type="2")/norm(theta_star,type="2")
  TPR = sum(nonzero_index<=s)/s
  FPR = sum(nonzero_index>s)/(p-s)
  est_error_refit = norm(theta_hat_refit_long[1:s]-theta_star[1:s],type="2")/norm(theta_star,type="2")
  est_error_refit_FP = norm(theta_hat_refit_long[(s+1):p],type="2")/norm(theta_star,type="2")
  TPR_refit = sum(which(theta_hat_refit_long!=0)<=s)/s
  FPR_refit = sum(which(theta_hat_refit_long!=0)>s)/(p-s)
  
  ###############Bootstrap##########
highdim = function(data,index){
 y = data[index,1]
 x = data[index,-1]  
 #fit quantile regression for estimating beta
  quan_model = conquer.reg(x,y,
                              lambda=lambda_beta,
                              tau=alpha,
                              h=max(0.05,sqrt(alpha*(1-alpha))*(log(p)/n)^(1/4)))
  beta_hat = quan_model$coeff[-1]
  beta0 = quan_model$coeff[1]

  ################## two step #########################
  Z_2step = apply(cbind(y,x),MARGIN = 1,FUN = Z_beta,beta=beta_hat,alpha=alpha,beta_0=beta0)
  ######### CV #########
  theta_model = glmnet(x, Z_2step,lambda=lambda_theta)
  theta_hat = theta_model$beta
  return (theta_hat[1])
}
theta0_boot = boot(cbind(y,x),highdim,R=100)
#theta0_hat = theta0_boot$t0
theta0_hat_boot = mean(theta0_boot$t)

#######################Inference on one parameter of interest####################
x_stand = sweep(x, 2, colMeans(x))
d = x_stand[,1]
x_tilde = x_stand[,-1]
### CV ###
cv_model <- cv.glmnet(x_tilde, d,intercept=FALSE,standardize=FALSE)
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

theta_0_debias = theta0_hat + t(resid1)%*%resid2/x[,1]%*%resid2
theta_0_debias_refit = theta0_hat + t(resid1)%*%resid2_refit/x[,1]%*%resid2_refit
theta_0_star = theta_star[1]
s_theta_hat = length(which(theta_hat!=0))

#########################Variance Estimation using RCV################################
### split sample into two
sample_index = sample(seq_len(nrow(x)),size = ceiling(n/2))
d_x_1 =x[sample_index,]
d_x_2 =x[-sample_index,]
y_1 = y[sample_index]
y_2 = y[-sample_index]
####d~x#######
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
### CV ###
cv_model <- cv.glmnet(d_x_1_stand[,-1], d_x_1_stand[,1],intercept=FALSE,standardize=FALSE)
lambda_cv <- cv_model$lambda.1se
dx1_model = glmnet(d_x_1_stand[,-1], d_x_1_stand[,1], lambda=lambda_cv,intercept=FALSE,standardize=FALSE)
gamma_hat = dx1_model$beta
nonzero_index = which(gamma_hat!=0)
print(nonzero_index)
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
sigma_w_hat_square_1 = sum(w_hat^2)/(s2-sm-1) #df: number of nonzero covariates + 1 intercept(standardization)
sigma_s_hat_square_1_empirical = t(w_hat^2)%*%(epsilon_minus-mean(epsilon_minus))^2/(s2-sm-sq-2)
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
### CV ###
cv_model <- cv.glmnet(d_x_2_stand[,-1], d_x_2_stand[,1],intercept=FALSE,standardize=FALSE)
lambda_cv <- cv_model$lambda.1se
dx2_model = glmnet(d_x_2_stand[,-1], d_x_2_stand[,1], lambda=lambda_cv,intercept=FALSE,standardize=FALSE)
gamma_hat = dx2_model$beta
nonzero_index = which(gamma_hat!=0)
print(nonzero_index)
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

sigma_s_hat_square_2_empirical = t(w_hat^2)%*%(epsilon_minus-mean(epsilon_minus))^2/(s1-sm-sq-2)
sigma_s_hat_2 = alpha^2*t(w_hat^2)%*%((res3)^2)/(s1-sm-sq-se-3)

sigma_w_hat_square = (sigma_w_hat_square_1+sigma_w_hat_square_2)/2
sigma_s_hat_square_empirical = (sigma_s_hat_square_1_empirical+sigma_s_hat_square_2_empirical)/2
sigma_s_hat = (sigma_s_hat_1+sigma_s_hat_2)/2

####### saving results ########
### estimation
estimation <-data.frame(Method=c("two-step","two-step+refitted","two-step oracle"),
                        alpha=rep(alpha,3),
                        n=rep(n,3),
                        p=rep(p,3),
                        seed=rep(seed,3),
                        s=rep(s,3),
                        est_error = c(est_error,est_error_refit,est_error_short),
                        est_error_FP = c(est_error_FP,est_error_refit_FP,0),
                        TPR = c(TPR,TPR_refit,NA),
                        FPR = c(FPR,FPR_refit,NA))
write.table(estimation,file=paste("results_ind/n",n,"p",p,"tau",alpha,"seed",seed,"s",s,".csv",sep=""),sep=",",
            row.names=FALSE)
########Inference: CI
# wid1 and wid2 are two slightly different formulations to estimate sigma_s_hat. We present res1 in the paper. 
wid1 = sqrt(sigma_s_hat_square_empirical)*qnorm(0.975)/sqrt(n)/alpha/sigma_w_hat_square
res1 = (theta_0_debias-wid1<theta_0_star && theta_0_star<theta_0_debias+wid1)
wid2 = sqrt(sigma_s_hat)*qnorm(0.975)/sqrt(n)/alpha/sigma_w_hat_square
res2 = (theta_0_debias-wid2<theta_0_star && theta_0_star<theta_0_debias+wid2)
# res 3 and 4 are baselines
res3 = (theta_0_debias-wid_true_emp<theta_0_star && theta_0_star<theta_0_debias+wid_true_emp)
res4 = (theta_0_debias-wid_true<theta_0_star && theta_0_star<theta_0_debias+wid_true)
inference <- data.frame(n=n,p=p,seed=seed,alpha=alpha,s=s,results1=res1,results2=res2,results3=res3,results4=res4,theta0_star=theta_0_star,theta0_hat=theta0_hat,theta0_debias=theta_0_debias[1],theta0_debias_refit=theta_0_debias_refit[1],theta0_boot=theta0_hat_boot,theta0_oracle=theta0_hat_lowdim,width1=wid1,width2=wid2,width3=wid_true_emp,width4=wid_true,sigma_w_hat_square=sigma_w_hat_square,sigma_s_hat_square_empirical=sigma_s_hat_square_empirical,sigma_s_hat=sigma_s_hat,sigma_w_hat_true=sigma_w_hat_true,sigma_s_hat_true=sigma_s_hat_true)
write.table(inference,file=paste("onepara_results_ind/n",n,"p",p,"tau",alpha,"seed",seed,"s",s,".csv",sep=""),sep=",",row.names=FALSE)
return(0)
}
