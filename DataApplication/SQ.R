dnames = '/Users/shushuz/Desktop/ES_highD/src/'
source(paste0(dnames,'highD_2step.R'))
library(haven)
library(quantreg)
library(xtable)
library(survey)

######### 1. High-dim ES Two Step Regression ######
alpha = 0.20
design_matrix = read.csv(paste0(dnames,'DataApplication/design_matrix_new.csv'))
y = -design_matrix[,"cotinine"]
x = subset(design_matrix,select = -c(cotinine))#,overnight2,mental2
x = x[,-grep("HSQ*",colnames(x))] 
x = as.matrix(x)
p = ncol(x)
n = nrow(x)

#####ES Estimation#####
####### Since we need intermediate values (i.e., beta_hat), we write the two-step procedure again here without using the "highdim_2step" function.
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
cbind(rownames(theta_hat)[which(theta_hat!=0)],theta_hat[which(theta_hat!=0)])


######## (1) Asian v.s. White ######
#x_stand = x - colMeans(x)
x_stand = sweep(x, 2, colMeans(x))
d = x_stand[,1]
x_tilde = x_stand[,-c(1,2,3)]
### CV ###
cv_model <- cv.glmnet(x_tilde, d,intercept=FALSE,standardize=FALSE)
plot(cv_model)
lambda_cv <- cv_model$lambda.1se
gamma_model = glmnet(x_tilde, d, lambda=lambda_cv,intercept=FALSE,standardize=FALSE)
gamma_hat = gamma_model$beta
resid2 = d - predict(gamma_model,newx=x_tilde)
theta_1_debias = theta_hat[1] + t(resid1)%*%resid2/x[,1]%*%resid2
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
#x_stand = sweep(x, 2, colMeans(x))
d = x_stand[,2]
#x_tilde = x_stand[,-c(1,2,3)]
### CV ###
cv_model <- cv.glmnet(x_tilde, d,lambda=lambdapath, standardize=FALSE,intercept=FALSE)
plot(cv_model)
lambda_cv <- cv_model$lambda.1se
gamma_model = glmnet(x_tilde, d, lambda=lambda_cv,intercept=FALSE,standardize=FALSE)
gamma_hat = gamma_model$beta
resid2 = d - predict(gamma_model,newx=x_tilde)
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
resid2 = d - predict(gamma_model,newx=x_tilde)
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

####Summary of Results
print(theta_1_debias);print(LB1);print(UB1)
print(theta_2_debias);print(LB2);print(UB2)
print(theta_3_debias);print(LB3);print(UB3)

#### 2. compare with low-dim estimation ####
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

###### 3. compare with LS lasso#####
### desparsified lasso
library(desla)
H = c(1,2,3)
despar = desla(x,y,H)#,alphas=0.01)
despar$intervals
