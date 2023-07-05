dnames = '/Users/shushuz/Desktop/ES_highD/src/'
source(paste0(dnames,'highD_2step.R'))
library(mvtnorm)
library(boot)
library(MultiRNG)

runsim_onepara = function(n,p,seed,alpha,s){
  
  ### generate 
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

################ two-step oracle method: only fit nonzero entries (low-dim) #######
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

########## two-step high-dim ES Estimation####################
  res_est = highdim_2step(x,y,alpha)
  theta0_hat = res_est$theta0_hat
  theta_hat = res_est$theta_hat
  resid1 = res_est$resid
  lambda_beta = res_est$lambda_beta
  lambda_theta = res_est$lambda_theta
  nonzero_index = which(theta_hat!=0)
  
  ####### refit model ####################
  Z_2step = apply(cbind(y,x),MARGIN = 1,FUN = Z_beta,beta=beta_star,alpha=alpha,beta_0=0)
  if(length(nonzero_index)==0){theta_hat_refit_long=theta_hat
  } else {
    Z_2step = apply(cbind(y,x),MARGIN = 1,FUN = Z_beta,beta=beta_star,alpha=alpha,beta_0=0)
    refit_model = lm(Z_2step~x[,nonzero_index])
    theta_hat_refit = refit_model$coefficients[-1]
    theta_hat_refit_long = rep(0,p)
    theta_hat_refit_long[nonzero_index] = theta_hat_refit}
  
  ##estimation outcome
  est_error = norm(theta_hat[1:s]-theta_star[1:s],type="2")/norm(theta_star,type="2")
  est_error_FP = norm(theta_hat[(s+1):p],type="2")/norm(theta_star,type="2")
  TPR = sum(nonzero_index<=s)/s
  FPR = sum(nonzero_index>s)/(p-s)
  est_error_refit = norm(theta_hat_refit_long[1:s]-theta_star[1:s],type="2")/norm(theta_star,type="2")
  est_error_refit_FP = norm(theta_hat_refit_long[(s+1):p],type="2")/norm(theta_star,type="2")
  TPR_refit = sum(which(theta_hat_refit_long!=0)<=s)/s
  FPR_refit = sum(which(theta_hat_refit_long!=0)>s)/(p-s)
  
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
  
  
  ###############two-step Bootstrap##########
highdim_boot = function(data,index){
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
theta0_boot = boot(cbind(y,x),highdim_boot,R=100)
#theta0_hat = theta0_boot$t0
theta0_hat_boot = mean(theta0_boot$t)

#######################Inference on one parameter of interest####################
res_inf = highdim_inf(x,y,alpha,col = 1,res_est=res_est)
theta_0_debias = res_inf$theta_debias
CI = res_inf$Conf.int
theta_0_star = theta_star[1]
res = (CI[1]<theta_0_star && theta_0_star<CI[2])
inference <- data.frame(n=n,p=p,seed=seed,alpha=alpha,s=s,results=res,theta0_star=theta_0_star,theta0_hat=theta0_hat,theta0_debias=theta_0_debias[1],theta0_boot=theta0_hat_boot,theta0_oracle=theta0_hat_lowdim)
write.table(inference,file=paste("onepara_results_ind/n",n,"p",p,"tau",alpha,"seed",seed,"s",s,".csv",sep=""),sep=",",row.names=FALSE)
return(0)
}
