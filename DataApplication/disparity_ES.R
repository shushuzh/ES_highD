dnames = '/Users/shushuzhang/Dropbox (University of Michigan)/ES_HighD/Empirical_studies/ES_highD/'
source(paste0(dnames,'src/highD_2step.R'))
library(haven)
library(quantreg)
library(xtable)
library(survey)
library(conquer)
library(ggplot2)
library(ggpubr)

######### 1. High-dim ES Two Step Regression ######
design_matrix = read.csv(paste0(dnames,'DataApplication/design_matrix_new.csv'))
y = -design_matrix[,"cotinine"]
x = subset(design_matrix,select = -c(cotinine))#,overnight2,mental2
x = x[,-grep("HSQ*",colnames(x))] 
x = as.matrix(x)
p = ncol(x)
n = nrow(x)


res = function(alpha){
res_est = highdim_2step(x,y,alpha)

## Asian v.s. White ##
res_inf1 = highdim_inf(x,y,alpha,col = 1,res_est=res_est)
theta_1_debias = res_inf1$theta_debias
CI1 = res_inf1$Conf.int

## Black v.s. White ##
res_inf2 = highdim_inf(x,y,alpha,col = 2,res_est=res_est)
theta_2_debias = res_inf2$theta_debias
CI2 = res_inf2$Conf.int

## Hispanic v.s. White ##
res_inf3 = highdim_inf(x,y,alpha,col = 3,res_est=res_est)
theta_3_debias = res_inf3$theta_debias
CI3 = res_inf3$Conf.int

####Summary of Results
return(list(mean=c(theta_1_debias,theta_2_debias,theta_3_debias),
            lower=c(CI1[1],CI2[1],CI3[1]),
            upper=c(CI1[2],CI2[2],CI3[2]),
            #theta0_hat = res_est$theta0_hat,
            theta_hat = res_est$theta_hat[1:3],
            #beta0_hat = res_est$beta0_hat,
            beta_hat = res_est$beta_hat[1:3]
            ))
}

# res_.05 = res(0.05)
# res_.1 = res(0.1)
# res_.15 = res(0.15)
# res_.2 = res(0.2)
#write.table(res_.2,file=paste0(dnames,'DataApplication/results.csv'))

results = read.csv(paste0(dnames,'DataApplication/results.csv'))
# results = - results
results$ci_group = rep(c("Asian", "Black", "Hispanic"), times = 6)
results$x = rep(c(0.95,0.9,0.85,0.8,0.75,0.7), each = 3)
## plot ##
#pdf(file=paste0(dnames,'DataApplication/CI_results.pdf'),width=8, height=4)
p1 <- ggplot(results, aes(x = x, y = mean, group = ci_group, color = ci_group)) +
  geom_point(aes(shape=ci_group),position = position_dodge(width=0.01),size=4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.01,size=1, position = position_dodge(width=0.01)) +
  scale_color_manual(values = c("Asian" = "red", "Black" = "black", "Hispanic" = "grey")) +
  scale_shape_manual(values = c(18, 17, 16)) +
  theme_minimal() +
  geom_abline(slope=0,intercept=0,linetype="dashed") + 
  labs(x = "Quantile level", y = "Disparity in Cotinine (ng/ml)",
       color = "Race",shape="Race")
#dev.off()

#### 2. compare with low-dim estimation ####
alpha = 0.2
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

###### 3. compare with mean lasso regression#####
### desparsified lasso ###
library(desla)
H = c(1,2,3)
despar = desla(x,y,H)#,alphas=0.01)
ci_lasso = -despar$intervals
ci_lasso = data.frame(ci_lasso)
ci_lasso$ci_group = c("Asian","Black","Hispanic")

####### Lasso plot #####
#pdf(file=paste0(dnames,'DataApplication/CI_Lasso.pdf'),width=8, height=8)
# Create a numeric variable for x-axis positioning
p2 <- ggplot(ci_lasso, aes(x = ci_group, y = bhat, color = ci_group)) +
  geom_point(aes(shape=ci_group),size=3) +
  geom_errorbar(aes(ymin = upper.0.05, ymax = lower.0.05), width = 0.05,size=1) +
  scale_color_manual(values = c("Asian" = "red", "Black" = "black", "Hispanic" = "grey")) +
  scale_shape_manual(values = c(18, 17, 16)) +
  theme_minimal() +
  geom_abline(slope=0,intercept=0,linetype="dashed") + 
  labs(x = "Race", y = "Disparity in Cotinine (ng/ml)",
       color = "Race",shape="Race") +
  theme(legend.position = "none") + 
  ylim(c(-214.1389,247.6286))
#dev.off()

##### Combine 
pdf(file=paste0(dnames,'DataApplication/CI_2.pdf'),width=16, height=8)
ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(3,1))
dev.off()

###### Lasso residual plot ######
resid = despar$residuals$init
category_labels <- c("White", "Asian", "Black", "Hispanic")
dummy_matrix = x[,1:4]
is_baseline <- rowSums(dummy_matrix) == 0
max_indices <- max.col(dummy_matrix, ties.method = "first") + 1
max_indices[is_baseline] <- 1
categorical_variable <- category_labels[max_indices]
categorical_variable <- factor(categorical_variable, levels = category_labels)
pdf(file=paste0(dnames,'DataApplication/Lasso_resid_plot.pdf'),width=8, height=4)
plot(categorical_variable,resid,xlab="Race",ylab="Residual")
dev.off()




