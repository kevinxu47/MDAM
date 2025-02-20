library(mvtnorm)
library(tidyverse)
library(survey)
library(mice)
library(nleqslv)

plogit = function(X){
  return(1/(1+exp(-X)))
}

logit = function(x){
  return(log(x)-log(1-x))
}

logit_inv = function(x){
  return(exp(x)/(1+exp(x)))
}


# Function of drawing samples from normal distribution
draw_samples_X <- function(T_X, V_X, dot_products_X_U0, max_sum_w_U1, n_simulations) {
  valid_samples <- c()
  while (length(valid_samples) < n_simulations) {
    sample <- rnorm(1, mean = T_X, sd = sqrt(V_X))
    if (sample > dot_products_X_U0 && sample < (dot_products_X_U0 + max_sum_w_U1)) {
      valid_samples <- c(valid_samples, sample)
    }
  }
  return(valid_samples)
}



# Define the function to solve
solve_beta0 <- function(beta0, beta1, A, B, C) {
  logit_inverse <- function(x) {
    return(1 / (1 + exp(-x)))
  }
  f <- logit_inverse(beta0 + beta1) * A + logit_inverse(beta0) * B - C
  return(f)
}

SAMPLE <- readRDS("SAMPLE.rds")
population <- readRDS("population.rds")


# MDAM-fix imputation
MDAM_fix = function(sub_dat, T_X1 = sum(population$X1), T_X2 = sum(population$X2), L = 10){
  
  n_unit0 = sum(sub_dat$U==0) # unit respondents
  n_unit1 = sum(sub_dat$U==1) # unit nonrespondents
  N = nrow(population)
  
  # adjust unit nonrespondents weights
  sub_dat$origin_W = sub_dat$W
  sub_dat$W = sub_dat$W*10
  sub_dat$pi = 1 / sub_dat$W
  
  # Calculate sum of weights
  sum_W = sum(sub_dat$W)
  # Weighting Factor WF
  WF = sum(sub_dat$W) / N
  
  # Delete data for Unit Nonresponse
  sub_dat[sub_dat$U==1, c('X1', 'X2', 'Y1', 'Y2', 'Y3', 'Y4','Rx1','Ry1','Ry2','Ry4')] = NA
  # Delete data for Item Nonresponse
  sub_dat[sub_dat$U==0 & sub_dat$Rx1==1,]$X1 = NA
  sub_dat[sub_dat$U==0 & sub_dat$Ry1==1,]$Y1 = NA
  sub_dat[sub_dat$U==0 & sub_dat$Ry2==1,]$Y2 = NA
  sub_dat[sub_dat$U==0 & sub_dat$Ry4==1,]$Y4 = NA
  
  subdat.completed1 = mice(sub_dat[,c("X1","X2","Y1","Y2","Y3","Y4",'Rx1','Ry1','Ry2','Ry4')], m=1, maxit = 5, seed=1)
  
  adj_dat = complete(subdat.completed1) %>% mutate(U=sub_dat$U, W=sub_dat$W, pi=sub_dat$pi, origin_W=sub_dat$origin_W)
  
  V_X1_pop_HT = sum((adj_dat$X1/adj_dat$pi)^2*(1-adj_dat$pi))
  V_X2_pop_HT = sum((adj_dat$X2/adj_dat$pi)^2*(1-adj_dat$pi))
  
  impdat = mice(sub_dat[,c("X1","X2","Y1","Y2","Y3","Y4",'Rx1','Ry1','Ry2','Ry4')],
                m=L, maxit = 5, seed=1)
  
  IMP_DAT = list()
  
  for(i in 1:L){
    adj_dat = complete(impdat, i) %>% 
      mutate(U=sub_dat$U, W=sub_dat$W, pi=sub_dat$pi, origin_W=sub_dat$origin_W)
    
    # Make Unit Nonresponse Missing values for imputation
    adj_dat[adj_dat$U==1, c('X1', 'X2', 'Y1', 'Y2', 'Y3', 'Y4', 'Rx1','Ry1','Ry2','Ry4')] = NA
    
    ###### Impute X1 (unit nonresponse only) ######
    
    dot_products_X1_U0 = sum( adj_dat[adj_dat$U == 0, ]$W * adj_dat[adj_dat$U == 0, ]$X1 )
    dot_products_X2_U0 = sum( adj_dat[adj_dat$U == 0, ]$W * adj_dat[adj_dat$U == 0, ]$X2 )
    sum_w_U1 = sum(adj_dat[adj_dat$U == 1, ]$W)
    
    # Generate T_X1_hat
    T_X1_hat = draw_samples_X(T_X1, V_X1_pop_HT, dot_products_X1_U0, sum_w_U1, 1)
      
    # Calculate X1_prob
    X1_prob = (T_X1_hat - sum(adj_dat[adj_dat$U == 0 & adj_dat$X1 == 1, ]$W)) / sum(adj_dat[adj_dat$U == 1, ]$W)
      
    # Check if X1_prob is within the required range
    if (X1_prob > 1) {
      X1_prob = 1
    }
      
    if (X1_prob < 0) {
      X1_prob = 0
    }
    
    adj_dat$X1[which(adj_dat$U == 1)] = rbinom(n_unit1, 1, prob = X1_prob)
    
    ###### Impute X2 (unit nonresponse only) ######
    n_U_1_X1_1 = nrow(adj_dat[which(adj_dat$U == 1 & adj_dat$X1 == 1),])
    n_U_1_X1_0 = nrow(adj_dat[which(adj_dat$U == 1 & adj_dat$X1 == 0),])    
    
    # Draw T_X2_hat
    T_X2_hat = draw_samples_X(T_X2, V_X2_pop_HT, dot_products_X2_U0, sum_w_U1, 1)
    p_X2 = (T_X2_hat - dot_products_X2_U0) / sum_w_U1
    
    # P(X2i = 1 | X1i = 1, Ui = 0)
    p_x2_1_given_x1_1 <- sum(adj_dat[adj_dat$U==0 & adj_dat$X1 == 1 & adj_dat$X2 == 1,]$W) /
      sum(adj_dat[adj_dat$U==0 & adj_dat$X1 == 1,]$W)
    # P(X2i = 1 | X1i = 0, Ui = 0)
    p_x2_1_given_x1_0 <- sum(adj_dat[adj_dat$U==0 & adj_dat$X1 == 0 & adj_dat$X2 == 1,]$W) / 
      sum(adj_dat[adj_dat$U==0 & adj_dat$X1 == 0,]$W)
    
    # P(X2 = 1 | X1 = 0, U = 1) is calculated by fixing log odd ratio from U = 0
    beta1_U_0 = logit(p_x2_1_given_x1_1) - logit(p_x2_1_given_x1_0)  
    beta0_U_0 = logit(p_x2_1_given_x1_1) - beta1_U_0
    
    beta1_U_1 = beta1_U_0    
    beta0_U_1 = nleqslv(x = 0, fn = solve_beta0, 
                        beta1 = beta1_U_1, 
                        A = sum(adj_dat[adj_dat$U==1 & adj_dat$X1 == 1,]$W), 
                        B = sum(adj_dat[adj_dat$U==1 & adj_dat$X1 == 0,]$W), 
                        C = T_X2_hat - dot_products_X2_U0)$x
    
    X2_prob_X1_1 = logit_inv(beta0_U_1 + beta1_U_1)
    X2_prob_X1_0 = logit_inv(beta0_U_1)
    
    if(X2_prob_X1_1 < 1 && X2_prob_X1_1 > 0){
      adj_dat[which(adj_dat$U == 1 & adj_dat$X1 == 1),]$X2 = rbinom(n_U_1_X1_1, 1, prob = X2_prob_X1_1)
    } else{
      if(X2_prob_X1_1 >= 1){
        adj_dat[which(adj_dat$U == 1 & adj_dat$X1 == 1),]$X2 = rep(1, n_U_1_X1_1)
      }
      
      if(X2_prob_X1_1 <= 0){
        adj_dat[which(adj_dat$U == 1 & adj_dat$X1 == 1),]$X2 = rep(0, n_U_1_X1_1)
      }
    }
    
    if(X2_prob_X1_0 < 1 && X2_prob_X1_0 > 0){
      adj_dat[which(adj_dat$U == 1 & adj_dat$X1 == 0),]$X2 = rbinom(n_U_1_X1_0, 1, prob = X2_prob_X1_0)
    } else{
      if(X2_prob_X1_0 >= 1){
        adj_dat[which(adj_dat$U == 1 & adj_dat$X1 == 0),]$X2 = rep(1, n_U_1_X1_0)
      }
      
      if(X2_prob_X1_0 <= 0){
        adj_dat[which(adj_dat$U == 1 & adj_dat$X1 == 0),]$X2 = rep(0, n_U_1_X1_0)
      }
    }
    
    ###### hot deck ###### 
    cell_00 = which(adj_dat$X1==0 & adj_dat$X2==0 & adj_dat$U==0)
    cell_01 = which(adj_dat$X1==0 & adj_dat$X2==1 & adj_dat$U==0)
    cell_10 = which(adj_dat$X1==1 & adj_dat$X2==0 & adj_dat$U==0)
    cell_11 = which(adj_dat$X1==1 & adj_dat$X2==1 & adj_dat$U==0)
    
    # Create a lookup table (all possible combinations of x1,x2)
    imp_lookup = list("00" = cell_00, "01" = cell_01, "10" = cell_10, "11" = cell_11)
    
    # impute unit nonrespondents for Y1,Y2,Y3,Y4
    rows_to_impute = which(adj_dat$U==1)
    for (k in rows_to_impute) {
      imp_key = paste0(adj_dat[k, "X1"], adj_dat[k, "X2"])
      adj_dat[k, c("Y1", "Y2", "Y3", "Y4")] = adj_dat[sample(imp_lookup[[imp_key]], 1), c("Y1", "Y2", "Y3", "Y4")]
    }
    IMP_DAT[[i]]=adj_dat
  }
  
  return(IMP_DAT)
}

# The following codes are copied from Yang and Reiter, 2025

RESULTS = list()
for(i in 1:length(SAMPLE)){
  RESULTS[[i]] = MDAM_fix(sub_dat=SAMPLE[[i]])
}


n_sim = length(SAMPLE)

T_est = T_premiss = T_SRS = matrix(NA, nrow=n_sim, ncol=6) 
cond_prob = cond_prob_premiss = matrix(NA, nrow=n_sim, ncol=8) 
joint_prob = joint_prob_premiss = matrix(NA, nrow=n_sim, ncol=4)
cond_prob_given12 = cond_prob_given12_premiss = matrix(NA, nrow=n_sim, ncol=8)


V_B_T = matrix(NA, nrow=n_sim, ncol=6) 
V_B_prob = matrix(NA, nrow=n_sim, ncol=8) 
V_B_given12 = matrix(NA, nrow=n_sim, ncol=8) 
V_B_joint = matrix(NA, nrow=n_sim, ncol=4) 


V_W_T = matrix(NA, nrow=n_sim, ncol=6) 
V_W_cond = matrix(NA, n_sim, 8) 
V_W_given12 = matrix(NA, n_sim, 8) 
V_W_joint = matrix(NA, n_sim, 4)


V_T = matrix(NA, nrow=n_sim, ncol=6)


N = nrow(population)
L = 10


for(i in 1:n_sim){
  sub_dat = SAMPLE[[i]]
  n = nrow(sub_dat)
  sub_dat$origin_W = sub_dat$W
  sub_dat$W = sub_dat$W*10
  sub_dat$pi = sub_dat$pi/10
  
  # Use original W and pi
  W = (sub_dat$origin_W)*10
  pi = 1 / W
  
  T_premiss[i,] = sapply(sub_dat[c("X1", "X2", "Y1", "Y2", "Y3", "Y4")], function(X) sum(X * sub_dat$W)) 
  
  premiss_dat = sub_dat %>% select(X1,X2,Y1,Y2,U) %>% 
    mutate(X1=as.factor(X1), X2=as.factor(X2), Y1=as.factor(Y1), Y2=as.factor(Y2)) %>%
    mutate(W = as.numeric(W), pi = as.numeric(pi), fpc = as.numeric(nrow(population)))
  mydesign_premiss = svydesign(id=~1, data = premiss_dat, weight = ~sub_dat$W, fpc=~fpc)
  
  cond_prob_premiss[i,1:2] = svyby(~as.factor(X2), ~as.factor(X1), mydesign_premiss, svymean)[1:2,2] # p(X2|X1)
  cond_prob_premiss[i,3:4] = svyby(~as.factor(X1), ~as.factor(X2), mydesign_premiss, svymean)[1:2,2] # p(X1|X2)
  cond_prob_premiss[i,5:6] = svyby(~as.factor(Y2), ~as.factor(Y1), mydesign_premiss, svymean)[1:2,2]  # p(Y2|Y1)
  cond_prob_premiss[i,7:8] = svyby(~as.factor(Y1), ~as.factor(Y2), mydesign_premiss, svymean)[1:2,2]  # p(Y1|Y2)
  
  joint_prob_premiss[i,1] = as.numeric(coef(svymean(~I(X1==0 & Y1==0), design=mydesign_premiss))[2]) # estimate for p(X1,Y1)
  joint_prob_premiss[i,2] = as.numeric(coef(svymean(~I(X1==1 & Y1==0), design=mydesign_premiss))[2]) # estimate for p(X1,Y1)
  joint_prob_premiss[i,3] = as.numeric(coef(svymean(~I(X1==0 & Y1==1), design=mydesign_premiss))[2]) # estimate for p(X1,Y1)
  joint_prob_premiss[i,4] = as.numeric(coef(svymean(~I(X1==1 & Y1==1), design=mydesign_premiss))[2]) # estimate for p(X1,Y1)
  
  cond_prob_given12_premiss[i, 1:4] = svyby(~as.factor(Y1), ~X1+X2, design=mydesign_premiss, svymean)[1:4,3]
  cond_prob_given12_premiss[i, 5:8] = svyby(~as.factor(Y2), ~X1+X2, design=mydesign_premiss, svymean)[1:4,3]
  
  T_EST = matrix(NA, nrow=L, ncol=6)  
  VW_T = matrix(NA, nrow=L, ncol=6)  
  COND_PROB = matrix(NA, nrow=L, ncol=8)
  VW_COND_PROB = matrix(NA, nrow=L, ncol=8)
  JOINT_PROB = matrix(NA, nrow=L, ncol=4)
  VW_JOINT_PROB = matrix(NA, nrow=L, ncol=4)
  GIVEN12 = matrix(NA, L, 8)
  VW_GIVEN12 = matrix(NA, L, 8)
  
  for(j in 1:L){
    IMP_dat = cbind(id = as.numeric(rownames(sub_dat)), 
                    RESULTS[[i]][[j]], fpc = as.numeric(nrow(population)))
    
    mydesign_imp = svydesign(id=~1, data = IMP_dat, weight = ~W) 
    T_EST[j,] = svytotal(~X1+X2+Y1+Y2+Y3+Y4, mydesign_imp)[1:6]
    VW_T[j, ] = SE(svytotal(~X1+X2+Y1+Y2+Y3+Y4, mydesign_imp))^2
    
    # cond_prob
    COND_PROB[j,1:2] = svyby(~as.factor(X2), ~as.factor(X1), mydesign_imp, svymean)[1:2,2] 
    COND_PROB[j,3:4] = svyby(~as.factor(X1), ~as.factor(X2), mydesign_imp, svymean)[1:2,2] 
    COND_PROB[j,5:6] = svyby(~as.factor(Y2), ~as.factor(Y1), mydesign_imp, svymean)[1:2,2]  
    COND_PROB[j,7:8] = svyby(~as.factor(Y1), ~as.factor(Y2), mydesign_imp, svymean)[1:2,2]  
    VW_COND_PROB[j,1:2] = (svyby(~as.factor(X2), ~as.factor(X1), mydesign_imp, svymean)[1:2,4])^2 
    VW_COND_PROB[j,3:4] = (svyby(~as.factor(X1), ~as.factor(X2), mydesign_imp, svymean)[1:2,4])^2 
    VW_COND_PROB[j,5:6] = (svyby(~as.factor(Y2), ~as.factor(Y1), mydesign_imp, svymean)[1:2,4])^2 
    VW_COND_PROB[j,7:8] = (svyby(~as.factor(Y1), ~as.factor(Y2), mydesign_imp, svymean)[1:2,4])^2 
    
    # joint_prob
    JOINT_PROB[j,1] = as.numeric(coef(svymean(~I(X1==0 & Y1==0), design=mydesign_imp))[2])
    JOINT_PROB[j,2] = as.numeric(coef(svymean(~I(X1==1 & Y1==0), design=mydesign_imp))[2])
    JOINT_PROB[j,3] = as.numeric(coef(svymean(~I(X1==0 & Y1==1), design=mydesign_imp))[2]) 
    JOINT_PROB[j,4] = as.numeric(coef(svymean(~I(X1==1 & Y1==1), design=mydesign_imp))[2]) 
    VW_JOINT_PROB[j,1] = vcov(svymean(~I(X1==0 & Y1==0), design=mydesign_imp))[1]
    VW_JOINT_PROB[j,2] = vcov(svymean(~I(X1==1 & Y1==0), design=mydesign_imp))[1]
    VW_JOINT_PROB[j,3] = vcov(svymean(~I(X1==0 & Y1==1), design=mydesign_imp))[1]
    VW_JOINT_PROB[j,4] = vcov(svymean(~I(X1==1 & Y1==1), design=mydesign_imp))[1]
    
    # cond prob given X1,X2
    GIVEN12[j,1:4] = svyby(~as.factor(Y1), ~X1+X2, design=mydesign_imp, svymean)[1:4,3]
    GIVEN12[j,5:8] = svyby(~as.factor(Y2), ~X1+X2, design=mydesign_imp, svymean)[1:4,3]
    VW_GIVEN12[j,1:4] = (svyby(~as.factor(Y1), ~X1+X2, design=mydesign_imp, svymean)[1:4, 5])^2
    VW_GIVEN12[j,5:8] = (svyby(~as.factor(Y2), ~X1+X2, design=mydesign_imp, svymean)[1:4, 5])^2
  }
  
  T_est[i,] = colMeans(T_EST)
  V_B_T[i,] = apply(T_EST, 2, var) * L / (L - 1) 
  V_W_T[i,] = colMeans(VW_T)
  
  cond_prob[i,] = colMeans(COND_PROB)
  V_B_prob[i,] = apply(COND_PROB, 2, var) * L / (L - 1) 
  V_W_cond[i,] = colMeans(VW_COND_PROB)
  
  joint_prob[i,] = colMeans(JOINT_PROB)
  V_B_joint[i,] = apply(JOINT_PROB, 2, var) * L / (L - 1) 
  V_W_joint[i,] = colMeans(VW_JOINT_PROB)
  
  cond_prob_given12[i,] = colMeans(GIVEN12)
  V_B_given12[i,] = apply(GIVEN12, 2, var) * L / (L - 1) 
  V_W_given12[i,] = colMeans(VW_GIVEN12)
}

colMeans(T_est)
colMeans(cond_prob)
colMeans(joint_prob)
colMeans(cond_prob_given12)

cond_prob_true = c(nrow(population %>% filter(X1==0,X2==0))/nrow(population %>% filter(X1==0)),
                   nrow(population %>% filter(X1==1,X2==0))/nrow(population %>% filter(X1==1)), 
                   nrow(population %>% filter(X1==0,X2==0))/nrow(population %>% filter(X2==0)), 
                   nrow(population %>% filter(X1==0,X2==1))/nrow(population %>% filter(X2==1)),
                   nrow(population %>% filter(Y1==0,Y2==0))/nrow(population %>% filter(Y1==0)), 
                   nrow(population %>% filter(Y1==1,Y2==0))/nrow(population %>% filter(Y1==1)), 
                   nrow(population %>% filter(Y1==0,Y2==0))/nrow(population %>% filter(Y2==0)), 
                   nrow(population %>% filter(Y1==0,Y2==1))/nrow(population %>% filter(Y2==1)) 
) 

joint_prob_true = c(nrow(population %>% filter(X1==0,Y1==0))/nrow(population), 
                    nrow(population %>% filter(X1==1,Y1==0))/nrow(population), 
                    nrow(population %>% filter(X1==0,Y1==1))/nrow(population), 
                    nrow(population %>% filter(X1==1,Y1==1))/nrow(population))

cond_prob_given12_true =
  c(nrow(population %>% filter(Y1==0,X1==0,X2==0))/nrow(population %>% filter(X1==0,X2==0)), 
    nrow(population %>% filter(Y1==0,X1==1,X2==0))/nrow(population %>% filter(X1==1,X2==0)), 
    nrow(population %>% filter(Y1==0,X1==0,X2==1))/nrow(population %>% filter(X1==0,X2==1)), 
    nrow(population %>% filter(Y1==0,X1==1,X2==1))/nrow(population %>% filter(X1==1,X2==1)), 
    nrow(population %>% filter(Y2==0,X1==0,X2==0))/nrow(population %>% filter(X1==0,X2==0)), 
    nrow(population %>% filter(Y2==0,X1==1,X2==0))/nrow(population %>% filter(X1==1,X2==0)), 
    nrow(population %>% filter(Y2==0,X1==0,X2==1))/nrow(population %>% filter(X1==0,X2==1)), 
    nrow(population %>% filter(Y2==0,X1==1,X2==1))/nrow(population %>% filter(X1==1,X2==1))  
  )


V_T = V_W_T + (1+1/L)*V_B_T
V_cond = V_W_cond + (1+1/L)*V_B_prob
V_given12 = V_W_given12 + (1+1/L)*V_B_given12
V_joint = V_W_joint+ (1+1/L)*V_B_joint

CI.T1 = 1 - sum(T_est[,1] - 2.009*sqrt(colMeans(V_T)[1]) >  sum(population$X1) | T_est[,1] + 2.009*sqrt((colMeans(V_T)[1]))  <  sum(population$X1))/n_sim
CI.T2 = 1 - sum(T_est[,2] - 2.009*sqrt(colMeans(V_T)[2]) >  sum(population$X2) | T_est[,2] + 2.009*sqrt(colMeans(V_T)[2]) <  sum(population$X2))/n_sim 
CI.T3 = 1 - sum(T_est[,3] - 2.009*sqrt(colMeans(V_T)[3]) >  sum(population$Y1) | T_est[,3] + 2.009*sqrt(colMeans(V_T)[3]) <  sum(population$Y1))/n_sim 
CI.T4 = 1 - sum(T_est[,4] - 2.009*sqrt(colMeans(V_T)[4]) >  sum(population$Y2) | T_est[,4] + 2.009*sqrt(colMeans(V_T)[4]) <  sum(population$Y2))/n_sim 
CI.T5 = 1 - sum(T_est[,5] - 2.009*sqrt(colMeans(V_T)[5]) >  sum(population$Y3) | T_est[,5] + 2.009*sqrt(colMeans(V_T)[5]) <  sum(population$Y3))/n_sim 
CI.T6 = 1 - sum(T_est[,6] - 2.009*sqrt(colMeans(V_T)[6]) >  sum(population$Y4) | T_est[,6] + 2.009*sqrt(colMeans(V_T)[6]) <  sum(population$Y4))/n_sim 

cond_prob_true1 = nrow(population %>% filter(X2==0, X1==0))/nrow(population %>% filter(X1==0))
CI.cond1 = 1 - sum(cond_prob[,1] - 2.009*sqrt(colMeans(V_cond)[1]) >  cond_prob_true1 | cond_prob[,1] + 2.009*sqrt(colMeans(V_cond)[1]) <  cond_prob_true1)/n_sim 

cond_prob_true2 = nrow(population %>% filter(X2==0, X1==1))/nrow(population %>% filter(X1==1))
CI.cond2 = 1 - sum(cond_prob[,2] - 2.009*sqrt(colMeans(V_cond)[2]) >  cond_prob_true2 | cond_prob[,2] + 2.009*sqrt(colMeans(V_cond)[2]) <  cond_prob_true2)/n_sim 

cond_prob_true3 = nrow(population %>% filter(X2==0, X1==0))/nrow(population %>% filter(X2==0))
CI.cond3 = 1 - sum(cond_prob[,3] - 2.009*sqrt(colMeans(V_cond)[3]) >  cond_prob_true3 | cond_prob[,3] + 2.009*sqrt(colMeans(V_cond)[3]) <  cond_prob_true3)/n_sim

cond_prob_true4 = nrow(population %>% filter(X2==1, X1==0))/nrow(population %>% filter(X2==1))
CI.cond4 = 1 - sum(cond_prob[,4] - 2.009*sqrt(colMeans(V_cond)[4]) >  cond_prob_true4 | cond_prob[,4] + 2.009*sqrt(colMeans(V_cond)[4]) <  cond_prob_true4)/n_sim 

cond_prob_true5 = nrow(population %>% filter(Y1==0, Y2==0))/nrow(population %>% filter(Y1==0))
CI.cond5 = 1 - sum(cond_prob[,5] - 2.009*sqrt(colMeans(V_cond)[5]) >  cond_prob_true5 | cond_prob[,5] + 2.009*sqrt(colMeans(V_cond)[5]) <  cond_prob_true5)/n_sim 

cond_prob_true6 = nrow(population %>% filter(Y1==1, Y2==0))/nrow(population %>% filter(Y1==1))
CI.cond6 = 1 - sum(cond_prob[,6] - 2.009*sqrt(colMeans(V_cond)[6]) >  cond_prob_true6 | cond_prob[,6] + 2.009*sqrt(colMeans(V_cond)[6]) <  cond_prob_true6)/n_sim 

cond_prob_true7 = nrow(population %>% filter(Y1==0, Y2==0))/nrow(population %>% filter(Y2==0))
CI.cond7 = 1 - sum(cond_prob[,7] - 2.009*sqrt(colMeans(V_cond)[7]) >  cond_prob_true7 | cond_prob[,7] + 2.009*sqrt(colMeans(V_cond)[7]) <  cond_prob_true7)/n_sim 

cond_prob_true8 = nrow(population %>% filter(Y1==0, Y2==1))/nrow(population %>% filter(Y2==1))
CI.cond8 = 1 - sum(cond_prob[,8] - 2.009*sqrt(colMeans(V_cond)[8]) >  cond_prob_true8 | cond_prob[,8] + 2.009*sqrt(colMeans(V_cond)[8]) <  cond_prob_true8)/n_sim 

CI_given12 = c()
for (k in 1:8){
  p_true = cond_prob_given12_true[k]
  var = colMeans(V_given12)[k]
  CI_given12[k] = 1 - sum(cond_prob_given12[,k] - 2.009*sqrt(var) > p_true | cond_prob_given12[,k] + 2.009*sqrt(var) < p_true)/n_sim
}

CI_joint = c()
joint_true = c(nrow(population%>%filter(X1==0, Y1==0))/N,
               nrow(population%>%filter(X1==1, Y1==0))/N,
               nrow(population%>%filter(X1==0, Y1==1))/N,
               nrow(population%>%filter(X1==1, Y1==1))/N)

for (k in 1:4){
  p_true = joint_true[k]
  var = V_joint[k]
  CI_joint[k] = 1 - sum(joint_prob[,k] - 2.009*sqrt(var) > p_true | joint_prob[,k] + 2.009*sqrt(var) < p_true)/n_sim
}

CI.list = list(CI.T = c(CI.T1, CI.T2, CI.T3, CI.T4, CI.T5, CI.T6), 
               CI.cond = c(CI.cond1, CI.cond2, CI.cond3, CI.cond4, CI.cond5, CI.cond6, CI.cond7, CI.cond8),
               CI.joint = CI_joint, CI.given12 = CI_given12)

VarMDAM.list = list(Premissvar.T = apply(T_premiss, 2, var), 
                    Premissvar.cond = apply(cond_prob_premiss, 2, var),
                    Premissvar.joint = apply(joint_prob_premiss, 2, var),
                    Premissvar.given12 = apply(cond_prob_given12_premiss, 2, var),
                    var.T = apply(T_est, 2, var), var.cond = apply(cond_prob, 2, var),
                    var.joint = apply(joint_prob, 2, var), var.given12 = apply(cond_prob_given12, 2, var),
                    AvgEstVar.T = colMeans(V_T), AvgEstVar.cond = colMeans(V_cond),
                    AvgEstVar.joint = colMeans(V_joint), AvgEstVar.given12 = colMeans(V_given12))

