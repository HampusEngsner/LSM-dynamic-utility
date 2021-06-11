
#RUN MODELS !!!!!

library(doParallel)

cl <- makeCluster(8)

registerDoParallel(cl)

# LIFE, 4 cohorts

timer1 = Sys.time()

T = 6
eta = 0.06
p=0.995
N_reg=0.5*10^5
N_VaR=10^5

param_stock = c(0.03, 0.1,100,0.03,0.1,100,0.4, 1)
#guarenteed death and survival benefits, plus terminal time
param_benifits = c(100,110,T)
#parameters for cohorts of form N_0, and p- vectors. In this case, T=2.
param_N0 = c(1000,1000,1000,1000)

param_pvectors = cbind(life_M90(50:(50+T)),life_M90(60:(60+T)),life_M90(70:(70+T)),life_M90(80:(80+T)))

LIFE_param = list(param_stock,param_benifits,param_N0,param_pvectors)

state_0 = c(0,LIFE_param[[1]][3],LIFE_param[[3]],LIFE_param[[1]][6],0)



valuation_list_life_4cohorts_201204 = LSM_COC_Valuation(N_reg, N_VaR,eta,p,T,state_0,LIFE_param, 
                                               sim_state_LIFE, sim_state_LIFE_unconditional,
                                               square_basis_vector_LIFE)



timer2 = Sys.time()

save.image()

thetamatrix_LIFE_4cohorts_201204 = valuation_list_life_4cohorts_201204[[3]]
  
phimatrix_LIFE_4cohorts_201204 = valuation_list_life_4cohorts_201204[[4]]

N_val = 10^4
N_VaR_val = 10^5



validation_list_life_4cohorts_201204<-LSM_COC_Validation(N_val, N_VaR_val,eta,p,T,state_0,LIFE_param, 
                                                         thetamatrix_LIFE_4cohorts_201204,
                                                         phimatrix_LIFE_4cohorts_201204, sim_state_LIFE, 
                                                         sim_state_LIFE_unconditional,square_basis_vector_LIFE)



timer3 = Sys.time()
save.image()


# LIFE, 10 cohorts

T = 6
eta = 0.06
p=0.995
N_reg=0.5*10^5
N_VaR=10^5
param_stock = c(0.03, 0.1,100,0.03,0.1,100,0.4, 1)
#guarenteed death and survival benefits, plus terminal time
param_benifits = c(100,110,T)
#parameters for cohorts of form N_0, and p- vectors. In this case, T=2.
param_N0 = rep(1000,10)
param_pvectors = cbind(life_M90(40:(40+T)), life_M90(45:(45+T)),life_M90(50:(50+T)),life_M90(55:(55+T)),
                       life_M90(60:(60+T)), life_M90(65:(65+T)),life_M90(70:(70+T)),  life_M90(75:(75+T)), 
                       life_M90(80:(80+T)), life_M90(85:(85+T)))
LIFE_param = list(param_stock,param_benifits,param_N0,param_pvectors)

state_0 = c(0,LIFE_param[[1]][3],LIFE_param[[3]],LIFE_param[[1]][6],0)

valuation_list_life_10cohorts_201204 = LSM_COC_Valuation(N_reg, N_VaR,eta,p,T,state_0,LIFE_param, 
                                                        sim_state_LIFE, sim_state_LIFE_unconditional,
                                                        square_basis_vector_LIFE)

timer4 = Sys.time()
save.image()

thetamatrix_LIFE_10cohorts_201204 = valuation_list_life_10cohorts_201204[[3]]

phimatrix_LIFE_10cohorts_201204 = valuation_list_life_10cohorts_201204[[4]]

N_val = 10^4
N_VaR_val = 10^5



validation_list_life_10cohorts_201204<-LSM_COC_Validation(N_val, N_VaR_val,eta,p,T,state_0,LIFE_param,thetamatrix_LIFE_10cohorts_201204, 
                                                          phimatrix_LIFE_10cohorts_201204, sim_state_LIFE, 
                                                         sim_state_LIFE_unconditional,square_basis_vector_LIFE)


timer5 = Sys.time()
save.image()


# One GARCH model

T = 6
eta = 0.06
p=0.995
N_reg=10^4
N_VaR=10^5
state_0=c(0,1)
GARCH_param=c(1,1,0.1,0.1,0.1)


valuation_list_one_GARCH_201204 = LSM_COC_Valuation(N_reg, N_VaR,eta,p,T,state_0,GARCH_param, 
                                                    sim_state_GARCH, sim_state_GARCH_unconditional,
                                                    square_basis_vector_GARCH)

N_reg=10^4
N_VaR=10^5

timer6 = Sys.time()

save.image()


thetamatrix_one_GARCH_201204 = valuation_list_one_GARCH_201204[[3]]

phimatrix_one_GARCH_201204 = valuation_list_one_GARCH_201204[[4]]

N_val = 10^4
N_VaR_val = 10^5

validation_list_one_GARCH_201204<-LSM_COC_Validation(N_val, N_VaR_val,eta,p,T,state_0,GARCH_param,thetamatrix_one_GARCH_201204, 
                                                      phimatrix_one_GARCH_201204, sim_state_GARCH, 
                                                 sim_state_GARCH_unconditional,square_basis_vector_GARCH)


timer7 = Sys.time()

save.image()



# Ten GARCH models

T = 6
eta = 0.06
p=0.995
N_reg=10^4
N_VaR=10^5
state_0=rep(c(0,1),10)
GARCH_param=c(1,1,0.1,0.1,0.1)


valuation_list_ten_GARCH_201204 =M=10^4$, $n = 10^5$  LSM_COC_Valuation(N_reg, N_VaR,eta,p,T,state_0,GARCH_param, 
                                                    sim_state_ten_GARCH, sim_state_ten_GARCH_unconditional,
                                                    square_basis_vector_ten_GARCH)


timer8 = Sys.time()

save.image()


thetamatrix_ten_GARCH_201204 = valuation_list_ten_GARCH_201204[[3]]

phimatrix_ten_GARCH_201204 = valuation_list_ten_GARCH_201204[[4]]

N_val = 10^4
N_VaR_val = 10^5

validation_list_ten_GARCH_201204<-LSM_COC_Validation(N_val, N_VaR_val,eta,p,T,state_0,GARCH_param,thetamatrix_ten_GARCH_201204, 
                                                     phimatrix_ten_GARCH_201204, sim_state_ten_GARCH, 
                                                     sim_state_ten_GARCH_unconditional,square_basis_vector_ten_GARCH)


timer9 = Sys.time()

save.image()






# Test run with dynamic basis function choices

#NOT READY YET!!!!!

T = 6
eta = 0.06
p=0.9
N_reg=10^4
N_VaR=10^4
p=0.9
param_stock = c(0.03, 0.1,100,0.03,0.1,100,0.4, 1)
#guarenteed death and survival benefits, plus terminal time
param_benifits = c(100,110,T)
#parameters for cohorts of form N_0, and p- vectors. In this case, T=2.
param_N0 = c(1000,1000,1000,1000)

param_pvectors = cbind(life_M90(50:(50+T)),life_M90(60:(60+T)),life_M90(70:(70+T)),life_M90(80:(80+T)))

LIFE_param = list(param_stock,param_benifits,param_N0,param_pvectors)

state_0 = c(0,LIFE_param[[1]][3],LIFE_param[[3]],LIFE_param[[1]][6],0)



valuation_list_life_4cohorts_fitReLucovariates_201204 = LSM_COC_Valuation_fitReLucovariates(N_reg, N_VaR,eta,p,T,state_0,LIFE_param, 
                                                        sim_state_LIFE, sim_state_LIFE_unconditional,
                                                        square_basis_vector_LIFE_fitReLucovariates)


save.image()

thetamatrix_LIFE_4cohorts_fitReLucovariates_201204 = valuation_list_life_4cohorts_fitReLucovariates_201204[[3]]

phimatrix_LIFE_4cohorts_fitReLucovariates_201204 = valuation_list_life_4cohorts_fitReLucovariates_201204[[4]]

N_val = 10^3
N_VaR_val = 10^3



validation_list_life_4cohorts_fitReLucovariates_201204<-LSM_COC_Validation(N_val, N_VaR_val,eta,p,T,state_0,LIFE_param,
                                                                           thetamatrix_LIFE_4cohorts_fitReLucovariates_201204, phimatrix_LIFE_4cohorts_fitReLucovariates_201204, sim_state_LIFE, 
                                                         sim_state_LIFE_unconditional,square_basis_vector_LIFE_fitReLucovariates)






# validation_list_one_GARCH_201204

# validation_list_ten_GARCH_201204

# validation_list_life_4cohorts_201204

#   validation_list_life_10cohorts_201204

temp_val_list = validation_list_life_10cohorts_201204

val_return_matrix = temp_val_list[[1]]
val_prob_matrix = temp_val_list[[2]]
val_VaR_matrix = temp_val_list[[3]]
val_Eval_matrix = temp_val_list[[4]]
val_VaR_est_matrix = temp_val_list[[5]]
val_Eval_est_matrix = temp_val_list[[6]]
val_V_est_matrix = val_VaR_est_matrix - (1/(1+eta))*val_Eval_est_matrix
val_V_matrix = val_VaR_matrix - (1/(1+eta))*val_Eval_matrix

V_RMSE = c()
V_RRMSE= c()
R_RMSE= c()
R_RRMSE= c()
E_RMSE= c()
E_RRMSE= c()

upper_eta = c()
lower_eta = c()
upper_alpha = c()
lower_alpha = c()



for(s in 1:5){
  V_RMSE[s]=sqrt(mean((val_V_matrix[s,]-val_V_est_matrix[s,])^2))
  V_RRMSE[s]=sqrt(mean((val_V_matrix[s,]-val_V_est_matrix[s,])^2))/sqrt(mean((val_V_matrix[s,])^2))
  
  E_RMSE[s]=sqrt(mean((val_Eval_matrix[s,]-val_Eval_est_matrix[s,])^2))
  E_RRMSE[s]=sqrt(mean((val_Eval_matrix[s,]-val_Eval_est_matrix[s,])^2))/sqrt(mean((val_Eval_matrix[s,])^2))
  
  R_RMSE[s]=sqrt(mean((val_VaR_matrix[s,]-val_VaR_est_matrix[s,])^2))
  R_RRMSE[s]=sqrt(mean((val_VaR_matrix[s,]-val_VaR_est_matrix[s,])^2))/sqrt(mean((val_VaR_matrix[s,])^2))
  
  upper_eta[s] = quantile(100*(val_return_matrix[s,]-1), 0.975)
  lower_eta[s] = quantile(100*(val_return_matrix[s,]-1), 0.025)
  upper_alpha[s] = quantile(1-val_prob_matrix[s,], 0.975)
  lower_alpha[s] = quantile(1-val_prob_matrix[s,], 0.025)
  }

round(V_RMSE,4)
round(100*V_RRMSE,4)

round(R_RMSE,4)
round(100*R_RRMSE,4)

round(E_RMSE,4)
round(100*E_RRMSE,4)

round(100*upper_alpha,3)
round(100*lower_alpha,3)
round(upper_eta,2)
round(lower_eta,2)

s=3
hist(100*(val_return_matrix[s,]-1),100,xlab="Return (%)",main = paste("Excess return at time" , s),xlim=c(0,15))

hist(100*(1.06*(val_Eval_matrix[s,])/val_Eval_est_matrix[s,]-1),50,xlab="Return (%)",main = paste("Excess return at time" , s),xlim=c(-10,25))

mean(100*(1.06*(val_Eval_matrix[s,] )/val_Eval_est_matrix[s,]-1))

mean(100*(val_return_matrix[s,]-1))

hist(1-val_prob_matrix[s,],100,xlab="Probability",main = paste("Probability of ruin at time" , s),xlim=c(0.002,0.009))

mean(1-val_prob_matrix[s,])

sqrt(mean((val_Eval_matrix[s,]-val_Eval_est_matrix[s,])^2))/sqrt(mean((val_Eval_matrix[s,])^2))

hist(val_Eval_matrix)

max(val_Eval_est_matrix)

hist(val_VaR_est_matrix)





# re- run, temporary code



retime1 = Sys.time()

T = 6
eta = 0.06
p=0.995

state_0=rep(c(0,1),10)
GARCH_param=c(1,1,0.1,0.1,0.1)
N_val = 10^4
N_VaR_val = 10^5

validation_list_ten_GARCH_201204<-LSM_COC_Validation(N_val, N_VaR_val,eta,p,T,state_0,GARCH_param,thetamatrix_ten_GARCH_201204, 
                                                     phimatrix_ten_GARCH_201204, sim_state_ten_GARCH, 
                                                     sim_state_ten_GARCH_unconditional,square_basis_vector_ten_GARCH)

retime2 = Sys.time()

param_stock = c(0.03, 0.1,100,0.03,0.1,100,0.4, 1)
#guarenteed death and survival benefits, plus terminal time
param_benifits = c(100,110,T)
#parameters for cohorts of form N_0, and p- vectors. In this case, T=2.
param_N0 = c(1000,1000,1000,1000)

param_pvectors = cbind(life_M90(50:(50+T)),life_M90(60:(60+T)),life_M90(70:(70+T)),life_M90(80:(80+T)))

LIFE_param = list(param_stock,param_benifits,param_N0,param_pvectors)

state_0 = c(0,LIFE_param[[1]][3],LIFE_param[[3]],LIFE_param[[1]][6],0)
N_val = 10^4
N_VaR_val = 10^5



validation_list_life_4cohorts_201204<-LSM_COC_Validation(N_val, N_VaR_val,eta,p,T,state_0,LIFE_param, 
                                                         thetamatrix_LIFE_4cohorts_201204,
                                                         phimatrix_LIFE_4cohorts_201204, sim_state_LIFE, 
                                                         sim_state_LIFE_unconditional,square_basis_vector_LIFE)







retime3 = Sys.time()

save.image()




