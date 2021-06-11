#FUNCTION-NAME-INPUT CODE

#function input for test run:


#Sim state function for GARCH model of the kind: y(t+1) = a0+a1*y(t) + sigma(t)*epsilon(t)
#                                                sigma(t)^2 = b0 +b1*sigma(t-1)^2 + b3*epsilon(t-1)^2
# State is y(t), sigma(t+1), epsilon(t)

sim_state_GARCH=function(state,GARCH_param,n,time=0){

  # this function requires either state to be a vector, or n to be equal to the number of rows of the
  # state matrix!
  if(is.matrix(state)==FALSE){
    state=matrix(state,nrow=1)
  }
  x=rnorm(n)
  y_tplus1 = GARCH_param[1] + GARCH_param[2]*state[,1]+state[,2]*rnorm(n)
  sigma_tplus2 = sqrt(GARCH_param[3] +GARCH_param[4]*state[,2]^2 +GARCH_param[5]*y_tplus1^2)
  
  return(cbind(y_tplus1,sigma_tplus2))
}



sim_state_GARCH_unconditional<-function(t,state_0,sim_state_GARCH,GARCH_param,n){
    state=state_0
    for(s in 1:t){
     state=sim_state_GARCH(state,GARCH_param,n)
    }
    return(state)
  }


#As the name indicates, a function which produces as covariates: 1, the state and all square- and cross-terms
square_basis_vector_GARCH <- function(state){
  if(is.matrix(state)==FALSE){
    state_inc_1 =c(1,state)
    output =c()
    for(i in seq(1:length(state_inc_1))){
      output = cbind(output,matrix(state_inc_1[i]*state_inc_1[1:i],nrow=1))
    }
    return(as.vector(output))
  }else{
    state_inc_1 = cbind(1,state)
    output =c()
    for(i in seq(1:length(state_inc_1[1,]))){
      output = cbind(output,state_inc_1[,i]*state_inc_1[,1:i])
    }
    return(output)
  }}

dummy_basis_vector <- function(state){
  l =length(state[,1])
  return(matrix(0,l,6))
}

#note that in the function below, the state is always a vector (i.e. the input is a single state)


#RUN MODEL

library(doParallel)

cl <- makeCluster(8)

registerDoParallel(cl)

T = 6
eta = 0.06
p=0.995
N_reg=10^4
N_VaR=10^5
state_0=c(0,1)
GARCH_param=c(1,1,0.1,0.1,0.1)
 


valuation_list_one_GARCH_201203 = LSM_COC_Valuation(N_reg, N_VaR,eta,p,T,state_0,GARCH_param, 
                                               sim_state_GARCH, sim_state_GARCH_unconditional,
                                               square_basis_vector_GARCH)

save.image()


thetamatrix_GARCH_201203 = valuation_list_GARCH_201203[[3]]

phimatrix_GARCH_201203 = valuation_list_GARCH_201203[[4]]

N_val = 10^4
N_VaR_val = 10^5

validation_list_GARCH_201203<-LSM_COC_Validation(N_val, N_VaR_val,eta,p,T,state_0,GARCH_param,thetamatrix_GARCH_201203, 
                                                phimatrix_GARCH_201203, sim_state_GARCH, 
                                                sim_state_GARCH_unconditional,square_basis_vector_GARCH)





val_return_matrix = validation_list_GARCH_201203[[1]]

val_prob_matrix = validation_list_GARCH_201203[[2]]

val_VaR_matrix = validation_list_GARCH_201203[[3]]

val_Eval_matrix = validation_list_GARCH_201203[[4]]

val_VaR_est_matrix = validation_list_GARCH_201203[[5]]

val_Eval_est_matrix = validation_list_GARCH_201203[[6]]




s=3

hist(100*(val_return_matrix[s,]-1),10,xlab="Return (%)",main = paste("Excess return at time" , s),xlim=c(0,15))

mean(100*(val_return_matrix[s,]-1))

hist(1-val_prob_matrix[s,],10,xlab="Probability",main = paste("Probability of ruin at time" , s),xlim=c(0.003,0.007))

mean(1-val_prob_matrix[s,],1000,xlab="Probability",main = paste("Probability of ruin at time" , s))

sqrt(mean((val_Eval_matrix-val_Eval_est_matrix)^2))

hist(val_Eval_matrix)

max(val_Eval_est_matrix)

hist(val_VaR_est_matrix)















































number_of_basis_functions=length(square_basis_vector_GARCH(state_0))


thetamatrix = matrix(rep(0,number_of_basis_functions*T),ncol=T)
phimatrix = matrix(rep(0,number_of_basis_functions*T),ncol=T)
 
ptime <- system.time({
 
for(i in seq(T-1,1)){
       #i=2
  
  print(i)
  theta = thetamatrix[,i+1]
  phi = phimatrix[,i+1]
       
   #parameters describing the value - linear combination of theta and phi
         
  beta = theta-(1/(1+eta))*phi
  print(beta)
  thetaphi = VaR_Value_regressionfn_GARCH(state_0,GARCH_param,beta,p,N_reg,N_VaR,time=i,
                                          sim_state_GARCH,sim_state_GARCH_uncond,VaR_simfn_GARCH,square_basis_vector_GARCH)
  thetamatrix[,i]=thetaphi[1:number_of_basis_functions]
  phimatrix[,i]=thetaphi[(number_of_basis_functions+1):(2*number_of_basis_functions)]
  thetamatrix[,i][is.na(thetamatrix[,i])] <- 0
  phimatrix[,i][is.na(phimatrix[,i])] <- 0
  }
})
ptime[[3]]
#Calculating time zero R and V
 
theta = thetamatrix[,1]
phi = phimatrix[,1]
  
beta = theta - (1/(1+eta))*phi
  
Capreq_Eval = VaR_Value_regressionfn_GARCH(state_0,GARCH_param,beta,p,N_reg,N_VaR*10,time=0,sim_state_GARCH,sim_state_GARCH_uncond,VaR_simfn_GARCH)
  
Value = Capreq_Eval[1]-(1/(1+eta))*Capreq_Eval[2]
Value
 
###########################################################################################################################
# PROPER VALIDATION OF MODEL ##############################################################################################
###########################################################################################################################
#saved from first run, june 25:
#thetamatrix_saved1
#phimatrix_saved1
#
# saved from june 27
# val_prob_matrix_saved1=val_prob_matrix
# val_prob_return_saved1=val_return_matrix
#
#saved matrices july 2 (actual GARCH model)
#val_prob_matrix_saved2 = val_prob_matrix
#val_return_matrix_saved2 = val_return_matrix
#val_VaR_est_matrix_saved2 =val_VaR_est_matrix
#val_VaR_matrix_saved2 =val_VaR_matrix
#val_Eval_est_matrix_saved2 = val_Eval_est_matrix
#val_Eval_est_matrix_saved2 =val_Eval_est_matrix
#
#thetamatrix_saved2 =thetamatrix
#phimatrix_saved2 =phimatrix
# 
# N_val=10^3
# N_VaR_val=10^4
# 
# val_prob_matrix =matrix(0,T-1,N_val)
# val_return_matrix =matrix(0,T-1,N_val)
# val_VaR_matrix =matrix(0,T-1,N_val)
# val_Eval_matrix =matrix(0,T-1,N_val)
# val_VaR_est_matrix =matrix(0,T-1,N_val)
# val_Eval_est_matrix =matrix(0,T-1,N_val)
# 
# 
# for(time in seq(1,(T-1))){
# 
# theta=thetamatrix[,time]
# phi=phimatrix[,time]
# theta_next=thetamatrix[,time+1]
# phi_next=phimatrix[,time+1]
# 
# 
# beta_next =theta_next - (1/(1+eta))*phi_next
# 
# states=sim_state_GARCH_uncond(time,c(0,1),sim_state_GARCH,GARCH_param,N_val)
# 
# VaR_vector_est = square_basis_vector_GARCH(states)%*%theta
# E_vector_est = square_basis_vector_GARCH(states)%*%phi
# 
# val_VaR_est_matrix[time,]=VaR_vector_est
# val_Eval_est_matrix[time,]=E_vector_est
# 
# #N_val=10
# val_returns = rep(0,N_val)
# val_probs = rep(0,N_val)
# val_VaR = rep(0,N_val)
# val_Eval = rep(0,N_val)
# 
# for(i in seq(1,N_val)){
#   if(i%%100==0){print(i)}
#   
#   
#   prob_return = VaR_Value_validationfn_GARCH(eta,p,VaR_vector_est[i],E_vector_est[i],N_VaR_val,states[i,],beta_next,GARCH_param
#                                              ,sim_state_GARCH,square_basis_vector_GARCH,time)
#   val_probs[i]=prob_return[1]
#   val_returns[i]=prob_return[2]
#   val_VaR[i] = prob_return[3]
#   val_Eval[i] = prob_return[4]
#   
#   
# }
# 
# val_return_matrix[time,]=val_returns
# val_prob_matrix[time,]=val_probs
# val_VaR_matrix[time,]=val_VaR
# val_Eval_matrix[time,]=val_Eval
# 
# }









list <- LSM_COC_Validation(100, 1000,eta,p,T,state_0,GARCH_param,thetamatrix, phimatrix, sim_state_GARCH, sim_state_GARCH_uncond,square_basis_vector_GARCH)

val_return_matrix = list[[1]]

val_prob_matrix = list[[2]]

val_VaR_matrix = list[[3]]

val_Eval_matrix = list[[4]]

val_VaR_est_matrix = list[[5]]

val_Eval_est_matrix = list[[6]]

s=3

hist(100*(val_return_matrix[s,]-1),100,xlab="Return (%)",main = paste("Excess return at time" , s))

hist(1-val_prob_matrix[s,],100,xlab="Probability",main = paste("Probability of ruin at time" , s))


#mean squared errors

#differences

VaR_diff=val_VaR_est_matrix-val_VaR_matrix

Evalue_diff=val_Eval_est_matrix-val_Eval_matrix

Value_diff = VaR_diff-(1/(1+eta))*Evalue_diff

val_Values=val_VaR_matrix-(1/(1+eta))*val_Eval_matrix

#Mean absolute error:

s=5
mean(abs(val_Values[s,]))

sqrt(mean((Value_diff[s,])^2))

hist(val_Values[s,])


hist(Value_diff[s,],50)


###########################################################################################################################
# NAIVE VALIDATION OF MODEL ###############################################################################################
###########################################################################################################################



