#FUNCTION-NAME-INPUT CODE

#function input for test run:


#Sim state function for GARCH model of the kind: y(t+1) = a0+a1*y(t) + sigma(t)*epsilon(t)
#                                                sigma(t)^2 = b0 +b1*sigma(t-1)^2 + b3*epsilon(t-1)^2
# State is y(t), sigma(t+1)

sim_state_one_GARCH=function(state,GARCH_param,n,time=0){

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

sim_state_ten_GARCH=function(state,GARCH_param,n,time=0){
  
  if(is.matrix(state)==FALSE){
    state=matrix(state,nrow=1)
  }
  
  
  output_state =sim_state_one_GARCH(state[,1:2],GARCH_param,n,time)
  
  for(i in 2:10){
    output_state = cbind(output_state,sim_state_one_GARCH(state[,(2*i-1):(2*i)],GARCH_param,n,time))
  }
  
  return(output_state)
}

sim_state_ten_GARCH_unconditional<-function(t,state_0,sim_state_ten_GARCH,GARCH_param,n){
    state=state_0
    for(s in 1:t){
     state=sim_state_ten_GARCH(state,GARCH_param,n)
    }
    return(state)
}

square_basis_vector_ten_GARCH <- function(state){
  if(is.matrix(state)==FALSE){
    indices_odd = seq(1,9,2)
    indices_even = seq(2,10,2)  
    sum_state=cbind(sum(state[indices_odd]), sqrt(sum(state[indices_even]^2)))
    state_inc_1 =c(1,sum_state)
    output =c()
    for(i in seq(1:length(state_inc_1))){
      output = cbind(output,matrix(state_inc_1[i]*state_inc_1[1:i],nrow=1))
    }
    return(c(as.vector(output)[-2],state))
  }else{
    #concatination:
    indices_odd = seq(1,9,2)
    indices_even = seq(2,10,2)  
    sum_state=cbind(rowSums(state[,indices_odd]), sqrt(rowSums(state[,indices_even]^2)))
    state_inc_1 = cbind(1,sum_state)
    output =c()
    for(i in seq(1:length(state_inc_1[1,]))){
      output = cbind(output,state_inc_1[,i]*state_inc_1[,1:i])
    }
    return(cbind(output[,-2],state))
  }}



#RUN MODEL

library(doParallel)

cl <- makeCluster(7)

registerDoParallel(cl)

T = 6
eta = 0.06
p=0.995 
N_reg=10^4
N_VaR=10^5
state_0=c(0,1)
GARCH_param=c(1,1,0.1,0.1,0.1)



valuation_list_GARCH_201203 = LSM_COC_Valuation(N_reg, N_VaR,eta,p,T,state_0,GARCH_param, 
                                                sim_state_GARCH, sim_state_GARCH_unconditional,
                                                square_basis_vector_GARCH)











