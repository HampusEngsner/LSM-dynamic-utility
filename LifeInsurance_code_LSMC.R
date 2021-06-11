#FUNCTION-NAME-INPUT CODE

#function input for test run:

#Sim state function for 

#state description

life_M90<-function(x){
  alpha=0.001
  beta= 0.000012
  gamma=0.101314
  
  l=exp(-alpha*x-(beta/gamma)*(exp(gamma*x)-1))
  
  return(l)
  
}

call_option_function<-function(S,K,mu,sigma,T,t){
  d1 = (1/sqrt(sigma^2*(T-t)))*(log(S/K)+(mu+sigma^2/2)*(T-t)) +sqrt(sigma^2*(T-t))/2
    
  d2=d1 - sqrt(sigma^2*(T-t))
  
  return(K+exp((mu+sigma^2/2)*(T-t))*pnorm(d1)*S - pnorm(d2)*K )
}


sim_state_LIFE=function(state,LIFE_param,n,time=0){

  # this function requires either state to be a vector, or n to be equal to the number of rows of the
  # state matrix!
  
  if(is.matrix(state)==FALSE){
    state=matrix(state,nrow=1)
  }
  
  cohorts = length(LIFE_param[[3]])
  
  brownian1 = rnorm(n)
  
  brownian2 = LIFE_param[[1]][7]*brownian1 + sqrt(1-LIFE_param[[1]][7]^2)*rnorm(n)
  
  ftplus1 = state[,2]*exp(LIFE_param[[1]][1]+ LIFE_param[[1]][2]*brownian1)
  
  ytplus1  = state[,3+cohorts]*exp(LIFE_param[[1]][4]+ LIFE_param[[1]][5]*brownian2)
  
  alive_tplus1 =  matrix(rbinom(n*cohorts, t(state[,3:(2+cohorts)]), LIFE_param[[4]][time+1,]),  ncol=cohorts, byrow = TRUE)
  
  
  if(nrow(state)==1){
    deaths = -t(apply(alive_tplus1, 1, '-',state[,3:(2+cohorts)]))
  }else{
    deaths =state[,3:(2+cohorts)] - alive_tplus1
  }
  
  deaths = rowSums(deaths)
  
  payments = deaths*(sapply(ftplus1,max,LIFE_param[[2]][1]) - LIFE_param[[1]][8]*ytplus1) #+ 100000*exp(rnorm(n,0,0.1))
  
  
  if(time == LIFE_param[[2]][3]-1){
    payments = payments + rowSums(alive_tplus1)*(sapply(ftplus1,max,LIFE_param[[2]][2]) - LIFE_param[[1]][8]*ytplus1)
  }
  return(cbind(payments,ftplus1,alive_tplus1,ytplus1,time+1))
}


#hist(sim_state_LIFE(state_0, LIFE_param, 100,5)[,1])

sim_state_LIFE_unconditional<-function(t,state_0,sim_state,LIFE_param,n){
  
  #simple code, simulate cohorts and stock price! No "last payment" will be simulated.
  cohorts = length(LIFE_param[[3]])
  
  brownian1 = rnorm(n)
  
  brownian2 = LIFE_param[[1]][7]*brownian1 + sqrt(1-LIFE_param[[1]][7]^2)*rnorm(n)
  
  ft= state_0[2]*exp(t*LIFE_param[[1]][1]+ sqrt(t)*LIFE_param[[1]][2]*brownian1)
  
  yt  = state_0[3+cohorts]*exp(t*LIFE_param[[1]][4]+  sqrt(t)*LIFE_param[[1]][5]*brownian2)
  
  survival_prob = apply(matrix((LIFE_param[[4]][1:t,]), ncol = cohorts),2,prod)
  
  alive = matrix(rbinom(n*cohorts, t(state_0[3:(2+cohorts)]), survival_prob),  ncol=cohorts, byrow = TRUE)
  return(cbind(rep(0,n), ft, alive,yt,t))
  }


#make sure first component of state is NOT included as it does not pertain to future payments!!!


square_basis_vector_LIFE<- function(state){
  
  cohorts = length(LIFE_param[[3]])
  
  # Empirically best K - values for time t= 5
  K_1 = 200
  K_2 = 162
  K_3 = 124
  K_4 = 103
  if(is.matrix(state)==FALSE){
   time=0 
  }else{
      time = tail(state[1,], n=1)
  }
  
  #modifying constants w.r.t. time
  time_factor = exp(LIFE_param[[1]][1]*(time-5) + 0.5*LIFE_param[[1]][2]^2*(time-5))
  K_1 = K_1*time_factor
  K_2 = K_2*time_factor
  K_3 = K_3*time_factor
  K_4 = K_4*time_factor
  
  
  
  if(is.matrix(state)==FALSE){
    
    state_inc_1 =c(1,state[-1])
    t =  state_inc_1[length(state_inc_1)]
    alive = state_inc_1[3:(length(state)-2)]
    
    alive_sums = sum(alive)
    
    means_deaths = sum(alive*(1-LIFE_param[[4]][t+1,])) 
    
    stdevs =  sqrt(sum(alive*(LIFE_param[[4]][t+1,]*(1-LIFE_param[[4]][t+1,]))))
    
    
    max1 = sapply(state_inc_1[2],max,K_1)-K_1
    max2 = sapply(state_inc_1[2],max,K_2)-K_2
    max3 = sapply(state_inc_1[2],max,K_3)-K_3
    max4 = sapply(state_inc_1[2],max,K_4)-K_4

    option1 = call_option_function(state_inc_1[2],LIFE_param[[2]][1],LIFE_param[[1]][1],LIFE_param[[1]][2],t+1,t)
    
    
    option2 = call_option_function(state_inc_1[2],LIFE_param[[2]][2],LIFE_param[[1]][1],LIFE_param[[1]][2],T,t)
    
    
    output =c(1,state_inc_1[2:(length(state)-1)])
    output = c(output,
               stdevs*max1,
               stdevs*max2,
               stdevs*max3,
               stdevs*max4,
               stdevs*max1*state_inc_1[3+cohorts], 
               stdevs*max2*state_inc_1[3+cohorts], 
               stdevs*max3*state_inc_1[3+cohorts], 
               stdevs*max4*state_inc_1[3+cohorts], 
               stdevs*state_inc_1[3+cohorts],
               stdevs*state_inc_1[2],
               stdevs*state_inc_1[2]^2,
               stdevs*state_inc_1[2]^3, 
               stdevs*state_inc_1[2]*state_inc_1[3+cohorts],
               stdevs*state_inc_1[2]^2*state_inc_1[3+cohorts],
               stdevs*state_inc_1[3+cohorts]^2,
               stdevs*option1, 
               stdevs*option2,
               stdevs*option1*state_inc_1[3+cohorts], 
               stdevs*option2*state_inc_1[3+cohorts], 
               means_deaths*max1, 
               means_deaths*max2,
               means_deaths*max3,
               means_deaths*max4,
               means_deaths*max1*state_inc_1[3+cohorts],
               means_deaths*max2*state_inc_1[3+cohorts],
               means_deaths*max3*state_inc_1[3+cohorts],
               means_deaths*max4*state_inc_1[3+cohorts],
               means_deaths*state_inc_1[3+cohorts],
               means_deaths*state_inc_1[2], 
               means_deaths*state_inc_1[2]^2, 
               means_deaths*state_inc_1[2]^3, 
               means_deaths*state_inc_1[2]*state_inc_1[3+cohorts],
               means_deaths*state_inc_1[2]^2*state_inc_1[3+cohorts],
               means_deaths*state_inc_1[3+cohorts]^2,
               means_deaths*option1,
               means_deaths*option2, 
               means_deaths*option1*state_inc_1[3+cohorts], 
               means_deaths*option2*state_inc_1[3+cohorts], 
               alive_sums*max1,
               alive_sums*max2,
               alive_sums*max3,
               alive_sums*max4,
               alive_sums*max1*state_inc_1[3+cohorts],
               alive_sums*max2*state_inc_1[3+cohorts], 
               alive_sums*max3*state_inc_1[3+cohorts],
               alive_sums*max4*state_inc_1[3+cohorts], 
               alive_sums*state_inc_1[3+cohorts],
               alive_sums*state_inc_1[2],
               alive_sums*state_inc_1[2]^2,
               alive_sums*state_inc_1[2]^3, 
               alive_sums*state_inc_1[2]*state_inc_1[3+cohorts],
               alive_sums*state_inc_1[2]^2*state_inc_1[3+cohorts],
               alive_sums*state_inc_1[3+cohorts]^2,
               alive_sums*option1,
               alive_sums*option2, 
               alive_sums*option1*state_inc_1[3+cohorts],
               alive_sums*option2*state_inc_1[3+cohorts] 
    )
    
    return(as.vector(output))
  }else{
    
    state_inc_1 =cbind(1,state[,-1])
    t =  state_inc_1[1,][length(state_inc_1[1,])]
    
    alive = state_inc_1[,3:(length(state[1,])-2)]
    
    alive_sums = rowSums(alive)
    
    means_deaths = alive%*%(1-LIFE_param[[4]][t+1,]) 
    
    stdevs =  sqrt(alive%*%(LIFE_param[[4]][t+1,]*(1-LIFE_param[[4]][t+1,])))
    
    
    max1 = sapply(state_inc_1[,2],max,K_1)-K_1
    max2 = sapply(state_inc_1[,2],max,K_2)-K_2
    max3 = sapply(state_inc_1[,2],max,K_3)-K_3
    max4 = sapply(state_inc_1[,2],max,K_4)-K_4
    
    
    option1 = call_option_function(state_inc_1[,2],LIFE_param[[2]][1],LIFE_param[[1]][1],LIFE_param[[1]][2],t+1,t)
    option2 = call_option_function(state_inc_1[,2],LIFE_param[[2]][2],LIFE_param[[1]][1],LIFE_param[[1]][2],T,t)
    
    
    
    
    output = cbind(1,state_inc_1[,2:(length(state[1,])-1)])
    
    output =cbind(output,
                  stdevs*max1,
                  stdevs*max2,
                  stdevs*max3,
                  stdevs*max4,
                  stdevs*max1*state_inc_1[,3+cohorts],
                  stdevs*max2*state_inc_1[,3+cohorts],
                  stdevs*max3*state_inc_1[,3+cohorts],
                  stdevs*max4*state_inc_1[,3+cohorts],
                  stdevs*state_inc_1[,3+cohorts],
                  stdevs*state_inc_1[,2],
                  stdevs*state_inc_1[,2]^2,
                  stdevs*state_inc_1[,2]^3, 
                  stdevs*state_inc_1[,2]*state_inc_1[,3+cohorts],
                  stdevs*state_inc_1[,2]^2*state_inc_1[,3+cohorts],
                  stdevs*state_inc_1[,3+cohorts]^2,                 
                  stdevs*option1,
                  stdevs*option2,
                  stdevs*option1*state_inc_1[,3+cohorts],
                  stdevs*option2*state_inc_1[,3+cohorts],
                  means_deaths*max1,
                  means_deaths*max2,
                  means_deaths*max3,
                  means_deaths*max4,
                  means_deaths*max1*state_inc_1[,3+cohorts],
                  means_deaths*max2*state_inc_1[,3+cohorts],
                  means_deaths*max3*state_inc_1[,3+cohorts],
                  means_deaths*max4*state_inc_1[,3+cohorts],
                  means_deaths*state_inc_1[,3+cohorts],
                  means_deaths*state_inc_1[,2],
                  means_deaths*state_inc_1[,2]^2,
                  means_deaths*state_inc_1[,2]^3, 
                  means_deaths*state_inc_1[,2]*state_inc_1[,3+cohorts],
                  means_deaths*state_inc_1[,2]^2*state_inc_1[,3+cohorts],
                  means_deaths*state_inc_1[,3+cohorts]^2,                 
                  means_deaths*option1,
                  means_deaths*option2,
                  means_deaths*option1*state_inc_1[,3+cohorts],
                  means_deaths*option2*state_inc_1[,3+cohorts],
                  alive_sums*max1,
                  alive_sums*max2,
                  alive_sums*max3,
                  alive_sums*max4,
                  alive_sums*max1*state_inc_1[,3+cohorts],
                  alive_sums*max2*state_inc_1[,3+cohorts],
                  alive_sums*max3*state_inc_1[,3+cohorts],
                  alive_sums*max4*state_inc_1[,3+cohorts],
                  alive_sums*state_inc_1[,3+cohorts],
                  alive_sums*state_inc_1[,2],
                  alive_sums*state_inc_1[,2]^2,
                  alive_sums*state_inc_1[,2]^3, 
                  alive_sums*state_inc_1[,2]*state_inc_1[,3+cohorts],
                  alive_sums*state_inc_1[,2]^2*state_inc_1[,3+cohorts],
                  alive_sums*state_inc_1[,3+cohorts]^2,                 
                  alive_sums*option1,
                  alive_sums*option2,
                  alive_sums*option1*state_inc_1[,3+cohorts],
                  alive_sums*option2*state_inc_1[,3+cohorts])
    return(output)
  }}

#RUN MODEL

library(doParallel)

cl <- makeCluster(7)

registerDoParallel(cl)

T = 6
eta = 0.06
p=0.995
N_reg=0.5*10^3
N_VaR=10^4


#mu and sigma and starting value of stock F, mu and sigma and starting value of stock Y
#, correlation between brownian motions and proportion held last.

param_stock = c(0.03, 0.1,100,0.03,0.1,100,0.4, 1)
#guarenteed death and survival benefits, plus terminal time
param_benifits = c(100,110,T)
#parameters for cohorts of form N_0, and p- vectors. In this case, T=2.
param_N0 = c(1000,1000,1000,1000)

param_pvectors = cbind(life_M90(50:(50+T)),life_M90(60:(60+T)),life_M90(70:(70+T)),life_M90(80:(80+T)))

LIFE_param = list(param_stock,param_benifits,param_N0,param_pvectors)

state_0 = c(0,LIFE_param[[1]][3],LIFE_param[[3]],LIFE_param[[1]][6],0)

#First Running of code starts here:
ptm <- proc.time()
valuation_list_life_201204 = LSM_COC_Valuation(N_reg, N_VaR,eta,p,T,state_0,LIFE_param, 
                                               sim_state_LIFE, sim_state_LIFE_unconditional,
                                               square_basis_vector_LIFE)


valuation_list_life_201204[[2]]

proc.time() - ptm

save.image()

thetamatrix_LIFE_201204  = valuation_list_life_201204[[3]]
  
phimatrix_LIFE_201204 = valuation_list_life_201204[[4]]

N_val = 10^3
N_VaR_val = 0.5*10^5

validation_list_life_201204<-LSM_COC_Validation(N_val, N_VaR_val,eta,p,T,state_0,LIFE_param,
                                                thetamatrix_LIFE_201204, phimatrix_LIFE_201204, 
                                                sim_state_LIFE, sim_state_LIFE_unconditional,
                                                square_basis_vector_LIFE)

proc.time() - ptm


save.image()


val_return_matrix = validation_list_life_201204[[1]]

val_prob_matrix = validation_list_life_201204[[2]]

val_VaR_matrix = validation_list_life_201204[[3]]

val_Eval_matrix = validation_list_life_201204[[4]]

val_VaR_est_matrix = validation_list_life_201204[[5]]

val_Eval_est_matrix = validation_list_life_201204[[6]]

s=1

hist(100*(val_return_matrix[s,]-1),100,xlab="Return (%)",main = paste("Excess return at time" , s),xlim=c(-10,25))

hist(100*(1.06*(val_Eval_matrix[s,]+1000)/val_Eval_est_matrix[s,]-1),50,xlab="Return (%)",main = paste("Excess return at time" , s),xlim=c(-10,25))

mean(100*(1.06*(val_Eval_matrix[s,] )/val_Eval_est_matrix[s,]-1))

mean(100*(val_return_matrix[s,]-1))

hist(1-val_prob_matrix[s,],100,xlab="Probability",main = paste("Probability of ruin at time" , s),xlim=c(0.002,0.008))

mean(1-val_prob_matrix[s,],1000,xlab="Probability",main = paste("Probability of ruin at time" , s))

sqrt(mean((val_Eval_matrix[s,]-val_Eval_est_matrix[s,])^2))/sqrt(mean((val_Eval_matrix[s,])^2))

hist(val_Eval_matrix)

max(val_Eval_est_matrix)

hist(val_VaR_est_matrix)


bfval = bruteforce_valuation(10000,eta,p,state_0,LIFE_param,beta =square_basis_vector_LIFE(state_0)*0,
                            basis_function = square_basis_vector_LIFE, sim_state_LIFE, regression_test = FALSE , time = 1)



valuation_list_life_200902[[1]]/(bfval[[1]][1] - 1/(1+eta)*bfval[[1]][2])


bfval_regtest = bruteforce_valuation(100000,eta,p,state_0,LIFE_param,beta =square_basis_vector_LIFE(state_0)*0,
                               basis_function = square_basis_vector_LIFE, sim_state_LIFE, regression_test = TRUE , time = 1)






















