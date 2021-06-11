#component functions


VaR_simfn <- function(state,param, p,beta,N_VaR,Evalue_calc = FALSE,sim_state,
                            basis_function,time=0,lasttime=FALSE){
  #matrix of future states, row-wise
  state_next = sim_state(state,param,N_VaR,time)
  
  #X_next =X_value(state,state_next,X_parameters)
  #V_next = square_basis_vector(state_next)%*%beta
    if(lasttime==FALSE){
        values = state_next[,1]+basis_function(state_next)%*%beta
    }else{
      values = state_next[,1]
      }
  

  VaR = quantile(values,p)
  if(Evalue_calc == TRUE){
    diff = VaR-values
    diff[(diff<0)]=0
    Evalue = mean(diff)
    return(c(VaR,Evalue))
  }
  return(VaR)
}



VaR_Value_regressionfn <- function(state_0,param,beta,p,T,N_reg=10^3,N_VaR=10^4,time=0,sim_state,
                                         sim_state_uncond,VaR_simfn,basis_function){
  
  if(time != 0){
    
    #simulate states
    
    states = sim_state_uncond(time,state_0,sim_state,param,N_reg)
    #N_reg =10000
    VaR_vector = seq(1,N_reg)
    Eval_vector = seq(1,N_reg)
    
    # for(i in seq(1,N_reg)){
    #   #i=1
    #   
    #   VaRsim=VaR_simfn(states[i,],GARCH_param, p,beta,N_VaR,TRUE,sim_state,basis_function)
    #   #VaRsim=VaR_simfn(states[i,],GARCH_param, p,beta,N_VaR,TRUE,sim_state,square_basis_vector)
    #   #VaRsim=VaR_simfn(states[i,],GARCH_param, p,beta,N_VaR,TRUE,sim_state,dummy_basis_vector
    #   VaR_vector[i] = VaRsim[1]
    #   Eval_vector[i] = VaRsim[2]
    # }
    
    #parallell programming code
    
    lasttime = (time==T-1)
    
    VaRsim <- foreach(i = 1:N_reg, .combine = 'cbind',.export=ls(envir=globalenv())) %dopar% {
      
      VaR_simfn(states[i,],param, p,beta,N_VaR,TRUE,sim_state,basis_function,time,lasttime)
    }
    
    
    VaR_vector = VaRsim[1,]
    Eval_vector = VaRsim[2,]
    
    #remove the column of ones from the basis function matrix!
    basis = basis_function(states)[,-1]
    
    modl_VaR <- lm(VaR_vector ~ basis)
    
    summary(modl_VaR)
    
    theta = modl_VaR[[1]]
    
    
    
    modl_Evalue <- lm(Eval_vector ~ basis)
    
    #summary(modl_Evalue)
    
    #hist(Eval_vector,50)
    #hist(modl_Evalue$fitted.values,50)
    
    phi = modl_Evalue[[1]]
    
    #predd = cbind(1,basis)%*%phi
    #hist(predd,50)
    parameters = cbind(theta, phi)
    return(parameters)
  }else{
    VaRsim=VaR_simfn(state_0,param, p,beta,N_VaR,TRUE,sim_state,basis_function,0,FALSE)
    VaR_0 = VaRsim[1]
    Eval_0 = VaRsim[2]
    return(c(VaR_0,Eval_0))
  }
}


VaR_Value_validationfn=function(eta,p,VaR,EVal,N_VaR,state,beta,param,sim_state,
                                      basis_function,time=0,terminaltime =2){
  
  #Recall V=VaR-(1/(1+eta))*EVal ==> C=(1/(1+eta))*EVal
  
  state_next = sim_state(state,param,N_VaR,time)
  values =  state_next[,1]
  if(time != terminaltime -1){
    values = values + basis_function(state_next)%*%beta
  }
  
  prob=ecdf(values)(VaR)
  VaR_part = quantile(values,p)
  values=VaR_part-values
  values[values<0]=0
  mean_part=mean(values)
  return_of_investment=(1+eta)*mean_part/EVal
  
  
  return(c(prob,return_of_investment,VaR_part,mean_part))
}




#aggregation functions that run valuation and validation

LSM_COC_Valuation = function(N_reg, N_VaR,eta,p,T,state_0,param, sim_state, sim_state_uncond,basis_function){
  
  
  number_of_basis_functions=length(basis_function(state_0))
  thetamatrix = matrix(rep(0,number_of_basis_functions*T),ncol=T)
  phimatrix = matrix(rep(0,number_of_basis_functions*T),ncol=T)
  
  for(i in seq(T-1,1)){
    print(i)
    theta = thetamatrix[,i+1]
    phi = phimatrix[,i+1]
    
    #parameters describing the value - linear combination of theta and phi
    
    beta = theta-(1/(1+eta))*phi
    print(beta)
    
    thetaphi = VaR_Value_regressionfn(state_0,param,beta,p,T,N_reg,N_VaR,time=i,
                                           sim_state,sim_state_uncond,VaR_simfn,basis_function)

    thetamatrix[,i]=thetaphi[1:number_of_basis_functions]
    phimatrix[,i]=thetaphi[(number_of_basis_functions+1):(2*number_of_basis_functions)]
    thetamatrix[,i][is.na(thetamatrix[,i])] <- 0
    phimatrix[,i][is.na(phimatrix[,i])] <- 0
  }
  
  #Calculating time zero R and V
  
  theta = thetamatrix[,1]
  phi = phimatrix[,1]
  
  beta = theta - (1/(1+eta))*phi
  
  Capreq_Eval = VaR_Value_regressionfn(state_0,param,beta,p,T,N_reg,10*N_VaR,time=0,
                                            sim_state,sim_state_unconditional,VaR_simfn,basis_function)
  
  Value = Capreq_Eval[1]-(1/(1+eta))*Capreq_Eval[2]
  
  
  return(list(Value, Capreq_Eval, thetamatrix,phimatrix))}


LSM_COC_Validation = function(N_val, N_VaR_val,eta,p,T,state_0,param,thetamatrix,
                              phimatrix, sim_state, sim_state_uncond,basis_function){
  
  
  val_prob_matrix =matrix(0,T-1,N_val)
  val_return_matrix =matrix(0,T-1,N_val)
  val_VaR_matrix =matrix(0,T-1,N_val)
  val_Eval_matrix =matrix(0,T-1,N_val)
  val_VaR_est_matrix =matrix(0,T-1,N_val)
  val_Eval_est_matrix =matrix(0,T-1,N_val)
  
  
  for(time in seq(1,(T-1))){
    
    theta=thetamatrix[,time]
    phi=phimatrix[,time]
    theta_next=thetamatrix[,time+1]
    phi_next=phimatrix[,time+1]
    
    
    beta_next =theta_next - (1/(1+eta))*phi_next
    
    states=sim_state_uncond(time,state_0,sim_state,param,N_val)
    
    VaR_vector_est = basis_function(states)%*%theta
    E_vector_est = basis_function(states)%*%phi
    
    val_VaR_est_matrix[time,]=VaR_vector_est
    val_Eval_est_matrix[time,]=E_vector_est
    
    #N_val=10
    val_returns = rep(0,N_val)
    val_probs = rep(0,N_val)
    val_VaR = rep(0,N_val)
    val_Eval = rep(0,N_val)
    
    # for(i in seq(1,N_val)){
    #   if(i%%100==0){print(i)}
    #   
    #   
    #   prob_return = VaR_Value_validationfn(eta,p,VaR_vector_est[i],E_vector_est[i],N_VaR_val,states[i,],beta_next,param
    #                                              ,sim_state,basis_function,time)
    #   val_probs[i]=prob_return[1]
    #   val_returns[i]=prob_return[2]
    #   val_VaR[i] = prob_return[3]
    #   val_Eval[i] = prob_return[4]
    #   
    #   
    # }
    
  val_matrix <- foreach(i = seq(1,N_val), .combine = cbind,.export=ls(envir=globalenv())) %dopar% {
        VaR_Value_validationfn(eta,p,VaR_vector_est[i],E_vector_est[i],N_VaR_val,states[i,],beta_next,param
                                           ,sim_state,basis_function,time,T)
    }
    val_return_matrix[time,]=val_matrix[2,]
    val_prob_matrix[time,]=val_matrix[1,]
    val_VaR_matrix[time,]=val_matrix[3,]
    val_Eval_matrix[time,]=val_matrix[4,]
    
    # val_return_matrix[time,]=val_returns
    # val_prob_matrix[time,]=val_probs
    # val_VaR_matrix[time,]=val_VaR
    # val_Eval_matrix[time,]=val_Eval
    
  }
  
  
  
  return(list(val_return_matrix,val_prob_matrix,val_VaR_matrix,val_Eval_matrix,val_VaR_est_matrix,val_Eval_est_matrix))
  
  
  }


#NOTE: State_0 may be a matrix of states!!!

bruteforce_valuation <- function(N0,N1,eta,p,state_0,param, beta=0, 
                                 basis_function, sim_state, regression_test = FALSE, 
                                 time =0,lasttime=TRUE){
  states_time1 = sim_state(state_0,param,N0,time)
  #values =  state_next[,1]
  VaRsim <- foreach(i = 1:N0, .combine = 'cbind',.export=ls(envir=globalenv())) %dopar% {
    VaR_simfn(states_time1[i,],param, p,beta,N1,TRUE,sim_state,basis_function,time+1,lasttime)
  }
  if(regression_test==TRUE){
    return(list(states_time1,VaRsim))
  } else{
    values = states_time1[,1]+VaRsim[1,] - (1/(1+eta))*VaRsim[2,]
    VaR = quantile(values,p)
    diff = VaR-values
    diff[(diff<0)]=0
    Evalue = mean(diff)
    return(list(c(VaR,Evalue),VaRsim))
  }
}








################################################################################################################
##################  UNDER CONSTRUCTION  ########################################################################
################################################################################################################

VaR_Value_regressionfn_fitReLucovariates <- function(state_0,param,beta,p,T,N_reg=10^3,N_VaR=10^4,time=0,
                                                     sim_state,sim_state_uncond,VaR_simfn,basis_function){
  
  #this function assumes existence of global matrix variable K
  
  if(time != 0){
    
    #simulate states
    
    states = sim_state_uncond(time,state_0,sim_state,param,N_reg)
    #N_reg =10000
    VaR_vector = seq(1,N_reg)
    Eval_vector = seq(1,N_reg)
    
    # for(i in seq(1,N_reg)){
    #   #i=1
    #   
    #   VaRsim=VaR_simfn(states[i,],GARCH_param, p,beta,N_VaR,TRUE,sim_state,basis_function)
    #   #VaRsim=VaR_simfn(states[i,],GARCH_param, p,beta,N_VaR,TRUE,sim_state,square_basis_vector)
    #   #VaRsim=VaR_simfn(states[i,],GARCH_param, p,beta,N_VaR,TRUE,sim_state,dummy_basis_vector
    #   VaR_vector[i] = VaRsim[1]
    #   Eval_vector[i] = VaRsim[2]
    # }
    
    #parallell programming code
    
    lasttime = (time==T-1)
    
    VaRsim <- foreach(i = 1:N_reg, .combine = 'cbind',.export=ls(envir=globalenv())) %dopar% {
      VaR_simfn(states[i,],param, p,beta,N_VaR,TRUE,sim_state,basis_function,time,lasttime)
    }
    
    
    VaR_vector = VaRsim[1,]
    Eval_vector = VaRsim[2,]
    
    
    # Code to fit relu covariates! 
    
    #V_vector = VaR_vector - 1/(1+eta)*Eval_vector
    
    basis = basis_function(states)[,-1]
    modl_preliminary <- lm(Eval_vector ~ basis)
    resids = residuals(modl_preliminary)
    resids2 = resids
    #Find 4 best K - values
    
    relu_values = seq(50,250,1)
    
    test_covariates  = mat.or.vec(N_reg,length(relu_values))
    
    for( i in 1:length(relu_values)){
      test_covariates[,i]  = sapply(states[,2],max,relu_values[i])
    }
    indices = rep(1,4)
    
    for(j in 1:4){
      
      rsqvec = 0*(1:length(relu_values))
      for( i in 1:length(relu_values)){
        
      temp_modl <- lm(resids2 ~  test_covariates[,i])

      rsqvec[i] =summary(temp_modl)$r.squared
      
      }      

      indices[j] =which(rsqvec == max(rsqvec))
      K[time,j] <<- relu_values[indices[j]]

      temp_modl2 <- lm(resids ~  test_covariates[,indices])
      
      resids2 = residuals(temp_modl2)
    }

    #remove the column of ones from the basis function matrix!
    basis = basis_function(states)[,-1]
    
    modl_VaR <- lm(VaR_vector ~ basis)
    
    summary(modl_VaR)
    #summary(modl_VaR)
    
    theta = modl_VaR[[1]]
    
    
    
    modl_Evalue <- lm(Eval_vector ~ basis)
    
    #summary(modl_Evalue)
    
    #hist(Eval_vector,50)
    #hist(modl_Evalue$fitted.values,50)
    
    phi = modl_Evalue[[1]]
    
    #predd = cbind(1,basis)%*%phi
    #hist(predd,50)
    parameters = cbind(theta, phi)
    return(parameters)
  }else{
    VaRsim=VaR_simfn(state_0,param, p,beta,N_VaR,TRUE,sim_state,basis_function,0,FALSE)
    VaR_0 = VaRsim[1]
    Eval_0 = VaRsim[2]
    return(c(VaR_0,Eval_0))
  }
}

LSM_COC_Valuation_fitReLucovariates = function(N_reg, N_VaR,eta,p,T,state_0,param, sim_state, sim_state_uncond,
                                               basis_function){
  
  #Initialize K matrix
  K<<-mat.or.vec(T-1,4)
  
  number_of_basis_functions=length(basis_function(state_0))
  thetamatrix = matrix(rep(0,number_of_basis_functions*T),ncol=T)
  phimatrix = matrix(rep(0,number_of_basis_functions*T),ncol=T)
  
  for(i in seq(T-1,1)){
    print(i)
    
    theta = thetamatrix[,i+1]
    phi = phimatrix[,i+1]
    
    #parameters describing the value - linear combination of theta and phi
    
    beta = theta-(1/(1+eta))*phi
    print(beta)
    
    thetaphi = VaR_Value_regressionfn_fitReLucovariates(state_0,param,beta,p,T,N_reg,N_VaR,time=i,
                                      sim_state,sim_state_uncond,VaR_simfn,basis_function)
    
    thetamatrix[,i]=thetaphi[1:number_of_basis_functions]
    phimatrix[,i]=thetaphi[(number_of_basis_functions+1):(2*number_of_basis_functions)]
    thetamatrix[,i][is.na(thetamatrix[,i])] <- 0
    phimatrix[,i][is.na(phimatrix[,i])] <- 0
  }
  
  #Calculating time zero R and V
  
  theta = thetamatrix[,1]
  phi = phimatrix[,1]
  
  beta = theta - (1/(1+eta))*phi
  
  Capreq_Eval = VaR_Value_regressionfn(state_0,param,beta,p,T,N_reg,10*N_VaR,time=0,
                                       sim_state,sim_state_unconditional,VaR_simfn,basis_function)
  
  Value = Capreq_Eval[1]-(1/(1+eta))*Capreq_Eval[2]
  
  
  return(list(Value, Capreq_Eval, thetamatrix,phimatrix))}

square_basis_vector_LIFE_fitReLucovariates<- function(state){
  
  cohorts = length(LIFE_param[[3]])
  
  # Empirically best K - values for time t= 5
  
  if(is.matrix(state)==FALSE){
    time=1
  }else{
    time = tail(state[1,], n=1)
  }
  
  K_1 =K[time,1]
  K_2 =K[time,1]
  K_3 =K[time,1]
  K_4 =K[time,1]
  
  #Note: depends on the global matrix variable K! 
  
  
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













