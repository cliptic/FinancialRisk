Heston_upout_barrier_call <- function(vector,S0, K, B,t, r, nbre_step, nbre_sim)
{ 
  v0 = vector[1]
  vbar = vector[2]
  rho = vector[3]
  lambda = vector[4]
  eta = vector[5]
  dt=t/nbre_step
  S=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  v=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  vect_of_ST = c()
  vect_of_payoffs = c()
  under_barrier = TRUE
  count= 0
  for (j in 1:nbre_sim)
  {
    S[j,1]=S0
    v[j,1]=v0
    for (i in 1:nbre_step) 
    {
      
      Z1 = rnorm(1, mean=0, sd=1)  # rnorm(n, mean, sd) generates random numbers from the normal distribution with mean parameter mu and standard deviation sd
      Z2 = rnorm(1, mean=0, sd=1)  # rnorm(n, mean, sd) generates random numbers from the normal distribution with mean parameter mu and standard deviation sd
      
      X1 = Z1
      X2= rho * Z1 + sqrt(1-rho^2) * Z2
      
      
      S[j,i+1] = S[j,i] +  r * S[j,i] * dt + S[j,i] * sqrt(v[j,i]) * sqrt(dt) * X1
      if (S[j,i+1]>B)
      {
        under_barrier = FALSE
        break
      }
      v[j,i+1] = abs(v[j,i] - lambda *( v[j,i] - vbar) * dt+ eta * sqrt(v[j,i]) * sqrt(dt) * X2)
      
    }
    if (under_barrier){count = count+1
    #vect_of_v[count] = v[j,nbre_step+1]
    vect_of_ST[count]=S[j,nbre_step+1]
    vect_of_payoffs[count] = max(vect_of_ST[count] - K,0)
    }else{under_barrier = TRUE}
    
  }  
  
  DF = exp(- r  * t)
  Heston_MC_Call = DF * 1/nbre_sim * sum(vect_of_payoffs)
  return(Heston_MC_Call)}

Heston_upin_barrier_call <- function(vector, S0, K, B,t, r, nbre_step, nbre_sim)
{ 
  v0 = vector[1]
  vbar = vector[2]
  rho = vector[3]
  lambda = vector[4]
  eta = vector[5]
  dt=t/nbre_step
  S=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  v=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  vect_of_ST = c()
  vect_of_payoffs = c()
  count= 0
  for (j in 1:nbre_sim)
  {
    S[j,1]=S0
    v[j,1]=v0
    touch = FALSE
    for (i in 1:nbre_step) 
    {
      
      Z1 = rnorm(1, mean=0, sd=1)  # rnorm(n, mean, sd) generates random numbers from the normal distribution with mean parameter mu and standard deviation sd
      Z2 = rnorm(1, mean=0, sd=1)  # rnorm(n, mean, sd) generates random numbers from the normal distribution with mean parameter mu and standard deviation sd
      
      X1 = Z1
      X2= rho * Z1 + sqrt(1-rho^2) * Z2
      
      
      S[j,i+1] = S[j,i] +  r * S[j,i] * dt + S[j,i] * sqrt(v[j,i]) * sqrt(dt) * X1
      if (S[j,i+1]>B)
      {
        touch = TRUE
      }
      v[j,i+1] = abs(v[j,i] - lambda *( v[j,i] - vbar) * dt+ eta * sqrt(v[j,i]) * sqrt(dt) * X2)
      
    }
    if (touch){count = count+1
    #vect_of_v[count] = v[j,nbre_step+1]
    vect_of_ST[count]=S[j,nbre_step+1]
    vect_of_payoffs[count] = max(vect_of_ST[count] - K,0)
    }
  }  
  
  DF = exp(- r  * t)
  Heston_MC_Call_Barrier_In = DF * 1/nbre_sim * sum(vect_of_payoffs)
  return(Heston_MC_Call_Barrier_In)
}

Heston_upout_delta <- function(vector, S0, K, B,t, r, delS, nbre_step, nbre_sim)
{ 
  v0 = vector[1]
  vbar = vector[2]
  rho = vector[3]
  lambda = vector[4]
  eta = vector[5]
  dt=t/nbre_step
  S1=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  S2=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  v=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  #vect_of_ST = c()
  vect_of_payoffs1 = c()
  vect_of_payoffs2 = c()
  existance1 = TRUE
  existance2 = TRUE
  count1= 0
  count2 = 0
  for (j in 1:nbre_sim)
  {
    S1[j,1]=S0
    S2[j,1]=S0+delS
    v[j,1]=v0
    for (i in 1:nbre_step) 
    {
      
      Z1 = rnorm(1, mean=0, sd=1)  # rnorm(n, mean, sd) generates random numbers from the normal distribution with mean parameter mu and standard deviation sd
      Z2 = rnorm(1, mean=0, sd=1)  # rnorm(n, mean, sd) generates random numbers from the normal distribution with mean parameter mu and standard deviation sd
      
      X1 = Z1
      X2= rho * Z1 + sqrt(1-rho^2) * Z2
      
      
      S1[j,i+1] = S1[j,i] +  r * S1[j,i] * dt + S1[j,i] * sqrt(v[j,i]) * sqrt(dt) * X1
      S2[j,i+1] = S2[j,i] +  r * S2[j,i] * dt + S2[j,i] * sqrt(v[j,i]) * sqrt(dt) * X1
      if (S1[j,i+1]>B)
      {
        existance1 = FALSE
      }
      if (S2[j,i+1]>B)
      {
        existance2 = FALSE
      }
      v[j,i+1] = abs(v[j,i] - lambda *( v[j,i] - vbar) * dt+ eta * sqrt(v[j,i]) * sqrt(dt) * X2)
      if(!(existance1|existance2))
      {break}
    }
    if (existance1){count1 = count1+1
    #vect_of_v[count] = v[j,nbre_step+1]
    #vect_of_ST[count]=S[j,nbre_step+1]
    vect_of_payoffs1[count1] = max(S1[j,nbre_step+1] - K,0)
    }else{existance1 = TRUE}
    if (existance2){count2 = count2+1
    #vect_of_v[count] = v[j,nbre_step+1]
    #vect_of_ST[count]=S[j,nbre_step+1]
    vect_of_payoffs2[count2] = max(S2[j,nbre_step+1] - K,0)
    }else{existance2 = TRUE}
    
  }  
  
  DF = exp(- r  * t)
  Heston_1 = DF * 1/nbre_sim * sum(vect_of_payoffs1)
  Heston_2 = DF * 1/nbre_sim * sum(vect_of_payoffs2)
  Delta = (Heston_2-Heston_1)/delS
  return(Delta)
}

Heston_upout_vega <- function(vector,S0, K, B, t,r, delv0, nbre_step, nbre_sim)
{ 
  v0 = vector[1]
  vbar = vector[2]
  rho = vector[3]
  lambda = vector[4]
  eta = vector[5]
  v0_1 = v0 + delv0
  dt=t/nbre_step
  S1=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  S2=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  v1=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  v2=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  #vect_of_ST = c()
  vect_of_payoffs1 = c()
  vect_of_payoffs2 = c()
  existance1 = TRUE
  existance2 = TRUE
  count1= 0
  count2 = 0
  for (j in 1:nbre_sim)
  {
    S1[j,1]=S0
    S2[j,1]=S0
    v1[j,1]=v0
    v2[j,1]=v0+delv0
    for (i in 1:nbre_step) 
    {
      
      Z1 = rnorm(1, mean=0, sd=1)  # rnorm(n, mean, sd) generates random numbers from the normal distribution with mean parameter mu and standard deviation sd
      Z2 = rnorm(1, mean=0, sd=1)  # rnorm(n, mean, sd) generates random numbers from the normal distribution with mean parameter mu and standard deviation sd
      
      X1 = Z1
      X2= rho * Z1 + sqrt(1-rho^2) * Z2
      
      
      S1[j,i+1] = S1[j,i] +  r * S1[j,i] * dt + S1[j,i] * sqrt(v1[j,i]) * sqrt(dt) * X1
      S2[j,i+1] = S2[j,i] +  r * S2[j,i] * dt + S2[j,i] * sqrt(v2[j,i]) * sqrt(dt) * X1
      if (S1[j,i+1]>B)
      {
        existance1 = FALSE
      }
      if (S2[j,i+1]>B)
      {
        existance2 = FALSE
      }
      v1[j,i+1] = abs(v1[j,i] - lambda *( v1[j,i] - vbar) * dt+ eta * sqrt(v1[j,i]) * sqrt(dt) * X2)
      v2[j,i+1] = abs(v2[j,i] - lambda *( v2[j,i] - vbar) * dt+ eta * sqrt(v2[j,i]) * sqrt(dt) * X2)
      if(!(existance1|existance2))
      {break}
    }
    if (existance1){count1 = count1+1
    #vect_of_v[count] = v[j,nbre_step+1]
    #vect_of_ST[count]=S[j,nbre_step+1]
    vect_of_payoffs1[count1] = max(S1[j,nbre_step+1] - K,0)
    }else{existance1 = TRUE}
    if (existance2){count2 = count2+1
    #vect_of_v[count] = v[j,nbre_step+1]
    #vect_of_ST[count]=S[j,nbre_step+1]
    vect_of_payoffs2[count2] = max(S2[j,nbre_step+1] - K,0)
    }else{existance2 = TRUE}
  }  
  
  DF = exp(- r  * t)
  Heston_1 = DF * 1/nbre_sim * sum(vect_of_payoffs1)
  Heston_2 = DF * 1/nbre_sim * sum(vect_of_payoffs2)
  Vega = (Heston_2-Heston_1)/delv0
  return(Vega)
}

Heston_upout_gamma <- function(vector, S0, K, B,t, r, delS, nbre_step, nbre_sim)
{ 
  v0 = vector[1]
  vbar = vector[2]
  rho = vector[3]
  lambda = vector[4]
  eta = vector[5]
  dt=t/nbre_step
  S1=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  S2=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  S3=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  v=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  #vect_of_ST = c()
  vect_of_payoffs1 = c()
  vect_of_payoffs2 = c()
  vect_of_payoffs3 = c()
  existance1 = TRUE
  existance2 = TRUE
  existance3 = TRUE
  count1= 0
  count2 = 0
  count3 = 0
  for (j in 1:nbre_sim)
  {
    S1[j,1]=S0-delS
    S2[j,1]=S0
    S3[j,1]=S0+delS
    v[j,1]=v0
    for (i in 1:nbre_step) 
    {
      
      Z1 = rnorm(1, mean=0, sd=1)  # rnorm(n, mean, sd) generates random numbers from the normal distribution with mean parameter mu and standard deviation sd
      Z2 = rnorm(1, mean=0, sd=1)  # rnorm(n, mean, sd) generates random numbers from the normal distribution with mean parameter mu and standard deviation sd
      
      X1 = Z1
      X2= rho * Z1 + sqrt(1-rho^2) * Z2
      
      
      S1[j,i+1] = S1[j,i] +  r * S1[j,i] * dt + S1[j,i] * sqrt(v[j,i]) * sqrt(dt) * X1
      S2[j,i+1] = S2[j,i] +  r * S2[j,i] * dt + S2[j,i] * sqrt(v[j,i]) * sqrt(dt) * X1
      S3[j,i+1] = S3[j,i] +  r * S3[j,i] * dt + S3[j,i] * sqrt(v[j,i]) * sqrt(dt) * X1
      if (S1[j,i+1]>B)
      {
        existance1 = FALSE
      }
      if (S2[j,i+1]>B)
      {
        existance2 = FALSE
      }
      if (S3[j,i+1]>B)
      {
        existance3 = FALSE
      }
      v[j,i+1] = abs(v[j,i] - lambda *( v[j,i] - vbar) * dt+ eta * sqrt(v[j,i]) * sqrt(dt) * X2)
      if(!(existance1|existance2|existance3))
      {break}
    }
    if (existance1){count1 = count1+1
    #vect_of_v[count] = v[j,nbre_step+1]
    #vect_of_ST[count]=S[j,nbre_step+1]
    vect_of_payoffs1[count1] = max(S1[j,nbre_step+1] - K,0)
    }else{existance1 = TRUE}
    if (existance2){count2 = count2+1
    #vect_of_v[count] = v[j,nbre_step+1]
    #vect_of_ST[count]=S[j,nbre_step+1]
    vect_of_payoffs2[count2] = max(S2[j,nbre_step+1] - K,0)
    }else{existance2 = TRUE}
    if (existance3){count3 = count3+1
    #vect_of_v[count] = v[j,nbre_step+1]
    #vect_of_ST[count]=S[j,nbre_step+1]
    vect_of_payoffs3[count3] = max(S3[j,nbre_step+1] - K,0)
    }else{existance3 = TRUE}
    
  }  
  
  DF = exp(- r  * t)
  H_1 = DF * 1/nbre_sim * sum(vect_of_payoffs1)
  H_2 = DF * 1/nbre_sim * sum(vect_of_payoffs2)
  H_3 = DF * 1/nbre_sim * sum(vect_of_payoffs3)
  Gamma = (H_1+H_3-2*H_2)/delS^2
  return(Gamma)
}

Heston_upout_survival <- function(vector,S0, K, B,t, r, nbre_step, nbre_sim)
{ 
  v0 = vector[1]
  vbar = vector[2]
  rho = vector[3]
  lambda = vector[4]
  eta = vector[5]
  dt=t/nbre_step
  S=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  v=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  under_barrier = TRUE
  count= 0
  for (j in 1:nbre_sim)
  {
    S[j,1]=S0
    v[j,1]=v0
    for (i in 1:nbre_step) 
    {
      
      Z1 = rnorm(1, mean=0, sd=1)  # rnorm(n, mean, sd) generates random numbers from the normal distribution with mean parameter mu and standard deviation sd
      Z2 = rnorm(1, mean=0, sd=1)  # rnorm(n, mean, sd) generates random numbers from the normal distribution with mean parameter mu and standard deviation sd
      
      X1 = Z1
      X2= rho * Z1 + sqrt(1-rho^2) * Z2
      
      
      S[j,i+1] = S[j,i] +  r * S[j,i] * dt + S[j,i] * sqrt(v[j,i]) * sqrt(dt) * X1
      if (S[j,i+1]>B)
      {
        under_barrier = FALSE
        break
      }
      v[j,i+1] = abs(v[j,i] - lambda *( v[j,i] - vbar) * dt+ eta * sqrt(v[j,i]) * sqrt(dt) * X2)
      
    }
    if (under_barrier){count = count+1
    } else {under_barrier = TRUE}
    
  }  
  
   return(count/nbre_sim)}
