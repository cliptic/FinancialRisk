
### Price ####

BS_MC_UO_Call_Price <- function(S0,K,B,T,sig,r, nbre_step, nbre_sim)
{ 
  dt = T/nbre_step
  S = matrix(nrow = nbre_sim, ncol = nbre_step+1)
  Sd = matrix(nrow = nbre_sim, ncol = nbre_step+1)
  survival = 0
  
  for (j in 1:nbre_sim)
  {
    under_barrier = TRUE
    S[j,1]=S0
    
    for (i in 1:nbre_step) 
    {
      Z = rnorm(1, mean=0, sd=1)  # rnorm(n, mean, sd) generates random numbers from the normal 
      #distribution with mean parameter mean and standard deviation sd
      S[j,i+1] = S[j,i] + r*S[j,i]*dt+sig*sqrt(dt)*S[j,i]*Z
      
      if (S[j,i+1]>B)
      {
        under_barrier = FALSE
        break
      }
      
    }
    if (under_barrier)
    {survival = survival+1} else {under_barrier = TRUE}
    
  }  
  
  VectofST = S[,nbre_step+1]
  vectofsteps = seq(from = 1, to = nbre_step+1, by = 1)
  
  DF = exp(-r*T)
  vectofpayoff = pmax(VectofST - K,0)
  
  BS_MC_Call_Price = DF * sum(vectofpayoff, na.rm = TRUE)/nbre_sim
  BS_MC_Call_survival_probability = survival/nbre_sim
  BSMC_Call = c(BS_MC_Call_Price, BS_MC_Call_survival_probability)
  return(BSMC_Call)
}


#### GREEKS ####

BS_MC_UO_Call_Delta <- function(S0,K,B,T,sig,r, nbre_step, nbre_sim, delS)
{ 
  dt = T/nbre_step
  S1 = matrix(nrow = nbre_sim, ncol = nbre_step+1)
  S2 = matrix(nrow = nbre_sim, ncol = nbre_step+1)
  vect_of_payoffs1 = c()
  vect_of_payoffs2 = c()
  under_barrier1 = TRUE
  under_barrier2 = TRUE
  count1 = 0
  count2 = 0
  sig1 = sig
  sig2 = sig
  
  for (j in 1:nbre_sim)
  {
    under_barrier1 = TRUE
    under_barrier2 = TRUE
    S1[j,1]=S0
    S2[j,1]=(S0+delS)
    
    for (i in 1:nbre_step) 
    {
      Z = rnorm(1, mean=0, sd=1)
      S1[j,i+1] = S1[j,i] + r*S1[j,i]*dt+sig1*sqrt(dt)*S1[j,i]*Z
      S2[j,i+1] = S2[j,i] + r*S2[j,i]*dt+sig2*sqrt(dt)*S2[j,i]*Z
      if (S1[j,i+1]>B)
      {
        under_barrier1 = FALSE
      }
      if (S2[j,i+1]>B)
      {
        under_barrier2 = FALSE
      }
      if(!(under_barrier1|under_barrier2))
      {break}
    }
    if (under_barrier1){count1 = count1+1
    vect_of_payoffs1[count1] = max(S1[j,nbre_step+1] - K,0)
    }
    if (under_barrier2){count2 = count2+1
    vect_of_payoffs2[count2] = max(S2[j,nbre_step+1] - K,0)
    }
  }  
  
  VectofST1 = S1[,nbre_step+1]
  VectofST2 = S2[,nbre_step+1]
  
  vectofsteps = seq(from = 1, to = nbre_step+1, by = 1)
  DF = exp(-r*T)
  
  vectofpayoff1 = pmax(VectofST1 - K,0)
  vectofpayoff2 = pmax(VectofST2 - K,0)
  
  BS_MC_Call_Price1 = DF * sum(vectofpayoff1, na.rm = TRUE)/nbre_sim
  BS_MC_Call_Price2 = DF * sum(vectofpayoff2, na.rm = TRUE)/nbre_sim
  BS_MC_Call_survival_probability1 = count1/nbre_sim
  BS_MC_Call_survival_probability2 = count2/nbre_sim
  BS_MC_Call_Delta = (BS_MC_Call_Price2-BS_MC_Call_Price1)/delS
  BSMC_Call <- c(BS_MC_Call_Price1, BS_MC_Call_Price2, BS_MC_Call_survival_probability1,BS_MC_Call_survival_probability2, BS_MC_Call_Delta)
  return(BSMC_Call)
}

BS_MC_UO_Call_Vega <- function(S0,K,B,T,sig,r, nbre_step, nbre_sim, del_sig)
{ 
  dt = T/nbre_step
  S1 = matrix(nrow = nbre_sim, ncol = nbre_step+1)
  S2 = matrix(nrow = nbre_sim, ncol = nbre_step+1)
  vect_of_payoffs1 = c()
  vect_of_payoffs2 = c()
  under_barrier1 = TRUE
  under_barrier2 = TRUE
  count1 = 0
  count2 = 0
  sig1 = sig
  sig2 = sig+del_sig
  
  for (j in 1:nbre_sim)
  {
    under_barrier1 = TRUE
    under_barrier2 = TRUE
    S1[j,1]=S0
    S2[j,1]=S0
    
    for (i in 1:nbre_step) 
    {
      Z = rnorm(1, mean=0, sd=1)
      S1[j,i+1] = S1[j,i] + r*S1[j,i]*dt+sig1*sqrt(dt)*S1[j,i]*Z
      S2[j,i+1] = S2[j,i] + r*S2[j,i]*dt+sig2*sqrt(dt)*S2[j,i]*Z
      if (S1[j,i+1]>B)
      {
        under_barrier1 = FALSE
      }
      if (S2[j,i+1]>B)
      {
        under_barrier2 = FALSE
      }
      if(!(under_barrier1|under_barrier2))
      {break}
    }
    if (under_barrier1){count1 = count1+1
    vect_of_payoffs1[count1] = max(S1[j,nbre_step+1] - K,0)
    }
    if (under_barrier2){count2 = count2+1
    vect_of_payoffs2[count2] = max(S2[j,nbre_step+1] - K,0)
    }
  }  
  
  VectofST1 = S1[,nbre_step+1]
  VectofST2 = S2[,nbre_step+1]
  
  vectofsteps = seq(from = 1, to = nbre_step+1, by = 1)
  DF = exp(-r*T)
  
  vectofpayoff1 = pmax(VectofST1 - K,0)
  vectofpayoff2 = pmax(VectofST2 - K,0)
  
  BS_MC_Call_Price1 = DF * sum(vectofpayoff1, na.rm = TRUE)/nbre_sim
  BS_MC_Call_Price2 = DF * sum(vectofpayoff2, na.rm = TRUE)/nbre_sim
  BS_MC_Call_survival_probability1 = count1/nbre_sim
  BS_MC_Call_survival_probability2 = count2/nbre_sim
  BS_MC_Call_Vega = (BS_MC_Call_Price1-BS_MC_Call_Price2)/del_sig
  BSMC_Call <- c(BS_MC_Call_Price1, BS_MC_Call_Price2, BS_MC_Call_survival_probability1,BS_MC_Call_survival_probability2, BS_MC_Call_Vega)
  return(BSMC_Call)
}

BS_MC_UO_Call_Gamma <- function(S0,K,B,T,sig,r, nbre_step, nbre_sim, delS)
{ 
  dt = T/nbre_step
  S1 = matrix(nrow = nbre_sim, ncol = nbre_step+1)
  S2 = matrix(nrow = nbre_sim, ncol = nbre_step+1)
  S3 = matrix(nrow = nbre_sim, ncol = nbre_step+1)
  vect_of_payoffs1 = c()
  vect_of_payoffs2 = c()
  vect_of_payoffs3 = c()
  under_barrier1 = TRUE
  under_barrier2 = TRUE
  under_barrier3 = TRUE
  count1 = 0
  count2 = 0
  count3 = 0
  
  for (j in 1:nbre_sim)
  {
    under_barrier1 = TRUE
    under_barrier2 = TRUE
    under_barrier3 = TRUE
    S1[j,1]=S0
    S2[j,1]=(S0+delS)
    S3[j,1]=(S0-delS)
    
    for (i in 1:nbre_step) 
    {
      Z = rnorm(1, mean=0, sd=1)
      
      # rnorm(n, mean, sd) generates random numbers from the normal 
      #distribution with mean parameter mean and standard deviation sd
      S1[j,i+1] = S1[j,i] + r*S1[j,i]*dt+sig*sqrt(dt)*S1[j,i]*Z
      S2[j,i+1] = S2[j,i] + r*S2[j,i]*dt+sig*sqrt(dt)*S2[j,i]*Z
      S3[j,i+1] = S3[j,i] + r*S3[j,i]*dt+sig*sqrt(dt)*S3[j,i]*Z
      
      if (S1[j,i+1]>B)
      {
        under_barrier1 = FALSE
      }
      if (S2[j,i+1]>B)
      {
        under_barrier2 = FALSE
      }
      if (S3[j,i+1]>B)
      {
        under_barrier3 = FALSE
      }
      if(!(under_barrier1|under_barrier2|under_barrier3))
      {break}
    }
    if (under_barrier1){count1 = count1+1
    vect_of_payoffs1[count1] = max(S1[j,nbre_step+1] - K,0)
    }
    if (under_barrier2){count2 = count2+1
    vect_of_payoffs2[count2] = max(S2[j,nbre_step+1] - K,0)
    }
    if (under_barrier3){count3 = count3+1
    vect_of_payoffs3[count3] = max(S3[j,nbre_step+1] - K,0)
    }
  }  
  
  VectofST1 = S1[,nbre_step+1]
  VectofST2 = S2[,nbre_step+1]
  VectofST3 = S3[,nbre_step+1]
  
  vectofsteps = seq(from = 1, to = nbre_step+1, by = 1)
  DF = exp(-r*T)
  
  vectofpayoff1 = pmax(VectofST1 - K,0)
  vectofpayoff2 = pmax(VectofST2 - K,0)
  vectofpayoff3 = pmax(VectofST3 - K,0)
  
  BS_MC_Call_Price1 = DF * sum(vectofpayoff1, na.rm = TRUE)/nbre_sim
  BS_MC_Call_Price2 = DF * sum(vectofpayoff2, na.rm = TRUE)/nbre_sim
  BS_MC_Call_Price3 = DF * sum(vectofpayoff3, na.rm = TRUE)/nbre_sim
  BS_MC_Call_survival_probability1 = count1/nbre_sim
  BS_MC_Call_survival_probability2 = count2/nbre_sim
  BS_MC_Call_survival_probability3 = count3/nbre_sim
  BS_MC_Call_Gamma = (BS_MC_Call_Price1+BS_MC_Call_Price3-2*BS_MC_Call_Price2)/delS^2
  BSMC_Call <- c(BS_MC_Call_Price1, BS_MC_Call_Price2, BS_MC_Call_Price3, BS_MC_Call_survival_probability1,BS_MC_Call_survival_probability2, BS_MC_Call_survival_probability3, BS_MC_Call_Gamma)
  return(BSMC_Call)
}


