## 3 ##

## Multivariate

#AAPL
S10 = 175
#sig1 = 0.241204
sig1 = 0.2421036
# vector1 = c(0.046710,  0.144576, -0.274333,  0.950736,  0.532452)
vector1 = c(0.2421036, 0.05033186, 0.12496464, -0.2544616, 0.96821798, 0.49965078)

#FB
S20 = 178
# sig2 = 0.266094
sig2 = 0.2663616
# vector2 = c(0.048647, 0.234835, -0.41129, 0.922655, 0.591163)
vector2 = c(0.2663616, 0.05890884, 0.19316360, -0.3990496, 0.92054186, 0.57198985)
            
  
nbre_step = 500
nbre_sim = 10000
rho = 0.451742
K = 0
r = 0.01

## functions

BS_MC_UO_Call_PriceMultivariate <- function(S10, S20, B, K, T, sig1, sig2, r, nbre_step, nbre_sim, rho, delrho)
{ 
  
  dt=T/nbre_step
  S1 = matrix(nrow= nbre_sim, ncol = nbre_step+1)
  S2 = matrix(nrow= nbre_sim, ncol = nbre_step+1)
  S3 = matrix(nrow= nbre_sim, ncol = nbre_step+1)
  for (j in 1:nbre_sim)
  {
    S1[j,1]=S10
    S2[j,1]=S20
    S3[j,1]=S20
    for (i in 1:nbre_step) 
    {
      Z1 = rnorm(1, mean=0, sd=1)
      Z2 = rnorm(1, mean=0, sd=1)
      X1 = Z1
      X2 = rho*Z1 + sqrt(1-rho^2)*Z2
      X3 = (rho+delrho)*Z1 + sqrt(1-(rho+delrho)^2)*Z2
      # rnorm(n, mean, sd) generates random numbers from the normal 
      #distribution with mean parameter mean and standard deviation sd
      S1[j,i+1] = S1[j,i] + r*S1[j,i]*dt+sig1*sqrt(dt)*S1[j,i]*X1
      S2[j,i+1] = S2[j,i] + r*S2[j,i]*dt+sig2*sqrt(dt)*S2[j,i]*X2
      S3[j,i+1] = S3[j,i] + r*S3[j,i]*dt+sig2*sqrt(dt)*S3[j,i]*X3
    }
  }  
  
  VectofS1T = S1[,nbre_step+1]
  VectofS2T = S2[,nbre_step+1]
  VectofS3T = S3[,nbre_step+1]
  vectofsteps = seq(from = 1, to = nbre_step+1, by = 1)
  DF = exp(-r*T)
  vectofpayoff1 = pmax(VectofS1T - K, VectofS2T - K,0)
  vectofpayoff2 = pmax(VectofS1T - K, VectofS3T - K,0)
  BS_MC_Call_Price1 = DF * sum(vectofpayoff1)/nbre_sim
  BS_MC_Call_Price2 = DF * sum(vectofpayoff2)/nbre_sim
  Correlation_sensitivity1 = (BS_MC_Call_Price2 - BS_MC_Call_Price1)/delrho
  Correlation_sensitivity2 = ((BS_MC_Call_Price2 - BS_MC_Call_Price1)/BS_MC_Call_Price1)/(delrho/rho)
  results = c(BS_MC_Call_Price1, Correlation_sensitivity1, Correlation_sensitivity2)
  return(results)
}


Heston_UO_Call_PriceMultivariate <- function(vector1, vector2, S10, S20, K, t, r, nbre_step, nbre_sim, rho, delrho)
{ 
#### HERE???
  correlation = rho
  delcorrelation = rho + delrho
####
  v10 = vector1[1]
  vbar1 = vector1[2]
  rho1 = vector1[3]
  lambda1 = vector1[4]
  eta1 = vector1[5]
  v20 = vector2[1]
  vbar2 = vector2[2]
  rho2 = vector2[3]
  lambda2 = vector2[4]
  eta2 = vector2[5]
  dt=T/nbre_step
  S1=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  S2=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  S3=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  v1=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  v2=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  v3=matrix(nrow= nbre_sim, ncol= nbre_step+1)
  vect_of_payoffs1 = c()
  vect_of_payoffs2 = c()
  vect_of_payoffs3 = c()
  for (j in 1:nbre_sim)
  {
    S1[j,1]=S10
    S2[j,1]=S20
    S3[j,1]=S20
    v1[j,1]=v10
    v2[j,1]=v20
    v3[j,1]=v20
    for (i in 1:nbre_step) 
    {
      Z1 = rnorm(1, mean=0, sd=1)  # rnorm(n, mean, sd) generates random numbers from the normal distribution with mean parameter mu and standard deviation sd
      Z2 = rnorm(1, mean=0, sd=1)  # rnorm(n, mean, sd) generates random numbers from the normal distribution with mean parameter mu and standard deviation sd
      Z3 = rnorm(1, mean=0, sd=1)
      Z4 = rnorm(1, mean=0, sd=1)
      X1_1 = Z1
      X1_2 = rho1 * Z1 + sqrt(1-rho1^2) * Z2
      X2_1 = correlation * Z1 + sqrt(1-correlation^2) * Z3
      X2_2 = rho2 * Z1 + sqrt(1-rho2^2) * Z4
      X3_1 = delcorrelation * Z1 + sqrt(1-delcorrelation^2) * Z3
      X3_2 = rho2 * Z1 + sqrt(1-rho2^2) * Z4
    
      S1[j,i+1] = S1[j,i] +  r * S1[j,i] * dt + S1[j,i] * sqrt(v1[j,i]) * sqrt(dt) * X1_1
      S2[j,i+1] = S2[j,i] +  r * S2[j,i] * dt + S2[j,i] * sqrt(v2[j,i]) * sqrt(dt) * X2_1
      S3[j,i+1] = S3[j,i] +  r * S3[j,i] * dt + S3[j,i] * sqrt(v3[j,i]) * sqrt(dt) * X3_1

      v1[j,i+1] = abs(v1[j,i] - lambda1 *( v1[j,i] - vbar1) * dt+ eta1 * sqrt(v1[j,i]) * sqrt(dt) * X1_2)
      v2[j,i+1] = abs(v2[j,i] - lambda2 *( v2[j,i] - vbar2) * dt+ eta2 * sqrt(v2[j,i]) * sqrt(dt) * X2_2) 
      v3[j,i+1] = abs(v3[j,i] - lambda2 *( v3[j,i] - vbar2) * dt+ eta2 * sqrt(v3[j,i]) * sqrt(dt) * X3_2) 
    }
  }  
  DF = exp(-r*T)
  vectofsteps = seq(from = 1, to = nbre_step+1, by = 1)  
  VectofS1T = S1[,nbre_step+1]
  VectofS2T = S2[,nbre_step+1]
  VectofS3T = S3[,nbre_step+1]
  vectofpayoff = pmax(VectofS1T - K, VectofS2T - K,0)
  vectofpayoff2 = pmax(VectofS1T - K, VectofS3T - K,0)
  Heston_Price = DF * 1/nbre_sim * sum(vectofpayoff)
  Heston_Price2 = DF * 1/nbre_sim * sum(vectofpayoff2)
  Correlation_sensitivity1 = (Heston_Price2 - Heston_Price)/(delcorrelation-correlation)
  Correlation_sensitivity2 = ((Heston_Price2 - Heston_Price)/Heston_Price)/((delcorrelation-correlation)/correlation)
  results = c(Heston_Price, Correlation_sensitivity1, Correlation_sensitivity2)
  return(results)
}


BSMC_Call_Multivariate = BS_MC_UO_Call_PriceMultivariate(S10, S20, B, K, T, sig1, sig2, r, nbre_step, nbre_sim, rho, rho/1000)
print(BSMC_Call_Multivariate)

Heston_Call_Multivariate = Heston_UO_Call_PriceMultivariate(vector1, vector2, S10, S20, K, t, r, nbre_step, nbre_sim, rho, rho/1000)
print(Heston_Call_Multivariate)
