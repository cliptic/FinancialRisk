WD = setwd('D:/GitHub/Financial_Risk')
source('2_BSMC_barrier_functions.R')
source('2_Heston_barrier_functions.R')

### 2 ####

nbre_step = 500
nbre_sim = 20000

r = 0.01
S0 = 175
t = 30/365
K = 177
B = 190
sig = 0.2412037 # from previous calibration
vector = c(0.0467102,  0.1445764, -0.2743338,  0.9507367,  0.5324521)

#price for up and out call option and survival probability:
print(BS_MC_UO_Call_Price(S0,K,B,T,sig,r, nbre_step, nbre_sim)[1])
print(Heston_upout_barrier_call(vector,S0,K,B,t,r, nbre_step, nbre_sim))

#delta, gamma, vega for up and out call option:

print(BS_MC_UO_Call_Delta(S0,K,B,T,sig,r, nbre_step, nbre_sim, S0/1000)[5])
print(Heston_upout_delta(vector,S0,K,B,t,r,S0/1000, nbre_step, nbre_sim))

print(BS_MC_UO_Call_Gamma(S0,K,B,T,sig,r, nbre_step, nbre_sim, S0/1000)[5])
print(Heston_upout_gamma(vector,S0,K,B,t,r,S0/1000, nbre_step, nbre_sim))

print(BS_MC_UO_Call_Vega(S0,K,B,T,sig,r, nbre_step, nbre_sim, sig/1000)[5])
print(Heston_upout_vega(vector,S0,K,B,t,r,sig/1000, nbre_step, nbre_sim))

print(BS_MC_UO_Call_Price(S0,K,B,T,sig,r, nbre_step, nbre_sim)[2])
Heston_upout_survival(vector,S0,K,B,t,r, nbre_step, nbre_sim)

## plots

strike_range = c(150:200)
delta = lapply(x, function(u) BS_upout_delta(sig,u,K,B,t,r,S0/1000))
plot(strike_range, delta, type = "l", col = "lightblue4",  ylab = "Delta",xlab = "Strike prices", main = "Delta of barrier  Call")
gamma = lapply(x, function(u) BS_upout_gamma(sig,u,K,B,t,r,S0/1000))
plot(strike_range, gamma, type = "l", col = "lightblue4",  ylab = "Gamma",xlab = "Strike prices", main = "Gamma of barrier  Call")
