WD = setwd('D:/GitHub/Financial_Risk')
source('data2018.R')
source('callHestoncf.R')

library(NMOF) 

#### ex1

#### Calibration  

#### first model - BS ####

#number of maturities from list
ma = 5

BS_Call_Price <- function(S0,K,t,sig)
{
  r = 0.01
  d1 = (log(S0/K) + (r+sig^2/2)*t)/(sig * sqrt(t))
  d2 = d1 - sig * sqrt(t)
  Nd1 <- pnorm(d1) 
  Nd2 <- pnorm(d2)
  DF = exp(-r*t)
  BSCall <- S0 * Nd1 - K * DF * Nd2
  
  return(BSCall)
  
}

Vega <- function(S0,K,T,sig)
{
  r=0.01
  d1 = (log(S0/K) + (r+sig^2/2)*T)/(sig * sqrt(T))
  nd1 <- dnorm(d1)   #dnorm(x) returns the Standard Gaussian Density of x
  Vega <- S0 * sqrt(T) * nd1  
}

AAPL_error<-function(sig)
{
  spot = Spot_AAPL
  error = 0
  for (i in 1:10)
  {
    for (j in 1:ma)
    {
      
      strike = Strikes_AAPL[i]
      t = Call_AAPL[[j]]$maturity
      market_price = Call_AAPL[[j]]$Average_Price[i]
      error = error + sqrt((market_price - BS_Call_Price(spot,strike,t,sig))^2)
    }
  }
  return(error)
}

FB_error<-function(sig)
{
  spot = Spot_FB
  error = 0
  for (i in 1:10)
  {
    for (j in 1:ma)
    {
      
      strike = Strikes_FB[i]
      t = Call_FB[[j]]$maturity
      market_price = Call_FB[[j]]$Average_Price[i]
      error= error + sqrt((market_price - BS_Call_Price(spot,strike,t,sig))^2)
    }
  }
  return(error)
}

sig_AAPL<-optim(par=0.1,fn=AAPL_error,method = 'L-BFGS-B',lower = 0.01,upper = 1)$par
print(sig_AAPL)

sig_FB<-optim(par=0.1,fn=FB_error,method = 'L-BFGS-B',lower = 0.01,upper = 1)$par
print(sig_FB)

#second model - Heston ####

Heston_Call_Price<-function(S,K,t,vector) # vector is a vector of 5 parameters for heston model
{
  r = 0.01
  q = 0
  v0 = vector[1]
  vT = vector[2]
  rho = vector[3]
  k = vector[4]
  sigma = vector[5]
  return(callHestoncf(S,K,t,r,q,v0,vT,rho,k,sigma))
}

AAPL_error<-function(vector)
{
  spot = Spot_AAPL
  error = 0
  for (i in 1:10)
  {
    for (j in 1:ma)
    {
      strike = Strikes_AAPL[i]
      t = Call_AAPL[[j]]$maturity
      market_price = Call_AAPL[[j]]$Average_Price[i]
      error= error + sqrt((market_price - Heston_Call_Price(spot,strike,t,vector))^2)
    }
  }
  return(error)
}

FB_error<-function(vector)
{
  spot = Spot_FB
  error = 0
  for (i in 1:10)
  {
    for (j in 1:ma)
    {
      strike = Strikes_FB[i]
      t = Call_FB[[j]]$maturity
      market_price = Call_FB[[j]]$Average_Price[i]
      error= error + sqrt((market_price - Heston_Call_Price(spot,strike,t,vector))^2)
    }
  }
  return(error)
}


vector_AAPL<-optim(par=c(0.2,0.2,-0.3,1,0.5),fn=AAPL_error,method = 'L-BFGS-B',
                   lower = c(0.01,0.01,-0.99,0.5,0.01),
                   upper = c(0.99,0.99,-0.1,5,1))$par

print(vector_AAPL) # result: c(0.0467102,  0.1445764, -0.2743338,  0.9507367,  0.5324521)

vector_FB<-optim(par=c(0.2,0.2,-0.3,1,0.5),fn=FB_error,method = 'L-BFGS-B',
                   lower = c(0.01,0.01,-0.99,0.5,0.01),
                   upper = c(0.99,0.99,-0.1,5,1))$par

print(vector_FB)  # result: c(0.04864663,  0.23483467, -0.41129397,  0.92265515,  0.59116282)

# plot results ####

#AAPL
plot_calls_AAPL<-function(maturity)
{
  t = Call_AAPL[[maturity]]$maturity
  y_ask = Call_AAPL[[maturity]]$Ask
  y_bid = Call_AAPL[[maturity]]$Bid
  x = Strikes_AAPL
  y_bs = lapply(x, function(x) BS_Call_Price(Spot_AAPL,x,t,sig_AAPL))
  y_heston = lapply(x, function(x) Heston_Call_Price(Spot_AAPL,x,t,vector_AAPL))
  dev.new() 
  plot(x,y_bs, type = "b", pch = 16, lty = 2, col = "coral2", ylab = "Call prices",xlab = "Strike prices")
  lines(x, y_heston, type = "b", pch = 16, lty = 2, col = "cyan3")
  lines(x, y_ask)
  lines(x,y_bid)
  legend('topright', c("Black-Scholes","Heston"), bty = "n", col= c("coral2","cyan3"), border = c("coral2","cyan3"), fill = c("coral2","cyan3"), main = "AAPL option prices")
}

plot_calls_AAPL(4)

#FB

plot_calls_FB<-function(maturity)
{
  t = Call_FB[[maturity]]$maturity
  y_ask = Call_FB[[maturity]]$Ask
  y_bid = Call_FB[[maturity]]$Bid
  x = Strikes_FB
  y_bs = lapply(x, function(x) BS_Call_Price(Spot_FB,x,t,sig_FB))
  y_heston = lapply(x, function(x) Heston_Call_Price(Spot_FB,x,t,vector_FB))
  dev.new()
  plot(x,y_bs, type = "b",pch = 16, lty = 2, col = "coral2", ylab = "Call prices",xlab = "Strike prices", main = "FB option prices")
  lines(x, y_heston, type = "b", pch = 16, lty = 2, col = "cyan3")
  lines(x, y_ask)
  lines(x, y_bid)
  legend('topright', c("Black-Scholes","Heston"), bty = "n", col= c("coral2","cyan3"), border = c("coral2","cyan3"), fill = c("coral2","cyan3"))
}

plot_calls_FB(1)

#implied volatility surfaces

Implied_volatility <- function(S0, K, t, Cmkt)
{ r = 0.01
  max_iter = 10000
precision = 0.00001
sig = 0.3  
for (i in 1:max_iter) 
{
  Call = BS_Call_Price(S0,K,t,sig)
  vegaBS = Vega(S0,K,t,sig) 
  error = abs(Call - Cmkt)
  
  if(error > precision){
      sig = sig - (Call - Cmkt)/vegaBS
      return(sig)}
  }
}

plot_impvol_AAPL<-function(maturity)
{
  t = Call_AAPL[[maturity]]$maturity
  
  x = Strikes_AAPL
  imp_ask = c()
  imp_bid = c()
  imp_bs = c()
  imp_heston = c()
  for (strike in 1:10)
  {
    imp_ask = c(imp_ask,Implied_volatility(Spot_AAPL,x[strike],t,
                                           Call_AAPL[[maturity]]$Ask[strike]))
    
    imp_bid = c(imp_bid,Implied_volatility(Spot_AAPL,x[strike],t,
                                           Call_AAPL[[maturity]]$Bid[strike]))
    
    imp_bs = c(imp_bs,Implied_volatility(Spot_AAPL,x[strike],t,
                                         BS_Call_Price(Spot_AAPL,x[strike],t,sig_AAPL)))
    
    imp_heston = c(imp_heston,Implied_volatility(Spot_AAPL,x[strike],t,
                                                 Heston_Call_Price(Spot_AAPL,x[strike],t,vector_AAPL)))
  }

 # dev.new() 
  plot(x,imp_bid, type="l", col = 1, ylim=c(0.23,0.27), ylab = "implied volatilities",xlab = "Strikes prices", main = "AAPL implied volatilities for Call")
  lines(x, imp_heston, type = "b",pch = 16, lty = 2, col = "cyan3")
  lines(x, imp_bs,type = "b",pch = 16, lty = 2, col = "coral2")
  lines(x,imp_ask, col = 1, type = "l")
  legend('topright', c("Black-Scholes","Heston"), bty = "n", col= c("coral2","cyan3"), border = c("coral2","cyan3"), fill = c("coral2","cyan3"))
}

plot_impvol_AAPL(4)

plot_impvol_FB<-function(maturity)
{
  t = Call_FB[[maturity]]$maturity
  
  x = Strikes_FB
  imp_ask = c()
  imp_bid = c()
  imp_bs = c()
  imp_heston = c()
  for (strike in 1:10)
  {
    imp_ask = c(imp_ask,Implied_volatility(Spot_FB,x[strike],t,
                                           Call_FB[[maturity]]$Ask[strike]))
    
    imp_bid = c(imp_bid,Implied_volatility(Spot_FB,x[strike],t,
                                           Call_FB[[maturity]]$Bid[strike]))
    
    imp_bs = c(imp_bs,Implied_volatility(Spot_FB,x[strike],t,
                                         BS_Call_Price(Spot_FB,x[strike],t,sig_FB)))
    
    imp_heston = c(imp_heston,Implied_volatility(Spot_FB,x[strike],t,
                                                 Heston_Call_Price(Spot_FB,x[strike],t,vector_FB)))
  }
  
  # dev.new()   
  plot(x,imp_bid, type = "l", col = 1, ylim=c(0.26,0.31), ylab = "Implied volatilities",xlab = "Strike prices", main = "FB implied volatilities for Call")
  lines(x, imp_heston, type = "b", pch = 16, lty = 2, col = "cyan3")
  lines(x, imp_bs, type = "b", pch = 16, lty = 2, col = "coral2")
  lines(x,imp_ask, col = 1, type = "l")
  legend('topright', c("Black-Scholes","Heston"), bty = "n", col= c("coral2","cyan3"), border = c("coral2","cyan3"), fill = c("coral2","cyan3"))
}

plot_impvol_FB(5)
