insulin_model_sim <- function(Time, State, Pars)
{
  #Fixed: Sg + p2 + Si + n + gamma + h + G0 + I0 ~ 1 
  #Sg        p2        Si         n     gamma         h        G0        I0 
  #-3.897740 -4.934301 -7.096229 -1.815087 -6.850169  4.139365
  coeff_gi = Pars[8]
  coeff_x = Pars[9]
  
  Gb =100
  Ib = 10
  Sg = exp(Pars[1])
  p2  = exp(Pars[2]) 
  Si  = exp(Pars[3]) 
  n = exp(Pars[4]+Pars[7]) 
  gamma = exp(Pars[5]) 
  h = exp(Pars[6]) 
  
  g = State[1]
  i = State[2]
  x = State[3]
  
  
  dgdt = coeff_gi*Sg*Gb - Sg*g -(x/coeff_x)*g
  didt =    -n*i+coeff_gi*n*Ib+gamma*Time*g -coeff_gi*gamma*Time*h
  dxdt =  -p2*x +coeff_x*p2*Si*(i/coeff_gi)-coeff_x*p2*Si*Ib
  
  return(list(c(dgdt, didt, dxdt)))
}