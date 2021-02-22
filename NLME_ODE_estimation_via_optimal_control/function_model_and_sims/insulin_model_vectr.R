insulin_model_vectr<- function(t,param_sub,param_pop,exo_par)
{
  coeff_gi = exo_par[[1]]
  coeff_x = exo_par[[2]]
  
  #Sg = exp(-3.897740)
  p2  = exp(-4.934301)
  #Si  = exp(-7.096229) 
  #n = exp(-1.815087) 
  gamma = exp(-6.850169)
  h =  exp(4.139365)
  Gb =100
  Ib = 10
  
  Sg = exp(param_pop[1])
  Si = exp(param_pop[2])
  n = exp(param_pop[3]+param_sub[1])
  
  
  
  #dgdt = coeff_gi*Sg*Gb - Sg*g -(x/coeff_x)*g
  #didt =    -n*i+coeff_gi*n*Ib+gamma*Time*g -coeff_gi*gamma*Time*h
  #dxdt =  -p2*x +coeff_x*p2*Si*(i/coeff_gi)-coeff_x*p2*Si*Ib
  
  res_t = matrix(0,3,1)
  
  res_t[1,1] = coeff_gi*Sg*Gb
  res_t[2,1] = coeff_gi*n*Ib-coeff_gi*gamma*t*h
  res_t[3,1] = -coeff_x*p2*Si*Ib
  
  return(res_t)
  
}