Lin_2dim_matA <- function(t,param_sub,param_pop){
  

  th1 = exp(param_sub[1] + param_pop[1])
  th2 =  exp(param_pop[2])
  
  res_t = matrix(0,2,2)
  
  # dX1 =  theta2*X2 -theta1*X1
  res_t[1,1] = -th1
  res_t[1,2] =  th2
  
  # dX2 =-theta2*X2
  res_t[2,2] = -th2
 
  return(res_t)
}