linmodel_2dim_sim <- function(Time, State, Pars)
{
  theta1 = exp(Pars[1]+Pars[3])
  theta2= exp(Pars[2])
  
  X1 = State[1]
  X2 = State[2]
  
  dX1 =  theta2*X2 -theta1*X1
  dX2 = -theta2*X2
  
  return(list(c(dX1, dX2)))
}
