Ebola_vaccine_model_sim <- function(Time, State, Pars)
{
  
  delta_S = exp(Pars[1])
  delta_L = exp(Pars[2])
  phi_S = exp(Pars[3] +Pars[6] )
  phi_L = exp(Pars[4] +Pars[7] )
  delta_Ab = exp(Pars[5]  +Pars[8])
  
  Ab = State[1]
  
  S_cell_part =  phi_S*exp(-delta_S*Time)
  L_cell_part =  phi_L*exp(-delta_L*Time)
  
  dAb = S_cell_part+L_cell_part -delta_Ab*Ab 
  
  return(list(c(dAb)))
  
}