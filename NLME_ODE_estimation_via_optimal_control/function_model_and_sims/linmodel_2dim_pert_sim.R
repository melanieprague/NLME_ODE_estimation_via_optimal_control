linmodel_2dim_pert_sim <- function(deb,h,fin,Pars,y0,sd_stoch_pert){
  
  Times_obs = seq(deb,fin, by=h) 
  nb_obs = length(Times_obs)
  
  dim_syst =  2
  
  
  Res_sim = matrix(rep(0,dim_syst*nb_obs),dim_syst,nb_obs)
  Res_sim[,1] = y0
  for (i in 1 : (nb_obs-1)){
    delta_i = Times_obs[i+1]-Times_obs[i]
    State = Res_sim[,i]
    theta1 = exp(Pars[1]+Pars[3])
    theta2= exp(Pars[2])
    
    X1 = State[1]
    X2 = State[2]
    
    dX1 =  theta2*X2 -theta1*X1
    dX2 = -theta2*X2
    
    part_det_i = matrix(c(0,0),dim_syst,1)
    part_det_i[1] =  dX1
    part_det_i[2] =  dX2
    
    Res_sim[1,i+1] = State[1] +delta_i*part_det_i[1]
    Res_sim[2,i+1] = State[2] +delta_i*part_det_i[2]+ sqrt(delta_i)*sd_stoch_pert*rnorm(1,0)
  }
  return(list(Times_obs,Res_sim))
}