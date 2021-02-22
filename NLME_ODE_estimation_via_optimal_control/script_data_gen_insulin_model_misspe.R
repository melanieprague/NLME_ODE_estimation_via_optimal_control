############################################
######## DATA GENERATION FOR SECTION 4.2.2
############################################

pathnames <- list.files(pattern="[.]R$", path="function_model_and_sims//", full.names=TRUE);
sapply(pathnames, FUN=source);
pathnames <- list.files(pattern="[.]R$", path="function_util_estimation//", full.names=TRUE);
sapply(pathnames, FUN=source);

library('deSolve')
library('optimx')
library('snow')
library('numDeriv')

# specification of population size, observation numbers & scaling parameter for state value & Times obs
nb_subject =20
nb_trial = 100

Times_ref <- seq(0, 180, by = 180/200)
Times_obs <- seq(0,180, by = 180/5)
coeff_mult_gi = 0.001 
coeff_mult_x = 10 

nb_obs <- length(Times_obs)

# data file name
rep_file = paste('data_gen_insulin_model_misspe_pop_',nb_subject,'.Rdata',sep="")
rep_file_trajectory = paste('true_trajectory_insulin_model_misspe_pop_',nb_subject,'.Rdata',sep="")

# specification of population parameters 
log_G0 = 5.52
log_I0 = 4.88
x0_mean <- c(coeff_mult_gi*exp(log_G0),coeff_mult_gi*exp(log_I0),coeff_mult_x *0)

log_Sg = -3.897740
log_p2  = -4.934301
log_Si  = -7.096229 
log_n = -1.815087 
log_gamma = -6.850169
log_h =  4.139365

std = coeff_mult_gi*3

std_log_n =  0.2632241
std_param = c(std_log_n)
dim_sub = length(std_param)

Psi = diag(std_param^2,1,1)
Chol_R= chol(Psi)
Delta_mat= chol((std^2)*solve(Psi))

theta_gen <- c(log_Sg,log_p2,log_Si,log_n,log_gamma,log_h)


list_Obs_subjects = list()
list_true_trajectory = list()
list_State_ini= list()
list_seq_true_bi = list()
list_seq_true_bi_dc = list()
list_seq_true_x0i = list()
list_seq_known_x0i = list()

for (nbt in 1:nb_trial){
  
  
  Obs_subjects = list()
  true_trajectory = list()
  State_ini= list()
  seq_true_bi = list()
  seq_true_bi_dc = list()
  seq_true_x0i = list()
  seq_known_x0i = list()
  
  for (nbs in 1:nb_subject){
    
    true_bi = Chol_R%*%rnorm(dim_sub ,0);
    x0i =c(  coeff_mult_gi*exp(log_G0 +0.17*rnorm(1,0)), coeff_mult_gi*exp(log_I0+0.1*rnorm(1,0)),coeff_mult_x*max(0.001+0.0001*rnorm(1,0),0))
    
    out_ref <- ode( y = x0i, times = Times_ref,func =insulin_model_sim, parms =c(theta_gen, true_bi,coeff_mult_gi,coeff_mult_x))
    out_pert_ref = insulin_model_pert_sim(deb=0,h=Times_obs[length(Times_obs)]/1000,fin=Times_obs[length(Times_obs)],
                                          Pars = c(theta_gen,true_bi,coeff_mult_gi,coeff_mult_x),x0i ,sd_stoch_pert = c(coeff_mult_gi*2,coeff_mult_gi*2,coeff_mult_x*0.0002))
    
    out_pert = t(out_pert_ref[[2]][,seq(1,length(out_pert_ref[[1]]),1001/5)])
    out_pert  =rbind( out_pert,  out_pert_ref[[2]][,1001])
    Time_pert  =  t(out_pert_ref[[1]][seq(1,length(out_pert_ref[[1]]),1001/5)])
    Time_pert = cbind( Time_pert,out_pert_ref[[1]][1001])
    Y_noise = t(out_pert[,1:2] +std*cbind(rnorm(n = nb_obs),rnorm(n=nb_obs)))
    
    if (nbt ==1){
      if (nbs==1){
        matplot(Times_ref ,out_ref[,2], type = "l", xlab = "Time", ylab = "",main = "Perturbed Insulin ODE model",ylim = c(0,1.4*max(out_ref[,2:3])))
        lines(Times_ref,out_ref[,3],lty=1) 
        lines(Times_ref,out_ref[,4],lty=1) 
        lines(out_pert_ref[[1]],out_pert_ref[[2]][1,],lty=2) 
        lines(out_pert_ref[[1]],out_pert_ref[[2]][2,],lty=2) 
        lines(out_pert_ref[[1]],out_pert_ref[[2]][3,],lty=2) 
        points(Times_obs,Y_noise[1,],pch=16,bg="black") 
        points(Times_obs,Y_noise[2,],pch=16,bg="black") 
        legend("topright",c("True stochastic curves","Deterministic approximation (ODE solution)","Observations"),pch=c(NA,NA,16),lty=c(2,1,NA),bg=c(NA,NA,"black"),col=c("black", "black", "black"), lwd = c(1,1,1))
      }
    }
    
    true_trajectory =  append(true_trajectory,list(out_pert_ref))
    
    Obs_subjects = append(Obs_subjects,list(list(Time_pert,Y_noise,list(coeff_mult_gi,coeff_mult_x))))
    State_ini = append( State_ini,list(rbind(Y_noise,matrix(0,1,nb_obs))))
    seq_true_bi  = append(seq_true_bi, list(true_bi))
    seq_true_bi_dc = append(seq_true_bi_dc, list(1.3*true_bi))
    seq_true_x0i  = append(seq_true_x0i,list(x0i))
    seq_known_x0i  = append(seq_known_x0i,list(list()))
  }
  
  
  list_true_trajectory = append(list_true_trajectory,list(true_trajectory))
  
  list_Obs_subjects  = append( list_Obs_subjects,list( Obs_subjects ))
  list_State_ini = append(list_State_ini,list(State_ini))
  list_seq_true_bi = append(list_seq_true_bi,list( seq_true_bi))
  list_seq_true_bi_dc  = append(list_seq_true_bi_dc ,list(seq_true_bi_dc))
  list_seq_true_x0i = append(list_seq_true_x0i, list(seq_true_x0i) )
  list_seq_known_x0i = append(list_seq_known_x0i, list(seq_known_x0i) )
  
}
save(file=rep_file,list_Obs_subjects,list_State_ini,list_seq_true_bi , list_seq_true_bi_dc, list_seq_true_x0i, list_seq_known_x0i )
save(file=rep_file_trajectory,list_true_trajectory)