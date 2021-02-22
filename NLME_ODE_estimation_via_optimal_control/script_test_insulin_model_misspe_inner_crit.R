##################################################################
######## SUBJECT SPECIFIC PARAMETERS ESTIMATION FOR SECTION 4.2.2
##################################################################

pathnames <- list.files(pattern="[.]R$", path="function_model_and_sims//", full.names=TRUE);
sapply(pathnames, FUN=source);
pathnames <- list.files(pattern="[.]R$", path="function_util_estimation//", full.names=TRUE);
sapply(pathnames, FUN=source);

library('deSolve')
library('optimx')
library('snow')
library('numDeriv')

# test data retrieval
nb_subject =20
load(paste('data_gen_insulin_model_misspe_pop_',nb_subject,'.Rdata',sep=""))
load(paste('true_trajectory_insulin_model_misspe_pop_',nb_subject,'.Rdata',sep=""))

#########  Specification of the population parameter value (theta and delta)

coeff_mult_gi = list_Obs_subjects[[1]][[1]][[3]][[1]]
coeff_mult_x = list_Obs_subjects[[1]][[1]][[3]][[2]]

log_G0 = 5.52
log_I0 = 4.88
x0_mean <- c(coeff_mult_gi*exp(log_G0),coeff_mult_gi*exp(log_I0),coeff_mult_x *0)

log_Sg = -3.897740
log_p2  = -4.934301
log_Si  = -7.096229 
log_n = -1.815087 
log_gamma = -6.850169
log_h =  4.139365


theta_gen <- c(log_Sg,log_p2,log_Si,log_n,log_gamma,log_h)
true_theta  <- c(log_Sg,log_Si,log_n)

std = coeff_mult_gi*3

std_log_n =  0.2632241
std_param = c(std_log_n)
Psi = diag(std_param^2,1,1)
Chol_R= chol(Psi)
true_Delta_mat= chol((std^2)*solve(Psi))


#number of tested trial
nb_trial = 10


# Specification of the required model element and dimension for parameter estimation
dim_sub = nrow(Psi)
dim_control = 3
dim_syst = 3
dim_obs = 2

mat_A_pop  =function(t,vx,param_sub,param_pop,exo_par){insulin_model_matA(t,vx,param_sub,param_pop,exo_par)}
vect_r_pop  =function(t,param_sub,param_pop,exo_par){insulin_model_vectr(t,param_sub,param_pop,exo_par)}

mat_B = diag(dim_syst)
mat_C = rbind(c(1,0,0),c(0,1,0))

log_prior_function = function(param_pop,Delta_mat){return(0)}

# Specification of the mesh size
mesh_iter=40

#Weighing parameter specification
lambda_seq = c(1,10,1000,10^6)
nb_lambda = length(lambda_seq)

#For each tested trial and subject in it: 
#1) Compute the estimation of the b_i  where the optimization algorithm starts from
# 1.1) the true parameter value
# 1.2) a wrongly chosen parameter value to test practical identifiability issues
for (nb_t  in 1:nb_trial){
  
  true_trajectory_nb_t = list_true_trajectory[[nb_t]]
  
  Obs_subjects_nb_t = list_Obs_subjects[[nb_t]]
  State_ini_nb_t = list_State_ini[[nb_t]]
  bi_list_est= list_seq_true_bi[[nb_t]]
  bi_list_est_dc= list_seq_true_bi_dc[[nb_t]]
  x0i_list_est = list_seq_true_x0i[[nb_t]]
  known_x0i_list = list_seq_known_x0i[[nb_t]] 
  
  
  for (nbs in 1:nb_subject){
    
    Time_real_true_trajectory_i =  true_trajectory_nb_t[[nbs]][[1]]
    True_trajectory_i = true_trajectory_nb_t[[nbs]][[2]]
    
    Times_obs_i =Obs_subjects_nb_t[[nbs]][[1]]
    Y_i =Obs_subjects_nb_t[[nbs]][[2]]
    State_ini_i = State_ini_nb_t[[nbs]]
    true_bi =  bi_list_est[[nbs]]
    bi_dc = bi_list_est_dc[[nbs]]
    true_x0i = x0i_list_est[[nbs]]
    known_x0i =  known_x0i_list[[nbs]]
    
    nb_time_integ_i = mesh_iter*(length(Times_obs_i)-1)+1;
    
    Times_integ_i =  matrix(0, nrow = 1, ncol = nb_time_integ_i )
    weight_integ_i = matrix(0, nrow = 1, ncol = nb_time_integ_i )
    
    pseudo_Y_i =  matrix(0, nrow = dim_obs, ncol = nb_time_integ_i )
    pseudo_State_ini_i =  matrix(0, nrow = dim_syst, ncol = nb_time_integ_i )
    
    for (obs_i in 1:(length(Times_obs_i)-1)){
      
      delta_i = Times_obs_i[obs_i+1] - Times_obs_i[obs_i] 
      
      refined_interval_i =  seq(Times_obs_i[obs_i],Times_obs_i[obs_i+1],by = (delta_i/mesh_iter))
      refined_interval_i[length( refined_interval_i)] = Times_obs_i[obs_i+1]
      Times_integ_i[seq(mesh_iter*(obs_i-1)+1,(mesh_iter*obs_i+1),by =1)] = refined_interval_i
      
      weight_integ_i[mesh_iter*(obs_i-1)+1] =  1/(refined_interval_i[2]-Times_obs_i[obs_i])
      
      pseudo_Y_i[,mesh_iter*(obs_i-1)+1]= Y_i[,obs_i]
      pseudo_Y_i[,mesh_iter*obs_i+1] =Y_i[,obs_i+1]
      
      pseudo_State_ini_i[,seq(mesh_iter*(obs_i-1)+1,(mesh_iter*obs_i+1),by =1)]= matrix(rep(State_ini_i[,obs_i],mesh_iter+1),dim_syst,mesh_iter+1)
      
    }
    Times_integ_i = c( Times_integ_i)
    weight_integ_i = c(weight_integ_i)
    
    for (nb_l  in 1:nb_lambda){
      
      mat_U = lambda_seq[nb_l]*diag(dim_control)
      print(mat_U)
      
      list_inner_arg = list(Times_integ_i, pseudo_Y_i,pseudo_State_ini_i, weight_integ_i,true_bi,true_Delta_mat,mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C, known_x0i,2,0
                            ,true_theta,list(coeff_mult_gi,coeff_mult_x),type_algorithm_optim = 1,c(200,1))
      res_est_inner = est_param_oca_inner_criteria_par_version(list_inner_arg)
      
      
      list_inner_arg_dc = list(Times_integ_i, pseudo_Y_i,pseudo_State_ini_i, weight_integ_i,1.3*true_bi,true_Delta_mat,mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C, known_x0i,2,0
                               ,true_theta,list(coeff_mult_gi,coeff_mult_x),type_algorithm_optim = 1,c(200,1))
      res_est_inner_dc = est_param_oca_inner_criteria_par_version(list_inner_arg_dc)
      
      
      
      print(true_bi)
      print("estimation from true initial guess ")
      print(  res_est_inner$bi_est)
      print("estimation from wrong initial guess ")
      print(  res_est_inner_dc$bi_est)
      
      main_title = paste("Solution for subject ",nbs,", lambda = ", lambda_seq[nb_l],sep="")
      matplot(Time_real_true_trajectory_i,t(True_trajectory_i),type="l" ,xlab = "Time", ylab = "X",main = main_title,col="red",ylim= c(0,1.2*max(True_trajectory_i)))
      lines(Times_integ_i,res_est_inner$Xi_opt[1,],col="green")
      lines(Times_integ_i,res_est_inner$Xi_opt[2,],col="green")
      
      points(Times_obs_i , Y_i[1,],pch=16,bg="black")
      points(Times_obs_i , Y_i[2,],pch=16,bg="black")
      legend("topright",c("True Trajectories","Est. Optimal Trajectories","Observations"),pch=c(NA,NA,16),lty=c(1,1,NA),bg=c(NA,NA,"black"),col=c("red", "green", "black"), lwd = c(1,1,1))
    }
    
    
  }
  
}
