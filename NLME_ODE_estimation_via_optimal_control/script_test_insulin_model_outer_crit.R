###########################################################
######## POPULATION PARAMETERs ESTIMATION FOR SECTION 4.2.1
###########################################################

pathnames <- list.files(pattern="[.]R$", path="function_model_and_sims//", full.names=TRUE);
sapply(pathnames, FUN=source);
pathnames <- list.files(pattern="[.]R$", path="function_util_estimation//", full.names=TRUE);
sapply(pathnames, FUN=source);

library('deSolve')
library('optimx')
library('snow')
library('numDeriv')



#########  Specification of the population parameter value (theta and delta) as first guess for the outer criteria optimization algorithm
# test data retrieval
nb_subject =20
load(paste('data_gen_insulin_model_pop_',nb_subject,'.Rdata',sep=""))


####### Specify if theta_Si has to be estimated (=1 yes, no otherwise)
Si_estimation = 1

### scaling parameter specification
coeff_mult_gi = list_Obs_subjects[[1]][[1]][[3]][[1]]
coeff_mult_x = list_Obs_subjects[[1]][[1]][[3]][[2]]

#########  Specification of the population parameter value (theta and delta)
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

if (Si_estimation==1){
  true_theta  <- c(log_Sg,log_Si,log_n)
}else{
  true_theta  <- c(log_Sg,log_n)
}

std = coeff_mult_gi*3

std_log_n =  0.2632241
std_param = c(std_log_n)
Psi = diag(std_param^2,1,1)
Chol_R= chol(Psi)
true_Delta_mat= chol((std^2)*solve(Psi))

#number of tested trial
nb_trial = 100


# Specification of the required model element and dimension for parameter estimation
dim_sub = nrow(Psi)
dim_control = 3
dim_syst = 3
dim_obs = 2


if (Si_estimation==1){
  mat_A_pop  =function(t,vx,param_sub,param_pop,exo_par){insulin_model_matA(t,vx,param_sub,param_pop,exo_par)}
  vect_r_pop  =function(t,param_sub,param_pop,exo_par){insulin_model_vectr(t,param_sub,param_pop,exo_par)}
}else{
  mat_A_pop  =function(t,vx,param_sub,param_pop,exo_par){insulin_model_matA_Siknown(t,vx,param_sub,param_pop,exo_par)}
  vect_r_pop  =function(t,param_sub,param_pop,exo_par){insulin_model_vectr_Siknown(t,param_sub,param_pop,exo_par)}
}


mat_B = diag(dim_syst)
mat_C = rbind(c(1,0,0),c(0,1,0))

log_prior_function = function(param_pop,Delta_mat){return(0)}

# precize the if the subject specific parameters have to be estimated or are fixed
est_pop_parameter_only=0

# precize the if the subject specific parameters variance has to be estimated 
delta_known =0

# precize if Asymptotic variance-covariance has to be estimated
est_var = 1

# Specification of the mesh size
mesh_iter=30

#Weighing parameter specification
lambda_seq = c(10^6,10^7,10^8)
nb_lambda = length(lambda_seq)


#For each tested trial and subject in it: 
#1) Compute the estimation of the population paramters where the optimization algorithm starts from
# 1.1) the true parameter value (registered in list_res_trial)
# 1.2) a wrongly chosen parameter value to test practical identifiability issues (registered in list_res_trial_dc)

list_res_trial =  list()
list_res_trial_dc =  list()

# specify if the inner criteria optimization is parallelized (using snow package) among the subjects
cl_cur <- makeSOCKcluster(rep("localhost",nb_subject))
#cl_cur <- list()

for (nb_t in 1:nb_trial){
  
  list_res_est= list()
  list_res_est_dc= list()
  
  
  Obs_subjects_nb_t = list_Obs_subjects[[nb_t]]
  State_ini_nb_t = list_State_ini[[nb_t]]
  bi_list_est= list_seq_true_bi[[nb_t]]
  bi_list_est_dc= list_seq_true_bi_dc[[nb_t]]
  x0i_list_est = list_seq_true_x0i[[nb_t]]
  known_x0i_list = list_seq_known_x0i[[nb_t]] 
  
  theta_pop_ini = true_theta
  delta_ini = log(diag(true_Delta_mat))
  
  theta_pop_ini_dc = 1.2*theta_pop_ini
  delta_ini_dc =  1.2*delta_ini 
  
  
  for (lambda  in lambda_seq){
    T1<-Sys.time() 
    
    mat_U = lambda*diag(dim_control)
    out_outer_crit = est_param_oca_outer_criteria_prof_par_version(Obs_subjects_nb_t,State_ini_nb_t,param_pop_ini=theta_pop_ini,log_Delta_vect_ini=delta_ini,
                                                                   mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C,
                                                                   mesh_iter,bi_list_est,known_x0i_list,
                                                                   est_pop_parameter_only,delta_known,log_prior=log_prior_function,
                                                                   type_optim =0,est_var =est_var,type_optim_inner=1,nb_iter_max=c(1500,200,0),cl_cur)
    
    
    out_outer_crit_dc = est_param_oca_outer_criteria_prof_par_version(Obs_subjects_nb_t,State_ini_nb_t,param_pop_ini=theta_pop_ini_dc,log_Delta_vect_ini=delta_ini_dc,
                                                                              mat_U,mat_A_pop,vect_r_pop,mat_B,mat_C,
                                                                              mesh_iter,bi_list_est,known_x0i_list,
                                                                              est_pop_parameter_only,delta_known,log_prior=log_prior_function,
                                                                              type_optim =0,est_var =est_var,type_optim_inner=1,nb_iter_max=c(1500,200,0),cl_cur)
    
    
        T2<-Sys.time() 
    time_parlapply = T2 -T1
    
    print("True population parameter values")
    print(c(true_theta,log(diag(true_Delta_mat))))
    
    print("parameter  estimation from true initial guess ")
    print(c(out_outer_crit$theta_est,out_outer_crit$delta_est))
    
    print("parameter  estimation from wrong initial guess ")
    print(c(out_outer_crit_dc$theta_est,out_outer_crit_dc$delta_est))
    
    if (est_var ==1){
      print("Variance estimation from true initial guess")
      print(diag(out_outer_crit$Variance_component$est_var_cov_Matrix))
      
      print("Variance estimation from wrong initial guess")
      print(diag(out_outer_crit_dc$Variance_component$est_var_cov_Matrix))
    }
    
    theta_pop_ini = out_outer_crit$theta_est
    delta_ini = out_outer_crit$delta_est
    
    theta_pop_ini_dc = out_outer_crit_dc$theta_est
    delta_ini_dc = out_outer_crit_dc$delta_est
    
    list_res_est = append(list_res_est,list(list(out_outer_crit,time_parlapply)))
    list_res_est_dc=  append(list_res_est_dc,list(list(out_outer_crit_dc,time_parlapply)))
    
  }
  list_res_trial = append(list_res_trial,list(list_res_est))
  list_res_trial_dc = append(list_res_trial_dc,list(list_res_est_dc))
}
stopCluster(cl_cur)