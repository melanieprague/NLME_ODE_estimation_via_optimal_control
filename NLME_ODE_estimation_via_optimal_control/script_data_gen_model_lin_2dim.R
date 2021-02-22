############################################
######## DATA GENERATION FOR SECTION 4.1.1
############################################

pathnames <- list.files(pattern="[.]R$", path="function_model_and_sims//", full.names=TRUE);
sapply(pathnames, FUN=source);
pathnames <- list.files(pattern="[.]R$", path="function_util_estimation//", full.names=TRUE);
sapply(pathnames, FUN=source);

library('deSolve')
library('optimx')
library('snow')
library('numDeriv')

# specification of population size and observation numbers
nb_subject =50
nb_trial = 100
Times_obs <- seq(0,10, by =10/10)
Times_ref<- seq(0,10, by =10/200)


# data file name
rep_file = paste('data_gen_linmodel_2dim_pop',nb_subject,'.Rdata',sep="")

# specification of population parameters 
x0_mean <- c(2,3)

log_theta1 = log(0.5)
log_theta2 = log(0.2)

theta_gen <- c(log_theta1,log_theta2)

nb_obs <- length(Times_obs)

std = 0.05

std_log_theta1  = 0.5

std_param = c(std_log_theta1)

Psi = diag(std_param^2,1,1)
Chol_R= chol(Psi)
dim_sub = nrow(Psi)
Delta_mat= chol((std^2)*solve(Psi))


list_Obs_subjects = list()
list_seq_true_bi = list()
list_seq_true_bi_dc = list()
list_seq_true_x0i = list()
list_seq_known_x0i = list()

# Data generation
for (nbt in 1:nb_trial){
  
  Obs_subjects = list()
  seq_true_bi = list()
  seq_true_bi_dc = list()
  seq_true_x0i = list()
  seq_known_x0i = list()
  
  for (nbs in 1:nb_subject){
    true_bi = Chol_R%*%rnorm(1,0)
    x0i = pmax(x0_mean*(1+0.5*rnorm(1,0)),0.5)
    #   x0i = pmax(x0_mean,0.5)
    
    out <- ode(func = linmodel_2dim_sim, y = x0i, parms =c(theta_gen,true_bi), times = Times_obs)
    out_ref <- ode(func = linmodel_2dim_sim, y = x0i, parms =c(theta_gen,true_bi), times = Times_ref)
    Times_obs <- out[,"time"]
    Y_noise = matrix(out[,2] +std*rnorm(n = nb_obs),1,nb_obs)
    
    Obs_subjects = append(Obs_subjects,list(list(Times_obs,Y_noise,list())))
    seq_true_bi  = append(seq_true_bi, list(true_bi))
    seq_true_bi_dc = append(seq_true_bi_dc, list(1.3*true_bi))
    seq_true_x0i  = append(seq_true_x0i,list(x0i))
    seq_known_x0i  = append(seq_known_x0i,list(list()))
    if (nbt ==1){
    if (nbs==1){
      matplot(Times_ref ,out_ref[,2], type = "l", xlab = "Time", ylab = "",main = "2 dimensional linear ODE",ylim = c(0,8))
      lines(Times_ref,out_ref[,3],lty=2) 
      points(Times_obs,Y_noise[1,],pch=16,bg="black") 
    }else{
      lines(Times_ref,out_ref[,2],lty=1) 
      lines(Times_ref,out_ref[,3],lty=2) 
      points(Times_obs,Y_noise[1,],pch=16,bg="black") 
    }
    legend("topright",c("X1","X2","Observations"),pch=c(NA,NA,16),lty=c(1,2,NA),bg=c(NA,NA,"black"),col=c("black", "black", "black"), lwd = c(1,1,1))
    
    }
  }
  
  list_Obs_subjects  = append( list_Obs_subjects,list( Obs_subjects ))
  list_seq_true_bi = append(list_seq_true_bi,list( seq_true_bi))
  list_seq_true_bi_dc  = append(list_seq_true_bi_dc ,list(seq_true_bi_dc))
  list_seq_true_x0i = append(list_seq_true_x0i, list(seq_true_x0i) )
  list_seq_known_x0i = append(list_seq_known_x0i, list(seq_known_x0i) )
}

save(file=rep_file,list_Obs_subjects,list_seq_true_bi , list_seq_true_bi_dc, list_seq_true_x0i, list_seq_known_x0i )