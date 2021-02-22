The code is primarly made to reproduce the examples given in section 4 "Results on simulated data" but it is generic enough to be easily adapted to
other NLME-ODE models and data set.

For the code to work you simply need to set this folder as the working directory, it is composed of:

Folders:

1) "function_model_and_sims" containing the functions specifying:
			     - the models used to generate data set, 
                             - the pseudo-linear representations of ODE models, 
                             - log-prior distribution for model (4.6).
2) "function_util_estimation" containing the generic functions used to estimate:
			     - the population parameters,
			     - their variance-covariance matrix asymptotic estimators,  
                             - the subject-specific parameters. 
All of these functions are commented to be used for new models.

Scripts:

For each subsection presented in section 4, 3 scripts are given:

I)   A first script to generate and register data which needs to be executed before estimation. 
     In addition, for perturbed models the true stochastic trajectories are registered for the sake of comparison.

II)  A second script to estimate the true subject specific parameters b_i_star according to equation (2.2) as well as the related optimal trajectories. 
     The estimation is done by using the true value b_i_star as initial guess for the optimization algorithm and a wrong value for checking identifiability issues. 

III) A third script to estimate the  population parameters (theta_star,delta_star) according to equation (2.5). 
     The parameter delta_star corresponds to the parametrization  for the variance  Psi_star(delta_star) given at the end of section 2.3. 
     The estimation is done by using the true value (theta_star,delta_star) as initial guess for the optimization algorithm and a wrong value for checking identifiability issues. 

The precise name of script for each subsections is given below:

For subsection 4.1.1:
- script_data_gen_model_lin_2dim.R:        generate trial of partially observed data set according to model (4.1)
- script_test_lin2dim_model_inner_crit.R:  estimate the subject-specific parameters b_i  from generated data via script_data_gen_model_lin_2dim.R
                                           for the true population parameters value (theta_star,Psi_star) by default (can be changed). 
- script_test_lin2dim_model_outer_crit.R:  estimate the mean population parameter theta_star as well as delta_star


For subsection 4.1.2:
- script_data_gen_model_lin_2dim_misspe.R:        generate trial of partially observed data set according to the pertubed model (4.2)
- script_test_lin2dim_misspe_model_inner_crit.R:  estimate the subject-specific parameters b_i  from generated data via script_data_gen_model_lin_2dim_misspe.R
                                                  for the true population parameters value (theta_star,Psi_star) by default (can be changed). 
- script_test_lin2dim_misspe_model_outer_crit.R:  estimate the mean population parameter theta_star as well as delta_star


For subsection 4.2.1:
- script_data_gen_insulin_model.R:         generate trial of partially observed data set according to model (4.3)
- script_test_insulin_model_inner_crit.R:  estimate the subject-specific parameters b_i  from generated data via script_data_gen_insulin_model.R
                                           for the true population parameters value (theta_star,Psi_star) by default (can be changed). 
- script_test_insulin_model_outer_crit.R:  estimate the mean population parameter theta_star as well as delta_star


For subsection 4.2.2:
- script_data_gen_insulin_model_misspe.R:         generate trial of partially observed data set according to the perturbed model (4.4)
- script_test_insulin_model_misspe_inner_crit.R:  estimate the subject-specific parameters b_i  from generated data via script_test_insulin_model_misspe_inner_crit.R
                                                  for the true population parameters value (theta_star,Psi_star) by default (can be changed). 
- script_test_insulin_model_misspe_outer_crit.R:  estimate the mean population parameter theta_star as well as delta_star


For subsection 4.3.1:
- script_data_gen_ebola_vaccine_model.R:              generate trial of observed data set according to model (4.6)
- script_test_ebola_vaccine_model_inner_crit.R:       estimate the subject-specific parameters b_i  from generated data via script_data_gen_ebola_vaccine_model.R
                                                      for the true population parameters value (theta_star,Psi_star) by default (can be changed). 
- script_test_ebola_vaccine_model_outer_crit.R:       estimate the mean population parameter theta_star as well as delta_star

For subsection 4.3.2:
- script_data_gen_ebola_vaccine_model_misspe.R:               generate trial of observed data set according to the perturbed model  (4.7)
- script_test_ebola_vaccine_model_misspe_inner_crit.R:        estimate the subject-specific parameters b_i  from generated data via script_data_gen_ebola_vaccine_model_misspe.R
                                                              for the true population parameters value (theta_star,Psi_star) by default (can be changed). 
- script_test_ebola_vaccine_model_misspe_outer_crit.R:        estimate the mean population parameter theta_star as well as delta_star
