
# -------------------------------------------------------------------------
#' Steffen Docken
#' 16-5-23
#' This code will run bootstrapping of the data and calculate 95% CIs. It also
#' visualizes the ln(Likelihood) function
#' 
# -------------------------------------------------------------------------

#' this code takes the infection data data frame with the desired predictor
#' variable defined as pred_var as an input
param_conf_bootstrap <-function(Infection_df, param.id, pred_var_name,
                                Model_ID_code){ 
  #bootstrapping animal data
  
  set.seed(3333) #setting the seed so bootstrapping is reproducible
  
# defining data structure -------------------------------------------------
  Infection_df = Infection_df %>% filter(!is.na(pred_var)) #only using
  # animals for whom current predictor variable is defined
  
  n_animals_per_arm = Infection_df %>% group_by(Arm) %>% summarise(n = n())
  #number of animals per study arm is in column n and Study Arm name is in column
  #Arm
  
  n_study_arms = length(n_animals_per_arm$Arm) #number of study arms
  
  bootstrapped_data_df = data.frame() #initializing data frame for bootstrapped 
  #data
  
# generate bootstrapped data ---------------------------------------------------------
  bootstrap_data_file = paste0("output/Bootstrap_Results/", pred_var_name, 
                               "_bootstrap_data.RData")
  
  if(!file.exists(bootstrap_data_file)){
    
    for (bootstrap_ind in 1:n_bootstrap) {
      bootstrapped_data_df_it = data.frame() #initializing data frame for this 
      #iteration of bootstrapped data
      
      for (ii in 1:n_study_arms) {
        study_arm_infection_data = Infection_df %>% 
          filter(Arm == n_animals_per_arm$Arm[ii])
        #collecting data for current study arm
        
        bootstrapped_ind_it = sample(x = n_animals_per_arm$n[ii], 
                                     size = n_animals_per_arm$n[ii],
                                     replace = TRUE)
        
        bootstrapped_data_df_it = rbind(bootstrapped_data_df_it,
                                        study_arm_infection_data[bootstrapped_ind_it,])
        
        
      }
      
      bootstrapped_data_df_it$boot_ind = bootstrap_ind #adding a variable to 
      #indicate which iteration of the bootstrapping is currently happening
      
      bootstrapped_data_df = rbind(bootstrapped_data_df,
                                   bootstrapped_data_df_it)
    }
    
    save(bootstrapped_data_df, file = bootstrap_data_file)
  }else{
    load(bootstrap_data_file)
  }
  
  
# fit model to bootstrapped data ------------------------------------------
  bootstrap_param_file = paste0("output/Bootstrap_Results/", pred_var_name,
                                "_", Model_ID_code, "_bootstrap_param.RData")
  
  if(!file.exists(bootstrap_param_file)){
  
    cl <- makeCluster(n_cores[1]-1) #not to overload your computer
    registerDoParallel(cl)
    
    Bootstrap_param_fit_df <- foreach(bootstrap_ind = 1:n_bootstrap, 
                                      .combine = rbind,
                                      .packages = c("tidyverse", "stats"),
                                      .export = c("n_LHS", "LHS_param_matrix",
                                                  "neg_logLike_infection",
                                                  "prob_infect")) %dopar% {
      
      bootstrapped_data_df_it = bootstrapped_data_df %>% filter(boot_ind == bootstrap_ind) %>% 
        mutate(infect_cens = if_else(week_of_infection <= 8, infect_cens, as.integer(0)),
               week_of_infection = if_else(week_of_infection <=8, week_of_infection, as.integer(8)))
      #censoring animals detected or censored after week 9 (this ensures we don't 
      #include animals detected more than a week after last inoculation, which
      #allows us to be consistent in our assumption that animals were infected
      #the week before detected infection.)
  
      neg_logLike_LHSvec = rep(1e10, n_LHS) #vector to record -logLikelihood for
      #each LHS parameter set
      
      for (kk in 1:n_LHS) {
        params_LHS_it = LHS_param_matrix[kk, param.id] #only using current parameters
        
        neg_logLike_LHSvec[kk] = neg_logLike_infection(params = params_LHS_it, 
                                                       Infection_df = bootstrapped_data_df_it,
                                                       param.id = param.id)
      }
      
      params_init = LHS_param_matrix[which.min(neg_logLike_LHSvec), param.id]
      
      for (opt_ind in 1:10) {
        neglogLike_min_list = nlm(neg_logLike_infection, params_init, 
                                  Infection_df = bootstrapped_data_df_it,
                                  param.id = param.id)
        
        params_init = neglogLike_min_list$estimate
        #browser()
      }
      
      opt_param_it = neglogLike_min_list$estimate #last parameter estimates
      
      Bootstrap_param_fit_df_it = 
        data.frame(x2 = 0, x = 0, c = 0, inv_x = 0, ln_x = 0, 
                   boot_ind = bootstrap_ind)
      
      Bootstrap_param_fit_df_it[,param.id] = opt_param_it
      
      Bootstrap_param_fit_df_it
    }
    
    stopCluster(cl) #stop the cluster
    
    save(Bootstrap_param_fit_df, file = bootstrap_param_file)
  }else{
    load(bootstrap_param_file)
  }
  
  Bootstrap_95_conf_int_df = 
    data.frame(x2 = unname(quantile(Bootstrap_param_fit_df$x2, probs = c(0.025, 0.975), type = 1)),
               x = unname(quantile(Bootstrap_param_fit_df$x, probs = c(0.025, 0.975), type = 1)),
               c = unname(quantile(Bootstrap_param_fit_df$c, probs = c(0.025, 0.975), type = 1)),
               inv_x = unname(quantile(Bootstrap_param_fit_df$inv_x, probs = c(0.025, 0.975), type = 1)),
               ln_x = unname(quantile(Bootstrap_param_fit_df$ln_x, probs = c(0.025, 0.975), type = 1)))
  
  Bootstrap_dfs = list(Bootstrap_param_fit_df, Bootstrap_95_conf_int_df)
  
  return(Bootstrap_dfs)
}


# log-likelihood profile plotting -----------------------------------------
#' This code takes the infection data data frame with the desired predictor
#' variable defined as pred_var as an input. It also takes the designator for 
#' parameters used in the model (Model_ID_code),
#' and the string name of the predictor variable (pred_var_name)
#' for look up of parameters
param_conf_loglikelihood <-function(Infection_df, Model_ID_code, pred_var_name){
  Model_fits <- read.csv('output/Survival_fit_outputs.csv') #loading model fit
  #parameters
  
  #getting parameters for this model
  params = filter(Model_fits, Pred_var == pred_var_name, 
                  Model_ID == Model_ID_code) %>%
    select(x2, x, c, inv_x, ln_x) %>% as.numeric()#getting 
  #best fit model params
  
  param.id = which(params != 0) #identifying non-zero parameter (i.e., parameters
  #included in the model)
  
  params = params[param.id] #removing non-used parameters and name
  
  n_param = length(param.id)
  
  param_names_it = c('x2', 'x', 'c', 'inv_x', 'ln_x')[param.id]
  
  logLike_profile_df = data.frame()
  logLike_opt_param_df = data.frame()
  norm_95CI_est_df = data.frame(param_name = rep(param_names_it, each = 2),
                                CI_lim_name = rep(c("lower", "upper"), times = n_param),
                                CI_lim = 0)
  
  for (ii in 1:n_param) {
    
    delta_param = rep(0, n_param)
    delta_param[ii] = params[ii]/100 #delta for estimating second 
    #derivative of log-likelihood is param value divided by 100. only the current
    #parameter is adjusted for calculation of second derivative
    
    d2_lnL = (-neg_logLike_infection(params = params + delta_param, 
                                    Infection_df = Infection_df,
                                    param.id = param.id) -
                2*(-neg_logLike_infection(params = params, 
                                        Infection_df = Infection_df,
                                        param.id = param.id)) +
                (-neg_logLike_infection(params = params - delta_param, 
                                        Infection_df = Infection_df,
                                        param.id = param.id)))/delta_param[ii]^2
    #estimation of second derivative of log-likelihood with respect to current 
    #parameter
    
    std_est = 1/(-d2_lnL)^(1/2) #estimation of standard error for parameter based
    #on assuming parameter estimates are normally distributed
    
    range_param = 2*1.96*std_est
    
    norm_95CI_est_df[norm_95CI_est_df$param_name == param_names_it[ii],]$CI_lim = 
      c(params[ii] - 1.96*std_est, params[ii] + 1.96*std_est)
    
    logLike_profile_df_it = data.frame(param_name = rep(param_names_it[ii], len = 101),
                                      param_val = seq(from =params[ii]-range_param,
                                        to = params[ii] + range_param,
                                        by = range_param/50),
                                      logLikelihood = rep(0, 101))
    
    for (jj in 1:101) {
      params_it = params
      params_it[ii] = logLike_profile_df_it$param_val[jj]
      
      logLike_profile_df_it$logLikelihood[jj] = 
        -neg_logLike_infection(params = params_it,
                               Infection_df = Infection_df,
                               param.id = param.id)
    }
    
    logLike_opt_param_df_it = filter(logLike_profile_df_it, param_val == params[ii])
    #collecting row of profile data frame corresponding to optimal parameter to 
    #be used in visualization plot
    
    logLike_profile_df = rbind(logLike_profile_df, logLike_profile_df_it) 
    
    logLike_opt_param_df = rbind(logLike_opt_param_df, logLike_opt_param_df_it)
    
  }
  
  logLike_profile_dfs = list(logLike_profile_df, logLike_opt_param_df, norm_95CI_est_df)
  
  return(logLike_profile_dfs)
}

