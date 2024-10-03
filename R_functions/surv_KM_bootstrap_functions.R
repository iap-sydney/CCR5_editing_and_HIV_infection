
# -------------------------------------------------------------------------
#' Steffen Docken
#' 29-5-23
#' Code to generate bootstrapped data for survival curves
#' 
# -------------------------------------------------------------------------

surv_95PI_range <- function(Infection_df, pred_var_name,
                            Model_ID_code){

  
  # defining data structure -------------------------------------------------
  Infection_df = Infection_df %>% filter(!is.na(pred_var)) #only using
  # animals for whom current predictor variable is defined
  
  n_animals_per_arm = Infection_df %>% group_by(Arm) %>% summarise(n = n())
  #number of animals per study arm is in column n and Study Arm name is in column
  #Arm
  
  n_study_arms = length(n_animals_per_arm$Arm) #number of study arms
  
# parametric bootstrapping -----------------------------------------------------------
  bootstrap_surv_data_file = paste0("output/Bootstrap_Results/", pred_var_name, 
                                   "_", Model_ID_code, "_bootstrap_surv_data.RData")
  
  if(!file.exists(bootstrap_surv_data_file)){
    param_bootstrap_data_df = data.frame()
    surv_curve_bootstrap = data.frame()
    
    #importing bootstrapped parameters
    bootstrap_param_file = paste0("output/Bootstrap_Results/", pred_var_name,
                                  "_", Model_ID_code, "_bootstrap_param.RData")
    
    load(bootstrap_param_file)
    
    #randomizing picking of parameter sets
    set.seed(3333)
    param_boot_ind = sample(x = length(Bootstrap_param_fit_df[,1]),
                            size = n_bootstrap, replace = TRUE)
    
    set.seed(2222) #setting the seed so bootstrapping is reproducible
    
    #for (boot_ind in 1:n_bootstrap) {
    cl <- makeCluster(n_cores[1]-1) #not to overload your computer
    registerDoParallel(cl)
    
    param_bootstrap_data_df <- foreach(boot_ind = 1:n_bootstrap, 
                                      .combine = rbind,
                                      .packages = c("tidyverse", "stats"),
                                      .export = c("prob_infect",
                                                  "prob_infect_week")) %dopar% {  
      
      # getting parameters for this iteration model
      params = Bootstrap_param_fit_df[param_boot_ind[boot_ind],]
      
      param.id = which(params != 0) #identifying non-zero parameter (i.e., parameters
      #included in the model)
      
      params = params[param.id] #removing non-used parameters and name
      
      n_param = length(param.id)
      
      param_bootstrap_data_df_it = Infection_df %>% rowwise () %>% 
        mutate(Animal = Animal,
               Arm = Arm,
               bootstrap = boot_ind,
               prob_inoc_infect =
                 prob_infect(params = as.numeric(params), x = pred_var, param.id = param.id),
               inf_wk1 = prob_infect_week(1, prob_inoc_infect, 1),
               inf_wk2 = inf_wk1 + prob_infect_week(2, prob_inoc_infect, 1),
               inf_wk3 = inf_wk2 + prob_infect_week(3, prob_inoc_infect, 1),
               inf_wk4 = inf_wk3 + prob_infect_week(4, prob_inoc_infect, 1),
               inf_wk5 = inf_wk4 + prob_infect_week(5, prob_inoc_infect, 1),
               inf_wk6 = inf_wk5 + prob_infect_week(6, prob_inoc_infect, 1),
               inf_wk7 = inf_wk6 + prob_infect_week(7, prob_inoc_infect, 1),
               inf_wk8 = inf_wk7 + prob_infect_week(8, prob_inoc_infect, 1),
               undet = 1,
               prob = runif(1),
               .keep = "none")
      
      param_bootstrap_data_df_it = param_bootstrap_data_df_it %>%
        mutate(week_of_infection = which(prob < c_across(inf_wk1:undet))[1],
               infect_cens = as.numeric(week_of_infection < 9)) %>% #if <9, 
        #infection detected, else, censoring
        mutate(week_of_infection = min(week_of_infection, 8)) %>% ungroup() 
      #changing undetected timing to week 8 generating the first bootstrap
      
      param_bootstrap_data_df_it
                                                  
      }
      # Survival curve for first bootstrapped data ------------------------------------
      surv_curve_bootstrap <- foreach(boot_ind = 1:n_bootstrap, 
                                       .combine = rbind,
                                       .packages = c("tidyverse", "stats",
                                                     "survival", "survminer")) %dopar% {
                                         
      surv_curve_bootstrap_full_it = 
        survfit(Surv(week_of_infection, infect_cens) ~ Arm, 
                data = filter(param_bootstrap_data_df, bootstrap == boot_ind))
      
      surv_curve_bootstrap_it = data.frame(week_of_infection = surv_curve_bootstrap_full_it$time,
                                                surv = surv_curve_bootstrap_full_it$surv,
                                                lower = surv_curve_bootstrap_full_it$lower,
                                                upper = surv_curve_bootstrap_full_it$upper,
                                                Arm = c(rep(names(surv_curve_bootstrap_full_it$strata[1]), surv_curve_bootstrap_full_it$strata[1]),
                                                        rep(names(surv_curve_bootstrap_full_it$strata[2]), surv_curve_bootstrap_full_it$strata[2]),
                                                        rep(names(surv_curve_bootstrap_full_it$strata[3]), surv_curve_bootstrap_full_it$strata[3]),
                                                        rep(names(surv_curve_bootstrap_full_it$strata[4]), surv_curve_bootstrap_full_it$strata[4]),
                                                        rep(names(surv_curve_bootstrap_full_it$strata[5]), surv_curve_bootstrap_full_it$strata[5])),
                                                bootstrap = rep(boot_ind, length(surv_curve_bootstrap_full_it$time)))
      
      surv_curve_bootstrap_it = rbind(surv_curve_bootstrap_it,
                                      data.frame(week_of_infection = rep(0, 5),
                                                 surv = rep(1,5),
                                                 lower = rep(1, 5),
                                                 upper = rep(1, 5),
                                                 Arm = names(surv_curve_bootstrap_full_it$strata),
                                                 bootstrap = rep(boot_ind, 5))) %>%
        arrange(Arm, week_of_infection)
      
      for (Arm_name in names(surv_curve_bootstrap_full_it$strata)) {
        for (week_ind in 0:8) {
          if(nrow(filter(surv_curve_bootstrap_it, 
                           Arm == Arm_name,
                           week_of_infection == week_ind)) == 0) { #no data for 
            #this week/arm combination
            prev_row_index = which((surv_curve_bootstrap_it$Arm == Arm_name)&
                                     (surv_curve_bootstrap_it$week_of_infection == week_ind-1))
            
            surv_curve_bootstrap_it = surv_curve_bootstrap_it %>% 
              add_row(week_of_infection = week_ind,
                      surv = surv_curve_bootstrap_it$surv[prev_row_index],
                      lower = surv_curve_bootstrap_it$lower[prev_row_index],
                      upper = surv_curve_bootstrap_it$upper[prev_row_index],
                      Arm = surv_curve_bootstrap_it$Arm[prev_row_index],
                      bootstrap = surv_curve_bootstrap_it$bootstrap[prev_row_index], 
                      .after = prev_row_index)
            #adding rows for weeks where events didn't happen in bootstrapped data
            
          }
        }
      }
      
      surv_curve_bootstrap_it
                                       
                                       }
    stopCluster(cl)
    
    
    save(param_bootstrap_data_df,
         surv_curve_bootstrap, file = bootstrap_surv_data_file)
  }else{
    load(bootstrap_surv_data_file)
  }
  
  surv_95PI_range_df = data.frame(week = rep(0:8, times = n_study_arms),
                                  Study_Arm = rep(unique(surv_curve_bootstrap$Arm), each = 9)) %>%
    rowwise() %>% mutate(lower = unname(quantile(filter(surv_curve_bootstrap, 
                                                        Arm == Study_Arm, 
                                                        week_of_infection == week)$surv, 
                                                 probs = 0.025, type = 1)),
                         upper = unname(quantile(filter(surv_curve_bootstrap, 
                                                        Arm == Study_Arm, 
                                                        week_of_infection == week)$surv, 
                                                 probs = 0.975, type = 1))) %>%
    mutate(Arm = Study_Arm, .keep = 'unused') %>% ungroup()

  output_list = list(param_bootstrap_data_df, surv_curve_bootstrap,
                     surv_95PI_range_df)
  
  return(output_list)
  
}