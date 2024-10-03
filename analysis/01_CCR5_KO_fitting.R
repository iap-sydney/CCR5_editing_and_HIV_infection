#Steffen Docken
#27-1-23
#This code will fits models for how effective various levels of CCR5 KO are in 
#protecting humanized mice from infection

Modeling_results_df = data.frame(Pred_var = rep('', num_pred_var*num_model_forms+2),
                                 Model_ID = rep(0, num_pred_var*num_model_forms+2),
                                 num_animals = rep(0, num_pred_var*num_model_forms+2),
                                 x2 = rep(0, num_pred_var*num_model_forms+2), #coefficient for x^2
                                 x = rep(0, num_pred_var*num_model_forms+2),#coefficient for x
                                 c = rep(0, num_pred_var*num_model_forms+2),#constant term
                                 inv_x = rep(0, num_pred_var*num_model_forms+2),#R in the (Rx-1)/(Rx) 
                                 #term of the proliferation risk model
                                 ln_x = rep(0, num_pred_var*num_model_forms+2),#coefficient for ln(x)
                                 max_logLike = rep(0, num_pred_var*num_model_forms+2),
                                 AIC = rep(0, num_pred_var*num_model_forms+2))

week_infect_vec = 0:8 #vector of weeks of inoculation for future use

Animal_wFull_infection_data_df = Animal_infection_data_df %>% 
  filter(!is.na(BM_CD19_CCR5_KO_mean)) #only including
#animals that have data for all measurements being considered as predictors
#browser()

Infection_df = Animal_wFull_infection_data_df %>% 
  mutate(Animal = Animal,
         week_of_infection = week_of_infection,
         infect_cens = infect_cens,
         Arm = Arm,
         pred_var = 1-BM_CD19_CCR5_KO_mean/100,#converting to fraction
                           #with CCR5
         max_pVL = max_pVL,
         .keep = "none")

Pred_var_it = 'CD19_fracCCR5'

#browser()
summary_Infection_df = Infection_df %>% group_by(Arm) %>% 
  summarise(pred_var = mean(pred_var)) 
#data frame with means of pred_var to be used in imaging of fit.

for (jj in 1:num_model_forms) {
  param.id = switch(jj, c(2), # rho = bx
                    c(2, 3), # rho = bx + c
                    c(1), # rho = ax^2
                    c(1, 2), # rho = ax^2 + bx
                    c(1, 2, 3), # rho = ax^2 + bx + c
                    c(2, 3, 4), # rho = (bx + c)*(Rx-1)/Rx
                    c(1, 2, 3, 4), # rho = (ax^2 + bx + c)*(Rx-1)/(Rx)
                    c(5), # rho = f*ln(x)
                    c(3, 5)) # rho = c + f*ln(x)
  
  neg_logLike_LHSvec = rep(1e10, n_LHS) #vector to record -logLikelihood for
  #each LHS parameter set
  
  for (kk in 1:n_LHS) {
    params_LHS_it = LHS_param_matrix[kk, param.id] #only using current parameters
    
    neg_logLike_LHSvec[kk] = neg_logLike_infection(params = params_LHS_it, 
                                                Infection_df = Infection_df,
                                                param.id = param.id)
    
  }

  params_init = LHS_param_matrix[which.min(neg_logLike_LHSvec), param.id]
  #browser()
  for (opt_ind in 1:10) {
    neglogLike_min_list = nlm(neg_logLike_infection, params_init, 
                           Infection_df = Infection_df,
                           param.id = param.id)
    
    params_init = neglogLike_min_list$estimate
    #browser()
  }
  
  opt_param_it = neglogLike_min_list$estimate #last parameter estimates
  
  opt_param_full = rep(0, 5)
  opt_param_full[param.id] = opt_param_it
  #browser()
  #recording results:
  param.id.record = switch(jj, 2, # rho = bx
                    23, # rho = bx + c
                    1, # rho = ax^2
                    12, # rho = ax^2 + bx
                    123, # rho = ax^2 + bx + c
                    234, # rho = (bx + c)*(Rx-1)/Rx
                    1234, # rho = (ax^2 + bx + c)*(Rx-1)/Rx
                    5, # rho = d*ln(x)
                    35) # rho = c + d*ln(x)
  
  Modeling_results_df$Pred_var[jj] = Pred_var_it
  Modeling_results_df[jj,2:length(Modeling_results_df[1,])] = 
    c(param.id.record, length(Infection_df$Animal), 
      opt_param_full[1], opt_param_full[2], opt_param_full[3], opt_param_full[4],
      opt_param_full[5], -neglogLike_min_list$minimum,
      AIC_generic(lnL = -neglogLike_min_list$minimum,
                  k = length(param.id),
                  n = length(Infection_df$Animal),
                  Corrected = 1))# calculating the corrected AIC
  
  print(glue::glue("{Pred_var_it}: {jj}")) #indicate status of code.

}

fit <- survfit(Surv(week_of_infection, infect_cens) ~ Arm,
                         data = Infection_df)

survival_plot_it <- ggsurvplot(fit, data = Infection_df, legend = "right")

print(survival_plot_it)

raw_data_plot_it = ggplot(Infection_df, 
                          aes(x = week_of_infection, y = pred_var, 
                              color = Arm, shape = factor(infect_cens))) + 
  geom_point() + scale_shape_manual(values = c(1, 16))+ 
  labs(y= Pred_var_it, title='', shape = 'Infection',
       x = 'Week of Infection')

print(raw_data_plot_it)


raw_data_v_pVL_plot_it = ggplot(Infection_df, 
                          aes(x = max_pVL, y = pred_var, 
                              color = Arm, shape = factor(infect_cens))) + 
  geom_point() + scale_shape_manual(values = c(1, 16))+ 
  labs(y= Pred_var_it, title='', shape = 'Infection',
       x = 'max pVL (log_10 copies/ml)')

print(raw_data_v_pVL_plot_it)


## fitting survival just by Arm
Infection_df = Animal_wFull_infection_data_df %>% 
  mutate(Animal = Animal,
         week_of_infection = week_of_infection,
         infect_cens = infect_cens,
         Arm = Arm,
         .keep = "none")

opt_Arm_param = fit_infection_prob_Arm(Infection_df = Infection_df)

logLike_max_Arm = -neg_logLike_infection_Arm(params = opt_Arm_param,
                                             Infection_df = Infection_df)

Modeling_results_df$Pred_var[num_pred_var*num_model_forms+1] = "Arm"
Modeling_results_df[num_pred_var*num_model_forms+1,2:length(Modeling_results_df[1,])] = 
  c(0, length(Infection_df$Animal), 
    opt_Arm_param[1], opt_Arm_param[2], opt_Arm_param[3], opt_Arm_param[4],
    opt_Arm_param[5], logLike_max_Arm,
    AIC_generic(lnL = logLike_max_Arm,
                k = length(opt_Arm_param),
                n = length(Infection_df$Animal),
                Corrected = 1))# calculating the corrected AIC


print("Arm") #indicate status of code.

## fitting overall average infection rate
Infection_df = Animal_wFull_infection_data_df %>% 
  mutate(Animal = Animal,
         week_of_infection = week_of_infection,
         infect_cens = infect_cens,
         Arm = Arm,
         .keep = "none")

opt_avg_param = fit_infection_prob_all(Infection_df = Infection_df)

logLike_max_all = -neg_logLike_infection_Arm(params = opt_avg_param,
                                             Infection_df = Infection_df)
#arm specific log-likelihood function works with all arms having the same 
#infection risk

Modeling_results_df$Pred_var[num_pred_var*num_model_forms+2] = "Overall Average"
Modeling_results_df[num_pred_var*num_model_forms+2,2:length(Modeling_results_df[1,])] = 
  c(0, length(Infection_df$Animal), 
    opt_avg_param[1], opt_avg_param[2], opt_avg_param[3], opt_avg_param[4],
    opt_avg_param[5], logLike_max_all,
    AIC_generic(lnL = logLike_max_all,
                k = 1,
                n = length(Infection_df$Animal),
                Corrected = 1))# calculating the corrected AIC


print("Overall Average") #indicate status of code.

Modeling_results_df = Modeling_results_df %>% mutate(delta_AIC = AIC - min(AIC))
#calculating delta AIC

write.csv(Modeling_results_df, file = 'output/Survival_fit_outputs.csv')

