
# -------------------------------------------------------------------------
#' Steffen Docken
#' 30-3-23
#' Code to examine bootstrapping of survival curves based on best fit models
#' 
# -------------------------------------------------------------------------

Model_fits <- read.csv('output/Survival_fit_outputs.csv')

frac_CCR5_vec = seq(from = 0, by = 0.01, to = 1)

pred_var_name = 'CD19_fracCCR5'

for (ii in c(1,3)) {
  
  Model_ID_code = switch(ii, 2, 23, 234)
  param.id = switch(ii, c(2), c(2,3), c(2,3,4))
  
  Model_name = switch(ii, pred_var_name,
                      paste0(pred_var_name, '_', Model_ID_code),
                      paste0(pred_var_name, '_prolif_model'))
  print(Model_name)
  
  pred_var_vec = frac_CCR5_vec
  pred_var_axis_name = expression(paste("% ", italic("CCR5"), "-edited"))
  bin_size = max(frac_CCR5_vec)/5
  
  
  params = filter(Model_fits, Pred_var == pred_var_name, 
                  Model_ID == Model_ID_code) %>%
    select(x2, x, c, inv_x, ln_x)#getting 
  #best fit model params
  
  if (params[4] != 0){
    print(paste0('The threshold to infection is ', 1/params[4], 
                 ' ', pred_var_name))
  }
  
  Infection_df = Animal_infection_data_df %>% #getting predictor variable info
    mutate(Animal = Animal,
           week_of_infection = week_of_infection,
           infect_cens = infect_cens,
           Arm = Arm,
           pred_var = 1-BM_CD19_CCR5_KO_mean/100,#converting to fraction with CCR5
           .keep = "none") %>% filter(!is.na(pred_var))
  
  Infection_df_summary = Infection_df %>% group_by(Arm) %>% 
    summarise(n = n(),
              mean_pred_var = mean(pred_var),
              median_pred_var = median(pred_var)) %>% ungroup()
  
  print(Infection_df_summary)
  
  Infection_df = Infection_df %>% rowwise() %>% #getting probability cut offs for each potential
    #week of infection
    mutate(prob_inoc_infect =
             prob_infect(params = as.numeric(params), x = pred_var, param.id = 1:5), #use all 
           #params because the full parameter set pulled in contains appropriate 0s
           inf_wk1 = prob_infect_week(1, prob_inoc_infect, 1),
           inf_wk2 = inf_wk1 + prob_infect_week(2, prob_inoc_infect, 1),
           inf_wk3 = inf_wk2 + prob_infect_week(3, prob_inoc_infect, 1),
           inf_wk4 = inf_wk3 + prob_infect_week(4, prob_inoc_infect, 1),
           inf_wk5 = inf_wk4 + prob_infect_week(5, prob_inoc_infect, 1),
           inf_wk6 = inf_wk5 + prob_infect_week(6, prob_inoc_infect, 1),
           inf_wk7 = inf_wk6 + prob_infect_week(7, prob_inoc_infect, 1),
           inf_wk8 = inf_wk7 + prob_infect_week(8, prob_inoc_infect, 1),
           undet = 1)%>% ungroup()
  
  n_animals = length(Infection_df$Animal)
  
  # Survival curve for data -------------------------------------------------
  surv_fit_data = survfit(Surv(week_of_infection, infect_cens) ~ Arm, 
                          data = Infection_df)
  
  surv_fit_data_extracted = data.frame(week_of_infection = surv_fit_data$time,
                                       surv = surv_fit_data$surv,
                                       lower = surv_fit_data$lower,
                                       upper = surv_fit_data$upper,
                                       Arm = c(rep(names(surv_fit_data$strata[1]), surv_fit_data$strata[1]),
                                               rep(names(surv_fit_data$strata[2]), surv_fit_data$strata[2]),
                                               rep(names(surv_fit_data$strata[3]), surv_fit_data$strata[3]),
                                               rep(names(surv_fit_data$strata[4]), surv_fit_data$strata[4]),
                                               rep(names(surv_fit_data$strata[5]), surv_fit_data$strata[5])))
  
  surv_fit_data_extracted = rbind(surv_fit_data_extracted,
                                  data.frame(week_of_infection = rep(0, 5),
                                             surv = rep(1,5),
                                             lower = rep(1, 5),
                                             upper = rep(1, 5),
                                             Arm = names(surv_fit_data$strata))) %>%
    mutate(Arm = factor(substring(Arm, 5), levels = c("0% CCR5-KO", "25% CCR5-KO", 
                                                        "50% CCR5-KO", "75% CCR5-KO", 
                                                        "100% CCR5-KO")))
  
  p_surv_data = ggplot(surv_fit_data_extracted, aes(x = week_of_infection,
                                                    y = surv, group = Arm,
                                                    color = Arm)) + 
    geom_step() + theme_classic() + labs(x = "Inoculation",
                                         y = "Remaining Fraction Uninfection") +
    Surv_arm_palette
  
  print(p_surv_data)

# Comparing KM curves for increasing %CCR5 KO -----------------------------

  for (jj in 1:4) {
    Arm1 = switch(jj, "0% CCR5-KO", "25% CCR5-KO", "50% CCR5-KO",
                  "75% CCR5-KO")
    Arm2 = switch(jj, "25% CCR5-KO", "50% CCR5-KO",
                  "75% CCR5-KO", "100% CCR5-KO")
    
    surv_fit_2arm = survfit(Surv(week_of_infection, infect_cens) ~ Arm, 
                              data = filter(Infection_df, 
                                            Arm == Arm1|Arm == Arm2))
    print(paste0(Arm1, " vs. ", Arm2))
    print(surv_pvalue(surv_fit_2arm))
    
  }

# bootstrap parameter confidence intervals --------------------------------

  Bootstrap_dfs = 
    param_conf_bootstrap(Infection_df = select(Infection_df, Animal,
                                               week_of_infection, infect_cens,
                                               Arm, pred_var),
                         param.id = param.id,
                         pred_var_name = pred_var_name, 
                         Model_ID_code = Model_ID_code)
  
  Bootstrap_param_fit_df = Bootstrap_dfs[[1]]
  Bootstrap_95_conf_int_df = Bootstrap_dfs[[2]]
  
  logLike_profile_dfs = 
    param_conf_loglikelihood(Infection_df = select(Infection_df, Animal,
                                                   week_of_infection, infect_cens,
                                                   Arm, pred_var),
                             Model_ID_code = Model_ID_code,
                             pred_var_name = pred_var_name)
  
  logLike_profile_df = logLike_profile_dfs[[1]]
  logLike_opt_param_df = logLike_profile_dfs[[2]]
  norm_95CI_est_df = logLike_profile_dfs[[3]]
  
  Bootstrap_95_conf_int_df_profile = Bootstrap_95_conf_int_df %>% 
    select(all_of(param.id)) %>% mutate(CI_lim_name = c("lower", "upper")) %>%
    melt(id.vars = "CI_lim_name", variable.name = "param_name", value.name = "CI_lim")

  write.csv(Bootstrap_95_conf_int_df_profile, 
            file = paste0('output/', Model_name, '_bootstrap_param_CIs.csv'))
    
  
  
  
  # Getting 95% CI for survival curves from parametric bootstrapping --------
  surv_95PI_range_list = surv_95PI_range(Infection_df = select(Infection_df, Animal,
                                                               week_of_infection, infect_cens,
                                                               Arm, pred_var), 
                                         pred_var_name = pred_var_name,
                                         Model_ID_code = Model_ID_code)
  
  param_bootstrap_df = surv_95PI_range_list[[1]]
  surv_bootstrap_df = surv_95PI_range_list[[2]]
  surv_95PI_range_df = surv_95PI_range_list[[3]]
  
  # plot to compare fraction of animals infected by week 8 to the prediction -----------------------------------------------------------------
  frac_infect_df = Infection_df %>% rowwise() %>%
    mutate(pred_var_bin = min(floor(pred_var/bin_size)*bin_size + bin_size/2,
                              max(pred_var_vec)-bin_size/2)) %>%#labeling
  # bins for animals by % with CCR5
    group_by(pred_var_bin) %>% summarise(frac_infect = mean(infect_cens),
                                         num_infect = sum(infect_cens),
                                         n = n()) %>% rowwise() %>% mutate(
                                           plot_var = 100*(1-pred_var_bin),
                                           frac_infect_CI_low = 0,
                                           frac_infect_CI_high = 0,
                                           prob_inoc_infect =
                                             prob_infect(params = as.numeric(params), x = pred_var_bin, param.id = 1:5),#use all 
                                           #params because the full parameter set pulled in contains appropriate 0s
                                           inf_wk1 = prob_infect_week(1, prob_inoc_infect, 1),
                                           inf_wk2 = inf_wk1 + prob_infect_week(2, prob_inoc_infect, 1),
                                           inf_wk3 = inf_wk2 + prob_infect_week(3, prob_inoc_infect, 1),
                                           inf_wk4 = inf_wk3 + prob_infect_week(4, prob_inoc_infect, 1),
                                           inf_wk5 = inf_wk4 + prob_infect_week(5, prob_inoc_infect, 1),
                                           inf_wk6 = inf_wk5 + prob_infect_week(6, prob_inoc_infect, 1),
                                           inf_wk7 = inf_wk6 + prob_infect_week(7, prob_inoc_infect, 1),
                                           pred_frac_infect_by_wk_8 = inf_wk7 + prob_infect_week(8, prob_inoc_infect, 1), #Expected fraction infected
                                           pred_n = pred_frac_infect_by_wk_8*n) %>%
    ungroup()#Expected number infected
  
  for (jj in 1:length(frac_infect_df$pred_var_bin)) {
    frac_infect_df_it = frac_infect_df[jj,] #only getting entries for this row
    prop_test_output = prop.test(x = frac_infect_df_it$num_infect, 
                                 n = frac_infect_df_it$n,
                                 conf.level = 0.95)
    
    frac_infect_df$frac_infect_CI_low[jj] = prop_test_output$conf.int[1]
    frac_infect_df$frac_infect_CI_high[jj] = prop_test_output$conf.int[2]
  }
  
  Expect_infect_df = data.frame(pred_var = pred_var_vec)
  
  Expect_infect_df = Expect_infect_df %>% rowwise() %>%
    mutate(plot_var = 100*(1-pred_var),
      prob_inoc_infect =
             prob_infect(params = as.numeric(params), x = pred_var, param.id = 1:5),#use all 
           #params because the full parameter set pulled in contains appropriate 0s
           inf_wk1 = prob_infect_week(1, prob_inoc_infect, 1),
           inf_wk2 = inf_wk1 + prob_infect_week(2, prob_inoc_infect, 1),
           inf_wk3 = inf_wk2 + prob_infect_week(3, prob_inoc_infect, 1),
           inf_wk4 = inf_wk3 + prob_infect_week(4, prob_inoc_infect, 1),
           inf_wk5 = inf_wk4 + prob_infect_week(5, prob_inoc_infect, 1),
           inf_wk6 = inf_wk5 + prob_infect_week(6, prob_inoc_infect, 1),
           inf_wk7 = inf_wk6 + prob_infect_week(7, prob_inoc_infect, 1),
           prob_infect_by_wk_8 = inf_wk7 + prob_infect_week(8, prob_inoc_infect, 1))
  
  fit_quality_plot = ggplot(data = frac_infect_df, 
                            aes(x = plot_var,
                                y = frac_infect)) +
    geom_bar(stat = "identity") + geom_errorbar(aes(ymin = frac_infect_CI_low,
                                                    ymax = frac_infect_CI_high)) + 
    geom_line(data = Expect_infect_df,
                           aes(x = plot_var, 
                               y = prob_infect_by_wk_8),
                           inherit.aes = FALSE) + theme_classic() +
    scale_y_continuous(labels = percent) +
    labs(y='% HIV-infected', title='Probability of infection within 8 challenges',
         x = pred_var_axis_name) + 
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  
  print(fit_quality_plot)
  
  ggsave(paste0('output/Figures/Fig7G_', Model_name, '_Fit_quality.eps'), fit_quality_plot, width = 8,
         height = 8, units = 'cm')
  
  identity_line_df = data.frame(x_frac = c(0,1),
                                y_frac = c(0,1),
                                x_num = c(0, 25),
                                y_num = c(0, 25))
  
  fit_quality_percent = ggplot(data = frac_infect_df, 
                               aes(x = pred_frac_infect_by_wk_8,
                                   y = frac_infect)) +
    geom_point() + geom_line(data = identity_line_df, 
                             aes(x = x_frac, y = y_frac), inherit.aes = FALSE) + 
    theme_classic() + scale_y_continuous(labels = percent) + 
    scale_x_continuous(labels = percent) +
    labs(y = 'Percent infected within 8 weeks',
         x = 'Expected Percent infected within 8 weeks')
  
  print(fit_quality_percent)
  
  fit_quality_n = ggplot(data = frac_infect_df, 
                               aes(x = pred_n,
                                   y = num_infect)) +
    geom_point() + geom_line(data = identity_line_df, 
                             aes(x = x_num, y = y_num), inherit.aes = FALSE) + 
    theme_classic() + labs(y ='Number infected within 8 weeks',
         x = 'Expected Number infected within 8 weeks')
  
  print(fit_quality_percent)
  

# Survival Curve ----------------------------------------------------------
  
  Arm_vec = unique(Infection_df$Arm)
  
  Infection_df = Infection_df %>% rowwise() %>% #getting probability cut offs for each potential
    #week of infection
    mutate(prob_inoc_infect_95low =
             prob_infect(params = as.numeric(Bootstrap_95_conf_int_df[1,]), x = pred_var, param.id = 1:5), 
           prob_inoc_infect_95high =
             prob_infect(params = as.numeric(Bootstrap_95_conf_int_df[2,]), x = pred_var, param.id = 1:5)) %>%
    ungroup()
  
  Expect_surv_df = data.frame(week = rep(0:8, times = length(Arm_vec)),
                              Arm = rep(Arm_vec, each = 9),
                              prob_surv = rep(0, times = 9*length(Arm_vec)),
                              prob_surv_95low = rep(0, times = 9*length(Arm_vec)),
                              prob_surv_95high = rep(0, times = 9*length(Arm_vec)))
  
  for (jj in 1:(9*length(Arm_vec))) {
    if (Expect_surv_df$week[jj] == 0) {
      Expect_surv_df$prob_surv[jj] = 1
      Expect_surv_df$prob_surv_95low[jj] = 1
      Expect_surv_df$prob_surv_95high[jj] = 1
    }else {
      prob_inoc_infect_vec = filter(Infection_df, Arm == Expect_surv_df$Arm[jj])$prob_inoc_infect
      Expect_surv_df$prob_surv[jj] = sum(prob_infect_week(week_of_infection =  Expect_surv_df$week[jj],
                           prob_inoc_infect= prob_inoc_infect_vec, 
                           infect_cens = 0))/length(prob_inoc_infect_vec)
      
      prob_inoc_infect_95low_vec = filter(Infection_df, Arm == Expect_surv_df$Arm[jj])$prob_inoc_infect_95low
      Expect_surv_df$prob_surv_95low[jj] = sum(prob_infect_week(week_of_infection =  Expect_surv_df$week[jj],
                                                          prob_inoc_infect= prob_inoc_infect_95low_vec, 
                                                          infect_cens = 0))/length(prob_inoc_infect_95low_vec)
      
      prob_inoc_infect_95high_vec = filter(Infection_df, Arm == Expect_surv_df$Arm[jj])$prob_inoc_infect_95high
      Expect_surv_df$prob_surv_95high[jj] = sum(prob_infect_week(week_of_infection =  Expect_surv_df$week[jj],
                                                          prob_inoc_infect= prob_inoc_infect_95high_vec, 
                                                          infect_cens = 0))/length(prob_inoc_infect_95high_vec)
    }
    
  }
  
  
  Expect_surv_plot_it = ggplot(data = Expect_surv_df, 
                               aes(x = week, y = prob_surv,
                                   color = Arm)) +
    geom_line() + labs(y='Suvival Probability',
                       title='', color = 'Arm')
  print(Expect_surv_plot_it)
  
  surv_95PI_range_step_df = filter(surv_95PI_range_df, week !=8) %>%
    mutate(week = week + 0.999) #creating additional points for step function, 
  #not adding a full 1 to make sure plotting looks right.
  
  surv_95PI_range_plot_df = rbind(surv_95PI_range_df, 
                                  surv_95PI_range_step_df) %>%
    mutate(Arm = factor(substring(Arm, 5), levels = c("0% CCR5-KO", "25% CCR5-KO", 
                                                      "50% CCR5-KO", "75% CCR5-KO", 
                                                      "100% CCR5-KO")))
  
  surv_95PI_range_plot_df = 
    surv_95PI_range_plot_df[order(surv_95PI_range_plot_df$Arm,
                                  surv_95PI_range_plot_df$week,
                                  -surv_95PI_range_plot_df$upper),]
  #generating data frame that will have 2 measurements for each week to get step
  #form for 95%CI ribbon
  
  
  Expect_surv_plot_comb_noband = ggplot(surv_fit_data_extracted, 
                                 aes(x = week_of_infection, y = surv, group = Arm,
                                     color = Arm)) + 
    geom_step() + 
    geom_line(data = Expect_surv_df, aes(x = week, y = prob_surv,
                                         color = Arm),
              linetype = "dashed", inherit.aes = FALSE) + 
    scale_y_continuous(labels = percent) + scale_x_continuous(breaks = 0:8) +
    Surv_arm_palette + theme_classic() + labs(x = "HIV challenge #",
                                              y = "% HIV-negative")
  
  ggsave(paste0('output/Figures/Fig7E_', Model_name, 
                '_surv_plot.eps'), 
         Expect_surv_plot_comb_noband + theme(legend.position = 'none'), 
         width = 8, height = 8, units = 'cm')
  
  Expect_surv_plot_comb2 = ggplot(surv_fit_data_extracted, aes(x = week_of_infection,
                                                              y = surv, group = Arm,
                                                              color = Arm)) + 
    geom_step() + 
    geom_ribbon(data = surv_95PI_range_plot_df, inherit.aes = FALSE, 
                aes(x = week, ymin = lower, ymax = upper, 
                    fill = as.factor(Arm)), alpha = 0.2)+ 
    geom_line(data = Expect_surv_df, aes(x = week, y = prob_surv,
                                         color = Arm),
              linetype = "dashed", inherit.aes = FALSE) + 
    scale_y_continuous(labels = percent) + scale_x_continuous(breaks = 0:8) +
    Surv_arm_palette +
    theme_classic() + labs(x = "Inoculation",
                                              y = "Percent Remaining Uninfected")
  
  
  print(Expect_surv_plot_comb2)
  
  ggsave(paste0('output/Figures/', Model_name, '_PI_Surv_plot_230419_wLegend.pdf'),
         Expect_surv_plot_comb2, width = 10, height = 8, units = 'cm')
  
  # New facet label names for %CCR5-KO
  Arm_labs <- c("0% CCR5-KO", "26% CCR5-KO", "54% CCR5-KO",
                "69% CCR5-KO","96% CCR5-KO")
  names(Arm_labs) <- c("0% CCR5-KO", "25% CCR5-KO", "50% CCR5-KO",
                       "75% CCR5-KO","100% CCR5-KO")
  
  ggsave(paste0('output/Figures/', Model_name, '_PI_Surv_plot_facet.pdf'),
         Expect_surv_plot_comb2 + theme(legend.position = 'none') +
           facet_wrap(~Arm, labeller = labeller(Arm = Arm_labs)), width = 16,
         height = 12, units = 'cm')
  
# Probability of infection and number of inoculation needed. --------------
  
  Expect_infect_df = Expect_infect_df %>% 
    mutate(inoc_50 = ceiling(log(0.5)/log(1-prob_inoc_infect)))
  
  Expect_infect_df$inoc_50[Expect_infect_df$prob_inoc_infect == 0] =Inf
  #setting number of inoculations needed to Infty (numerical calculation sets it
  #to -Infty)
  
  prop_const = 1/4
  
  p_prob_infect_v_KO_plot2 = ggplot(data = Expect_infect_df,
                                   aes(x = plot_var))+
    geom_line(aes(y = 100*prob_inoc_infect), color = "red") +
    geom_line(aes(y = inoc_50*prop_const), color = "blue") +
    scale_y_continuous(name = "Risk of infection per challenge",
                       sec.axis = sec_axis(~./prop_const, name = "Challenges to get >50% infected"),
                       limits = c(NA, ifelse(pred_var_name == 'CD19_fracCCR5',
                                             30, 60)))+ 
    theme_classic() + labs(x = pred_var_axis_name) + 
    theme(axis.title.y = element_text(color = "red"),
          axis.title.y.right = element_text(color = "blue"))
  
  print(p_prob_infect_v_KO_plot2)
  ggsave(paste0('output/Figures/Fig7F_', Model_name, '_Prob_infect.eps'), p_prob_infect_v_KO_plot2, width = 8,
         height = 8, units = 'cm')
  
  All_plot = plot_grid(Expect_surv_plot_comb2, p_prob_infect_v_KO_plot2,
                       fit_quality_plot, fit_quality_percent, fit_quality_n,
                       ncol =2)
  
  ggsave(paste0('output/Figures/', Model_name, '_all_plots.pdf'), 
         All_plot, width = 16,
         height = 22, units = 'cm')

}

