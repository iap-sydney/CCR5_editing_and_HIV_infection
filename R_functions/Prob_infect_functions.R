#Steffen Docken
#14-3-23
#Script to contain codes used in CCR5_KO_fitting.R



neg_logLike_infection <- function(params, Infection_df, param.id){

  Infection_df = Infection_df %>% rowwise() %>% mutate(prob_inoc_infect =
    prob_infect(params = params, x = pred_var, param.id = param.id)) %>% ungroup()
  
  Infection_df = Infection_df %>%
    mutate(loglike_animal = log((1-prob_inoc_infect)^(week_of_infection -infect_cens)) +
             log((prob_inoc_infect)^infect_cens)) #calculating log-likelihood of 
  #each animal getting infected at the observed time
  
  logLike_val = max(sum(Infection_df$loglike_animal), -1e8) #summing the animal specific
  #log-likelihoods. If logLikelihood evaluates to -Inf, this also sets it to -1e8
  #to avoid errors.
  
  return(-logLike_val)
  
}

neg_logLike_infection_Arm <- function(params, Infection_df){
  Arm_ind_link <- c("100% CCR5-KO" = 1, "75% CCR5-KO" = 2, "50% CCR5-KO" = 3,
                      "25% CCR5-KO" = 4, "0% CCR5-KO" = 5)
  
  Infection_df =  Infection_df %>% mutate(Arm_ind = Arm_ind_link[Arm])
  
  Infection_df$prob_inoc_infect =
    prob_infect_Arm(params = params, Arm_ind = Infection_df$Arm_ind)
  #browser()
  Infection_df = Infection_df %>%
    mutate(loglike_animal = log((1-prob_inoc_infect)^(week_of_infection -infect_cens)) +
             log((prob_inoc_infect)^infect_cens)) #calculating log-likelihood of 
  #each animal getting infected at the observed time
  #browser()
  logLike_val = max(sum(Infection_df$loglike_animal), -1e8) #summing the animal specific
  #log-likelihoods. If logLikelihood evaluates to -Inf, this also sets it to -1e8
  #to avoid errors.
  if ((max(params)>1)|(min(params)<0)){
    logLike_val=-1e8 #ensuring probability of infection is between 0 and 1
  }
  
  return(-logLike_val)
  
}

fit_infection_prob_Arm <- function(Infection_df){
  
  Infection_df =  Infection_df %>% mutate(Arm_ind = Arm_ind_link[Arm])
  
  params = rep(0,5) #will contain probability of infection per inoculation for 
  #100% KO, 75%, etc.
  
  for (ii in 1:5) {
    Infection_df_it = filter(Infection_df, Arm_ind == ii)
    
    params[ii] = ifelse(sum(Infection_df_it$infect_cens) == 0,
                        0, #if no animals were infected
      1/(1 + sum(Infection_df_it$week_of_infection - 
                              Infection_df_it$infect_cens)/
                      sum(Infection_df_it$infect_cens)))
  }
  
  return(params)
  
}

fit_infection_prob_all <- function(Infection_df){
    
  avg_infect_risk = ifelse(sum(Infection_df$infect_cens) == 0,
                      0, #if no animals were infected
                      1/(1 + sum(Infection_df$week_of_infection - 
                                   Infection_df$infect_cens)/
                           sum(Infection_df$infect_cens)))
  
  params = rep(avg_infect_risk, 5) #repeating average infection risk 5 times to 
  # output optimization in consistent form.
  
  return(params)
  
}


prob_infect_week <- function(week_of_infection, prob_inoc_infect, infect_cens){
  #function to calculate the probability of getting infected on a specific week
  #(i.e., after 'week_of_infection' inoculations)
  #prob_inoc_infect is the probability of a single inoculation causing an infection
  #infect_cens: 1 infection detected, 0 censoring occured at current week
  
  
  Prob_infect_week = (1-prob_inoc_infect)^(week_of_infection -infect_cens)*
    (prob_inoc_infect)^infect_cens
  
  
  return(Prob_infect_week)
  
}

prob_infect <-function(params, x, param.id){
  #x is either the concentration of target cells or the fraction of
  #cell with CCR5
  theta = rep(0, 5) #will contain the parameters for this round (including those
  #set to 0). In order, the parameters contained in theta are the coefficients of
  #x^2, x, constant term, 1/x, ln(x)
  
  theta[param.id] = params #updating the relevant elements of theta for the
  #parameters used here.
  
  if(theta[4]==0){
    exp_term = sum(c(theta[1]*x^2, theta[2]*x, theta[3],
                     theta[5]*log(x)), na.rm = TRUE) # calculating the general 
    #exponential term ax^2 + bx + c + f*ln(x). NaNs are ignored to evaluate
    #0*log(0) to 0
    
    exp_term = pmax(0, exp_term) #if the exponential term is <0, that would imply
    #exp(-exp_term) > 1, and prob_infect is <0, which is nonsensical.
  }else if(theta[4]*x <= 1){
    exp_term = 0 #This means R_0 <= 1, so sustained infection is not possible. 
    #Therefore, prob_infect will be 0
  }else{
    exp_term = sum(c(theta[1]*x^2, theta[2]*x, theta[3]))*
      (1-1/(theta[4]*x))# calculating the general 
    #exponential term (ax^2 + bx + c)(Rx-1)/Rx.
    
    exp_term = pmax(0, exp_term) #if the exponential term is <0, that would imply
    #exp(-exp_term) > 1, and prob_infect is <0, which is nonsensical.
  }
  
  prob_infect = 1-exp(-exp_term)
  
  return(prob_infect)
}

prob_infect_Arm <-function(params, Arm_ind){
  #params is the list of probability of getting infected after inoculation for
  # each Arm (100% KO, 75% KO, ...) and Arm is the index for which Arm the 
  # animal was part of (1=100% KO, 2 = 75% KO, ....)
  
  prob_infect_it = params[Arm_ind] #probability of getting infected for this
  #animal
  
  return(prob_infect_it)
}
