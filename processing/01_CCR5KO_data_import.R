
# -------------------------------------------------------------------------
#'Steffen Docken
#'17-3-23
#'This code imports and cleans the data

# -------------------------------------------------------------------------


# load raw data -----------------------------------------------------------

if(!file.exists(data_file)){
  stop(glue("{data_file} does not exist."))
}

Animal_infection_data <- read.csv(data_file)

Animal_infection_data_df <- data.frame(Animal = Animal_infection_data$Hu.HSC.ID,
                                       week_of_infection = Animal_infection_data$Time.of.infection..wks.post.1st.challenge.,
                                       infect_cens = Animal_infection_data$X0...lost.to.follow.up..1...infected, 
                                       #1 = infection detected, 0 = censored at week_of_infection time
                                       Arm = factor(Animal_infection_data$Experimental.CCR5.KO.Cat,
                                                    levels = c("0% CCR5-KO", "25% CCR5-KO", 
                                                               "50% CCR5-KO", "75% CCR5-KO", 
                                                               "100% CCR5-KO")),
                                       Infection_wave = factor(Animal_infection_data$Infection.wave),
                                       BM_CD19_CCR5_KO_min = Animal_infection_data$CD19..cells_Total.ccr5.KO_.Min.editing., 
                                       #minimum estimate % without CCR5
                                       BM_CD19_CCR5_KO_max = Animal_infection_data$CD19..cells_Total.ccr5.KO_.Max.editing.,
                                       #maximum estimate % without CCR5
                                       max_pVL = Animal_infection_data$max_pVL)

Animal_infection_data_df = Animal_infection_data_df %>% 
  mutate(infect_cens = if_else(week_of_infection <= 8, infect_cens, as.integer(0)),
         week_of_infection = if_else(week_of_infection <=8, week_of_infection, as.integer(8)))
#censoring animals detected or censored after week 9 (not including animals 
#detected more than a week after last inoculation)

Animal_infection_data_df = Animal_infection_data_df %>% rowwise() %>%
  mutate(BM_CD19_CCR5_KO_mean = 
           min(max((BM_CD19_CCR5_KO_min + BM_CD19_CCR5_KO_max)/2,0), 100)) %>% 
  ungroup()
         #ensuring that the mean CD19 CCR5 KO is >=0% and <=100%
