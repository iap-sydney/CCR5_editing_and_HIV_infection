
# -------------------------------------------------------------------------
#'Steffen Docken
#'17-3-23
#'Code to plot VL timecourses

# -------------------------------------------------------------------------


# Max pVL vs.  percent CCR5 KO --------------------------------------------

max_pVL_vs_raw_data_plot_it = ggplot(filter(Animal_wFull_infection_data_df, infect_cens == 1), 
                                aes(y = max_pVL, x = BM_CD19_CCR5_KO_mean, 
                                    color = Arm)) + 
  geom_point()  + 
  geom_point(data = filter(Animal_wFull_infection_data_df, infect_cens == 0),
             aes(y = 0.5, x = BM_CD19_CCR5_KO_mean), shape = 1) + xlim(0, 100) +
  theme_classic() + Surv_arm_palette +
  labs(y= TeX(r"(max pVL ($log_{10}$ copies/ml))"), title='', shape = 'Infection',
       x = expression(paste("% ", italic("CCR5"), "-edited"))) + 
  scale_y_continuous(breaks=c(0.5,1,2, 3, 4), 
                     labels=c("Uninfected", "1", "2", "3", "4"),
                     limits = c(0.5,NA)) +
  theme(axis.title.y = element_text(vjust = -10))

print(max_pVL_vs_raw_data_plot_it)
ggsave('output/Figures/Fig7H_peakVL_v_CCR5KO_CD19.eps', max_pVL_vs_raw_data_plot_it +
         theme(legend.position = 'none'), 
       width = 8, height = 8, units = 'cm')


Peak_Spear_Corr <- cor.test(x = filter(Animal_wFull_infection_data_df, infect_cens == 1)$BM_CD19_CCR5_KO_mean,
                            y = filter(Animal_wFull_infection_data_df, infect_cens == 1)$max_pVL,
                            method = 'spearman')

print(paste0("The Spearman rank correlation between CCR5 KO and max pVL is ",
             Peak_Spear_Corr$p.value))
