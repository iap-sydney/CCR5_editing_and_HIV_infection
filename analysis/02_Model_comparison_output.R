
# -------------------------------------------------------------------------
#' Steffen Docken
#' 3-4-23
#' Code to visualize delta AICs of model fits
#' 
# -------------------------------------------------------------------------

Model_fits <- read.csv('output/Survival_fit_outputs.csv')

AIC_1 = ggplot(data = Model_fits, aes(x = interaction(Pred_var, Model_ID),
                                      y = pmin(delta_AIC, 30))) +
  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(),
                                      axis.title.x = element_blank()) + 
  ylim(0, 30) + labs(y = "delta AIC", x = "")

print(AIC_1)

AIC_2 = ggplot(data = filter(Model_fits, delta_AIC <=10), aes(x = interaction(Pred_var, Model_ID),
                                      y = delta_AIC)) +
  geom_bar(stat = "identity") +  ylim(0, 10) + labs(y = "delta AIC", x = "") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(AIC_2)

full_AIC_plot = plot_grid(AIC_1, AIC_2, ncol = 1, rel_heights = c(1,3))

print(full_AIC_plot)

ggsave('output/Figures/AIC_comparison.pdf', full_AIC_plot)

