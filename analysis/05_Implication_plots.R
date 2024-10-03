
# -------------------------------------------------------------------------
#' Steffen Docken
#' 19-4-23
#' Code to plot expected implications of kCCR5 KO
#' 
# -------------------------------------------------------------------------

human_R_0 = 8

Implication_df = data.frame(CCR5_KO_perc = seq(from = 0, to = 100, by = 2))

Implication_df = Implication_df %>% rowwise() %>% 
  mutate(frac_CCR5 = 1 - CCR5_KO_perc/100,
         Expected_react_delay = (1-exp(-b_h))/(1-exp(-b_h*frac_CCR5)),
         rel_SP = pmax((human_R_0*frac_CCR5 -1)/(human_R_0 -1),0)) %>%
  ungroup()

prop_const = 10

Implication_plot = ggplot(data = Implication_df, aes(x = CCR5_KO_perc)) +
  geom_line(aes(y = Expected_react_delay), color = "red") + 
  geom_line(aes(y = rel_SP*prop_const), color = "blue") +
  scale_y_continuous(name = "Expected time to rebound (weeks)", 
                     sec.axis = sec_axis(~./prop_const, name = "Fold change in set-point"),
                     limits = c(NA, 12)) +
  theme_classic() + labs(x = expression(paste("% ", italic("CCR5"), "-edited"))) + 
  theme(axis.title.y = element_text(color = "red"),
        axis.title.y.right = element_text(color = "blue"))

print(Implication_plot)

ggsave('output/Figures/Fig8B_Implication_Figure.eps', Implication_plot, width = 8,
       height = 8, units = 'cm')
