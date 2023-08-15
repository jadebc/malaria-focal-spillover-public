# Plot incidence TMLE estimates 
# For all sensitivity analyses 

rm(list=ls())

source(paste0(here::here(), "/0-config.R"))


# load results ------------------------------------------------------------
estimates_primary <- readRDS(paste0(results_path, "namibia_htmle_inc_all.RDS")) %>% 
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito")),
         analysis = "1km")

estimates_2km <- readRDS(paste0(results_path,
                                "namibia_htmle_inc_all_sens_spzone_2km.RDS")) %>% 
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito")),
         analysis = "2km")

estimates_3km <- readRDS(paste0(results_path,
                                "namibia_htmle_inc_all_sens_spzone_3km.RDS")) %>% 
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito")),
         analysis = "3km")


# combine results ------------------------------------------------------------
all_estimates <- bind_rows(
  estimates_primary,
  estimates_2km, estimates_3km
) %>% filter(parameter == "spillover") %>% 
  mutate(reservoir = fct_rev(reservoir),
         analysis = fct_rev(analysis))
  

palette = c("#96d2a4","#4da284","#266b6e")
# palette = c("#c4e6c3","#96d2a4","#6dbc90","#4da284","#36877a","#266b6e","#1d4f60")

# make plot - 5 weeks observation period ------------------------------------------------------------
plot = ggplot(all_estimates %>% filter(Psi_type == "RR" ), 
                   aes(x = reservoir, y = Psi_hat)) + 
  geom_hline(yintercept = 1, col = "#A4A5A5") +
  geom_point(aes(col = analysis),
             position = position_dodge(width=0.5), size = 3) +
  geom_linerange(aes(ymin = CI_l, ymax = CI_u, col = analysis),
                 position = position_dodge(width=0.5)) + 
  scale_y_continuous(trans = "log", limits = c(0.2,2),
  breaks = c(0.25, 0.5, 1,2)) +
  scale_color_manual(values = palette) +
  # scale_shape_manual(values = c(16, 17, 17, 15)) + 
  facet_wrap(~period) + 
  coord_flip() +
  ylab("Cumulative incidence ratio (95% CI)") + 
  theme_bw() + 
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 11),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  ) +   guides(col = guide_legend(reverse = TRUE))


plot

ggsave(plot, filename = paste0(figure_path, "plot-inc-est-tmle-short-spzone.png"),
       width= 6, height=3)




