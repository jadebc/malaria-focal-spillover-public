################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Plot of prevalence results for spillover
# Stratified by distance 
################################################ 

rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

# process results --------------------------------------------------
res = readRDS(paste0(results_path, "prevalence-tmle-results-dist.RDS"))

plot_data = bind_rows(
  res$human$sp1km$res_est %>% filter(type=="RR") %>% mutate(reservoir="A) Human intervention\nrfMDA vs. RACD", distance="500m-1km"),
  res$human$sp2km$res_est %>% filter(type=="RR") %>% mutate(reservoir="A) Human intervention\nrfMDA vs. RACD", distance="1-2km"),
  res$human$sp3km$res_est %>% filter(type=="RR") %>% mutate(reservoir="A) Human intervention\nrfMDA vs. RACD", distance="2-3km"),
  
  res$mosq$sp1km$res_est %>% filter(type=="RR") %>% mutate(reservoir="B) Mosquito intervention\nRAVC vs. No RAVC", distance="500m-1km"),
  res$mosq$sp2km$res_est %>% filter(type=="RR") %>% mutate(reservoir="B) Mosquito intervention\nRAVC vs. No RAVC", distance="1-2km"),
  res$mosq$sp3km$res_est %>% filter(type=="RR") %>% mutate(reservoir="B) Mosquito intervention\nRAVC vs. No RAVC", distance="2-3km"),
  
  res$both$sp1km$res_est %>% filter(type=="RR") %>% mutate(reservoir="C) Combined intervention\nrfMDA + RAVC vs. RACD only", distance="500m-1km"),
  res$both$sp2km$res_est %>% filter(type=="RR") %>% mutate(reservoir="C) Combined intervention\nrfMDA + RAVC vs. RACD only", distance="1-2km"),
  res$both$sp3km$res_est %>% filter(type=="RR") %>% mutate(reservoir="C) Combined intervention\nrfMDA + RAVC vs. RACD only", distance="2-3km")
) %>% 
  mutate(reservoir = factor(reservoir, 
                            levels = c("A) Human intervention\nrfMDA vs. RACD", 
                                       "B) Mosquito intervention\nRAVC vs. No RAVC", 
                                       "C) Combined intervention\nrfMDA + RAVC vs. RACD only")),
         distance = factor(distance, levels = c("500m-1km", "1-2km", "2-3km"))) %>% 
  mutate(efficacy = (psi_transformed-1)*100,
         efficacy_lb = (lower_transformed-1)*100,
         efficacy_ub = (upper_transformed-1)*100)

# make plot --------------------------------------------------
plot <- ggplot(plot_data, aes(x = distance, y = efficacy+100)) + 
  geom_hline(yintercept = 100, size=0.3) +
  geom_point(aes(col=distance)) + 
  geom_linerange(aes(col=distance, ymin = efficacy_lb+100, ymax = efficacy_ub+100)) +
  theme_bw() +
  ylab("% Effectiveness (95% CI)") + 
  xlab("Distance to nearest intervention recipient") + 
  facet_grid(~reservoir) +
  scale_color_manual(values = c("#6c94cc", "#204a85", "#042452")) + 
  scale_y_continuous(trans= "log10",
                     breaks = c(1, 5, 15, 25, 40,  65, 100, 200, 300, 500),
                     labels = c(-99, -95, -85, -75, -60, -35, 0, 100, 200, 400),
                     limits = c(1, 600), 
                     expand = c(0, 0)
  ) +
  theme(legend.position = "none",
        strip.text = element_text(size=11, face="bold"),
        strip.background = element_blank()) 
  
ggsave(plot, filename = paste0(figure_path, "plot-prev-est-tmle-dist.png"),
       width=8, height=3)
