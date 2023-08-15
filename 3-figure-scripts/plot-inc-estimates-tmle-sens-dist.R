################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Plot incidence TMLE estimates 
# For sensitivity analyses for different distance radii

# The primary analyses ran both Qlevel indiv and cohort, so
# subsetting to the level picked by the adaptive algorithm. 
################################################ 
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))
library(tidyverse)

# load primary results ------------------------------------------------------------
estimates_primary <- readRDS(paste0(results_path, "namibia_htmle_inc_all.RDS")) %>% 
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito")),
         analysis = "Primary")

Qlevel = readRDS(paste0(results_path, "namibia_htmle_inc_adaptQ.RDS")) 

estimates_primary_sub <- left_join(estimates_primary, Qlevel, by = c("reservoir", "parameter")) 

drops_short <- which((estimates_primary_sub$Short != estimates_primary_sub$Qlevel) & 
                       estimates_primary_sub$period=="5 weeks")
drops_long <- which((estimates_primary_sub$Long != estimates_primary_sub$Qlevel) & 
                      estimates_primary_sub$period=="6 months")

# dropping indiv total hm even though adaptive Q didn't choose because all other
# long obs analyses chose cohort level
drops_te <- which(estimates_primary_sub$Qlevel=="individual" & 
                    estimates_primary_sub$period=="6 months" &
                      estimates_primary_sub$parameter=="total" & 
                    estimates_primary_sub$reservoir=="Human & mosquito")

estimates_primary_sub <- estimates_primary_sub[-c(drops_short, drops_long, drops_te),]

# load sensitivity results ------------------------------------------------------------
estimates_2km <- readRDS(paste0(results_path,
                                "namibia_htmle_inc_all_sens_spzone_2km.RDS")) %>%
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito")),
         analysis = "Spillover zone = 2km")

estimates_3km <- readRDS(paste0(results_path,
                                "namibia_htmle_inc_all_sens_spzone_3km.RDS")) %>%
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito")),
         analysis = "Spillover zone = 3km")


# combine results ------------------------------------------------------------
all_estimates <- bind_rows(
  estimates_primary_sub,
  estimates_2km, estimates_3km
) %>% 
  filter(parameter=="spillover") %>% 
  mutate(analysis = ifelse(analysis=="Primary", 1, analysis)) %>% 
  mutate(analysis = ifelse(analysis=="Spillover zone = 2km", 2, analysis)) %>% 
  mutate(analysis = ifelse(analysis=="Spillover zone = 3km", 3, analysis)) %>% 
  # mutate(analysis = as.numeric(as.character(analysis))) %>% 
  mutate(comparison = case_when(
    reservoir == "Human" ~ "A) Human intervention\nrfMDA vs. RACD",
    reservoir == "Mosquito" ~ "B) Mosquito intervention\nRAVC vs. No RAVC",
    reservoir == "Human & mosquito" ~ "C) Combined intervention\nrfMDA + RAVC vs. RACD")) %>%
  mutate(comparison = factor(comparison, levels = c("A) Human intervention\nrfMDA vs. RACD", 
                                                    "B) Mosquito intervention\nRAVC vs. No RAVC",
                                                    "C) Combined intervention\nrfMDA + RAVC vs. RACD")))


# subset to observation period that align with arm
drops_human <- which(all_estimates$period=="6 months" & all_estimates$reservoir=="Human")
drops_mosq <- which(all_estimates$period=="5 weeks" & all_estimates$reservoir=="Mosquito")
drops_hm <- which(all_estimates$period=="5 weeks" & all_estimates$reservoir=="Human & mosquito")


all_estimates <- all_estimates[-c(drops_human, drops_mosq, drops_hm),]

# make plot ------------------------------------------------------------
plot = ggplot(all_estimates %>% filter(Psi_type == "RR"), 
                   aes(x = analysis, y = Psi_hat)) + 
  geom_hline(yintercept = 1, col = "#A4A5A5") +
  geom_point(aes(col = analysis),
             position = position_dodge(width=0.65), size = 2) +
  geom_linerange(aes(ymin = CI_l, ymax = CI_u, col = analysis),
                 position = position_dodge(width=0.65)) + 
  scale_y_continuous(trans = "log", limits = c(0.25,4),
  breaks = c(0.4, 0.6, 0.8, 1,  2, 4)) +
  scale_color_manual(values = c("#0e495c", "#3e7c8f", "#8bc7da")) +
  facet_wrap(~comparison) +
  xlab("Radius around index case (km)") + 
  ylab("Cumulative incidence ratio (95% CI)") + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(size=10, face="bold"),
        strip.background = element_blank()
  )
plot

ggsave(plot, filename = paste0(figure_path, "plot-inc-est-tmle-sens-dist.png"),
       width= 6, height=3)


# percent of int recips in target area for figure caption  ------------------------------------------------------------
df_long <- readRDS(namibia_df_long_path)

df_long %>% group_by(intarm) %>% 
  filter(int_recip==1) %>% 
  summarise(percent_intarget = mean(target_area))


