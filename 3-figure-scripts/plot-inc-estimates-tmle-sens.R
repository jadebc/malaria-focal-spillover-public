################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Plot incidence TMLE estimates 
# For sensitivity analyses 

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
estimates_nooverlap_target <- readRDS(paste0(results_path, 
                                    "namibia_htmle_inc_all_sens_nooverlap_target.RDS")) %>% 
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito")),
         analysis = "No overlap of target areas")

estimates_nooverlap_spill <- readRDS(paste0(results_path, 
                                    "namibia_htmle_inc_all_sens_nooverlap_spill.RDS")) %>% 
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito")),
         analysis = "No overlap of spillover zones")

estimates_obs <- readRDS(paste0(results_path, 
                                  "namibia_htmle_inc_all_sens_obs.RDS")) %>% 
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito")),
         analysis = "Different observation period")


# combine results ------------------------------------------------------------
all_estimates <- bind_rows(
  estimates_primary_sub,
  estimates_nooverlap_target, 
  estimates_nooverlap_spill,
  estimates_obs

) %>% 
  mutate(analysis = ifelse(analysis=="Primary", "Primary analysis", analysis)) %>% 
  mutate(analysis = ifelse(analysis=="Different observation period", "Alternative observation period", analysis)) %>% 
  mutate(analysis = fct_relevel(analysis, "Primary analysis")) %>% 
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

all_estimates <- all_estimates %>% mutate(
 parameter = case_when(parameter == "direct" ~ "Direct\neffect",
                       parameter == "spillover" ~ "Spillover\neffect",
                       parameter == "total" ~ "Total\neffect") 
) %>% 
  mutate(analysis = factor(analysis, levels = c(
    "Primary analysis", "Alternative observation period", 
    # "Spillover zone = 2km", "Spillover zone = 3km",
    "No overlap of spillover zones", "No overlap of target areas"
  )))


palette = c("#000000", "#1D6996","#73AF48",
            "#E17C05")

# make plot ------------------------------------------------------------
plot = ggplot(all_estimates %>% filter(Psi_type == "RR"), 
                   aes(x = parameter, y = Psi_hat)) + 
  geom_hline(yintercept = 1, col = "#A4A5A5") +
  geom_point(aes(col = analysis),
             position = position_dodge(width=0.65), size = 2) +
  geom_linerange(aes(ymin = CI_l, ymax = CI_u, col = analysis),
                 position = position_dodge(width=0.65)) + 
  scale_y_continuous(trans = "log", limits = c(0.05,8),
  breaks = c(0.1, 0.2,0.4, 0.6, 0.8, 1,  2, 4, 8), expand= c(0,0)) +
  scale_color_manual(values = palette) +
  facet_wrap(~comparison) + 
  ylab("Cumulative incidence ratio (95% CI)") + 
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.title = element_blank(),
        legend.box="vertical",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        strip.background = element_blank()
  )  +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) 
plot

ggsave(plot, filename = paste0(figure_path, "plot-inc-est-tmle-sens.png"),
       width= 7, height=4)




