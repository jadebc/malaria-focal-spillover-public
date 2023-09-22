# Plot incidence TMLE estimates 
# For all sensitivity analyses 
library(forcats)
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))


# load results ------------------------------------------------------------
all_estimates <- readRDS(paste0(results_path, "namibia_htmle_inc_all_em.RDS")) 

estimates <- all_estimates %>% 
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito")),
         analysis = "Primary") %>% 
  mutate(comparison = case_when(
    reservoir == "Human" ~ "Human intervention\nrfMDA vs. RACD",
    reservoir == "Mosquito" ~ "Mosquito intervention\nRAVC vs. No RAVC",
    reservoir == "Human & mosquito" ~ "Combined intervention\nrfMDA + RAVC vs. RACD only")) %>%
  mutate(comparison = factor(comparison, levels = c("Human intervention\nrfMDA vs. RACD", 
                                                    "Mosquito intervention\nRAVC vs. No RAVC",
                                                    "Combined intervention\nrfMDA + RAVC vs. RACD only"))) %>% 
  mutate(emlevel = as.factor(case_when(
    modifier_level ==1 & modifier!="sex" ~ "Above median",
    modifier_level ==0 & modifier!="sex" ~ "Below median",
    modifier_level ==1 & modifier=="sex" ~ "Male",
    modifier_level ==0 & modifier=="sex" ~ "Female"
  ))) %>% 
  mutate(emvar = as.factor(case_when(
    modifier == "tx_cov_cohort_abovemed" ~ "Cohort-level treatment coverage",
    modifier == "surface_temp_abovemed" ~ "Surface temperature",
    modifier == "pre_spray_cover_abovemed" ~ "Pre-trial IRS coverage",
    modifier == "pre_rainfall_abovemed" ~ "Rainfall",
    modifier == "pre_incidence_abovemed" ~ "Pre-trial incidence",
    modifier == "pre_evi_abovemed" ~ "Enhanced vegetative index",
    modifier == "ea_elevation_abovemed" ~ "Elevation",
    modifier == "sex" ~ "Gender",
    modifier == "min_dist_hf_abovemed" ~ "Distance to nearest health facility"
  )))%>% 
  mutate(emvar_cat = as.factor(case_when(
    emvar == "Cohort-level treatment coverage" ~ "Trial factors",
    emvar == "Elevation" ~ "Environmental factors",
    emvar == "Enhanced vegetative index" ~ "Environmental factors",
    emvar == "Pre-trial incidence" ~ "Epidemiologic factors",
    emvar == "Pre-trial IRS coverage" ~ "Epidemiologic factors",
    emvar == "Rainfall" ~ "Environmental factors",
    emvar == "Surface temperature" ~ "Environmental factors",
    emvar == "Gender" ~ "Epidemiologic factors",
    emvar == "Distance to nearest health facility" ~ "Epidemiologic factors"
  ))) %>% 
  mutate(emvar = factor(emvar, levels = c(
    "Cohort-level treatment coverage",
    "Elevation",
    "Enhanced vegetative index",
    "Rainfall",
    "Surface temperature" ,
    "Pre-trial IRS coverage",
    "Pre-trial incidence",
    "Distance to nearest health facility",
    "Gender"
  ))) 


# plotting function ------------------------------------------------------------
make_plot <- function(data, parameter_name){
  data = data %>% 
    filter(parameter == parameter_name)
  
  if(parameter_name=="total") param = "Total effect"
  if(parameter_name=="spillover") param = "Spillover effect"
  if(parameter_name=="direct") param = "Direct effect"
  
  # truncate upper bound at 16
  data$CI_u = ifelse(data$CI_u>16, 15.9, data$CI_u)
  
  ggplot(data %>% filter(Psi_type == "RR"), aes(x = emvar, y = Psi_hat)) + 
  geom_hline(yintercept = 1, col = "#A4A5A5") +
  geom_point(aes(col = emvar_cat, shape = emlevel),
             position = position_dodge(width=0.4), size = 3) +
  geom_linerange(aes(ymin = CI_l, ymax = CI_u, col = emvar_cat, 
                     shape = emlevel),
                 position = position_dodge(width=0.4)) + 
    scale_y_continuous(trans = "log10",
                       breaks = c(0.02, 0.1,0.3, 0.6, 1,  2, 4, 8, 16),
                       labels = c("0.02", "0.1","0.3","0.6", "1", "2","4","8","16"),
                       limits = c(-0.1, 16)) +
    scale_color_manual(values = c("#119E57","#DE9D03","#128BCB")) +
  scale_shape_manual(values = c(19,1, 15, 17)) +
  coord_flip() + 
  facet_wrap(~comparison) +
  ylab("Cumulative incidence ratio (95% CI)") + 
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 11),
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(size=11)
  ) +   guides(col = guide_legend(reverse = TRUE)) 
}

# make plots ------------------------------------------------------------
plot_spillover <- make_plot(data = estimates,
                                  parameter_name = "spillover")
plot_spillover

# save plots ------------------------------------------------------------
ggsave(plot_spillover, filename = paste0(figure_path, "plot-inc-est-tmle-em.png"),
       width= 11, height=6.5)



