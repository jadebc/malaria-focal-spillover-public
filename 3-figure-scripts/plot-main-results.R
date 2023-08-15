################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Plot with key results
################################################ 

rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

# load and process tmle results ------------------------------------------------------
estimates <- readRDS(paste0(results_path, "namibia_htmle_inc_all.RDS")) %>%
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito")))


# drop observation periods that are not relevant --------------------------------
drops_human <- which(estimates$reservoir=="Human" & estimates$period=="6 months")
drops_mosq <- which(estimates$reservoir=="Mosquito" & estimates$period=="5 weeks")
drops_hm <- which(estimates$reservoir=="Human & mosquito" & estimates$period=="5 weeks")

estimates <- estimates[-c(drops_human, drops_mosq, drops_hm),] %>% 
  filter(Psi_type=="RR" )

# subset to either individual or cohort level model results selected data-adaptively --------------------------------
adaptQ <- readRDS(paste0(results_path, "namibia_htmle_inc_adaptQ.RDS")) %>% 
  mutate(selection = case_when(
    reservoir=="Human" ~ Short,
    reservoir=="Mosquito" ~ Long,
    reservoir=="Human & mosquito" ~ Long
  )) %>% dplyr::select(reservoir, parameter, selection)


estimates_all <- left_join(estimates, adaptQ, by = c("parameter", "reservoir")) %>% 
  filter(Qlevel==selection)

estimates_lb <- reshape2::melt(estimates_all %>% dplyr::select(Psi_hat, CI_l, 
                                                           CI_l_unadj, 
                                                           reservoir, parameter),
                               id.vars = c("reservoir", "parameter", "Psi_hat")) %>% 
  rename(lower = variable,
         lower_value = value)

estimates_ub <- reshape2::melt(estimates_all %>% dplyr::select(Psi_hat, CI_u, 
                                                           CI_u_unadj, 
                                                           reservoir, parameter),
                               id.vars = c("reservoir", "parameter", "Psi_hat")) %>% 
  rename(upper = variable,
         upper_value = value)

estimates_plot <- cbind(estimates_lb, 
                        estimates_ub %>% dplyr::select(upper, upper_value)) %>% 
  mutate(analysis = case_when(
    lower=="CI_l" ~ "Adjusted for interference",
    lower=="CI_l_unadj" ~ "Traditional"
  ))

inc_data <- estimates_plot %>% 
  dplyr::select(-c(lower, upper)) %>% 
  mutate(
    PR = Psi_hat,
    lb = lower_value, 
    ub = upper_value
  ) %>% 
  mutate(Analysis = ifelse(analysis == "Traditional", "Standard", "Interference")) %>% 
  mutate(parameter = case_when(
    parameter == "direct" ~ "Direct effect",
    parameter == "spillover" ~ "Spillover effect",
    parameter == "total" ~ "Total effect"
  ))

prev_data <- readRDS(paste0(namibia_process_path, "plot-data-prev-primary-qpcr.RDS")) %>% 
  dplyr::select(-c(result_PR, result_PR_adj)) %>% 
  mutate(Analysis = "Prevalence") %>% 
  rename(PR = PR_unadj, 
         lb = lb_unadj, 
         ub = ub_unadj)

sero_data <- readRDS(paste0(namibia_process_path, "plot-data-prev-primary-sero.RDS")) %>% 
  dplyr::select(-c(result_PR, result_PR_adj)) %>% 
  mutate(Analysis = "Seroprevalence") %>% 
  rename(PR = PR_unadj, 
         lb = lb_unadj, 
         ub = ub_unadj)


all_data <- bind_rows(inc_data, prev_data, sero_data) %>% 
  dplyr::select(reservoir, parameter, Analysis, PR, lb, ub)
  

plot_data <- bind_rows(all_data) %>% 
  mutate(Analysis = factor(Analysis, levels = c(
                                                "Standard", 
                                                "Interference",
                                                "Prevalence",
                                                "Seroprevalence"))) 

plot_data <- plot_data %>% 
  mutate(reservoir_plot = case_when(
    reservoir == "Human" ~ "A) Human intervention", 
    reservoir == "Mosquito" ~ "B) Mosquito intervention", 
    reservoir == "Human & mosquito" ~ "C) Combined intervention"
  )) %>% 
  mutate(reservoir_plot = factor(reservoir_plot, levels = c(
    "A) Human intervention", 
    "B) Mosquito intervention", 
    "C) Combined intervention"
  ))) %>% 

  # manually imputing the upper bound of of the CI for direct effect combined
  mutate(
    ub = ifelse(reservoir_plot=="C) Combined intervention" & parameter == "Direct effect" &
                  Analysis == "Standard", 3.5, ub)
  ) %>%
  mutate(efficacy = (PR-1)*100,
         efficacy_lb = (lb-1)*100,
         efficacy_ub = (ub-1)*100)


plot_data = plot_data %>% mutate(parameter_f = case_when(
  parameter =="Direct effect" ~ "Direct\neffect",
  parameter =="Spillover effect" ~ "Spillover\neffect",
  parameter =="Total effect" ~ "Total\neffect"
))

colors = c("#1d7e96", "#789c5b", "#75634f")

plot <- ggplot(plot_data %>% filter(Analysis %in% c(
  "Standard"
) ), aes(x = parameter_f, y = efficacy+100)) +
  geom_hline(yintercept = 100, size=0.3) +
  geom_point(aes(col = parameter),size=2) +
  geom_linerange(aes(col = parameter, ymin = efficacy_lb+100, ymax = efficacy_ub+100)) +
  scale_y_continuous(trans= "log10",
                     breaks = c(15, 25, 40,  65, 100, 200, 300),
                     labels = c(-85, -75, -60, -35, 0, 100, 200),
                     limits = c(10, 350), 
                     expand = c(0, 0)
                     ) +
  scale_color_manual(values = colors) +
  ylab("% Effectiveness (95% CI)") + 
  theme_bw() + 
  facet_grid(~reservoir_plot) +
  theme(legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=9),
    axis.text.x = element_text(size=8),
    strip.text = element_text(size=8, face="bold"),
    strip.background = element_blank()) 
  
plot

ggsave(plot, filename = paste0(figure_path, "plot-inc-results.png"),
       width = 5, height=2)





plot_prev <- ggplot(plot_data %>% filter(Analysis %in% c(
  "Prevalence",
  "Seroprevalence"
) ), aes(x = parameter_f, y = efficacy+100)) +
  geom_hline(yintercept = 100, size=0.3) +
  geom_point(aes(col = parameter_f, shape = Analysis), position = position_dodge(width=0.5),
             size=2) +
  geom_linerange(aes(ymin = efficacy_lb+100, ymax = efficacy_ub+100, col = parameter_f, 
                     shape = Analysis), 
                 position = position_dodge(width=0.5)) +
  scale_y_continuous(trans= "log10", 
                     breaks = c(5, 15, 25, 40,  65, 100, 200),
                     labels = c(-95, -85, -75, -60, -35, 0, 100),
                     limits = c(4, 200)) +
  scale_color_manual(values = colors, guide=F) +
  scale_shape_manual(values = c(19, 3)) +
  ylab("% Effectiveness (95% CI)") + 
  theme_bw() + 
  facet_grid(~reservoir_plot) +
  theme(
    legend.position = c(0.13,0.15),
    legend.title = element_blank(),
    legend.spacing.y = unit(0.0, 'cm'), # space between rows
    legend.key.size = unit(0.4, "cm"),
    legend.box.background = element_rect(color = "black"),
    legend.text = element_text(size=8),
    legend.margin = margin(t = 0.0, r=1, l = 0.0, b=0),
    axis.title.x = element_blank(),
    strip.text = element_text(size=8, face="bold"),
    axis.title.y = element_text(size=9),
    axis.text.x = element_text(size=7),
    strip.background = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE))

plot_prev


ggsave(plot_prev, filename = paste0(figure_path, "plot-prev-results.png"),
       width = 5, height=2)

