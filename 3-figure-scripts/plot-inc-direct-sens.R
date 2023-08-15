################################################
# Spillover effects of reactive, focal malaria 
# interventions

# Namibia trial
# Plot Sensitivity analysis
# Exclude 500m boundary for direct effects

################################################ 

rm(list=ls())

source(paste0(here::here(), "/0-config.R"))

# load and process primary tmle results ------------------------------------------------------
estimates <- readRDS(paste0(results_path, "namibia_htmle_inc_all.RDS")) %>%
  mutate(reservoir = factor(reservoir, levels = c("Human", "Mosquito", "Human & mosquito")))

## drop observation periods that are not relevant --------------------------------
drops_human <- which(estimates$reservoir=="Human" & estimates$period=="6 months")
drops_mosq <- which(estimates$reservoir=="Mosquito" & estimates$period=="5 weeks")
drops_hm <- which(estimates$reservoir=="Human & mosquito" & estimates$period=="5 weeks")

estimates <- estimates[-c(drops_human, drops_mosq, drops_hm),] %>% 
  filter(Psi_type=="RR" )

## subset to either individual or cohort level model results selected data-adaptively --------------------------------
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
  )) %>% 
  filter(parameter=="direct") %>% 
  mutate(pop = "Intervention recipients <500m from index case")

# load and process sensitivity analysis results -------------------------------------------
sens_human <- readRDS(paste0(results_path, "namibia_htmle_inc_direct_human_sens_direct.RDS"))
sens_mosq <- readRDS(paste0(results_path, "namibia_htmle_inc_direct_mosq_sens_direct.RDS"))
sens_hm <- readRDS(paste0(results_path, "namibia_htmle_inc_direct_hm_sens_direct.RDS"))

sens <- data.frame(rbind(
  sens_human$res_est$estimates,
  sens_mosq$res_est$estimates,
  sens_hm$res_est$estimates
)) %>% filter(Psi_type=="RR")

sens_lb <- reshape2::melt(sens %>% dplyr::select(Psi_hat, CI_l, 
                                                               CI_l_unadj, 
                                                               reservoir, parameter),
                               id.vars = c("reservoir", "parameter", "Psi_hat")) %>% 
  rename(lower = variable,
         lower_value = value)

sens_ub <- reshape2::melt(sens %>% dplyr::select(Psi_hat, CI_u, 
                                                               CI_u_unadj, 
                                                               reservoir, parameter),
                               id.vars = c("reservoir", "parameter", "Psi_hat")) %>% 
  rename(upper = variable,
         upper_value = value)

sens_plot <- cbind(sens_lb, sens_ub %>% dplyr::select(upper, upper_value)) %>% 
  mutate(analysis = case_when(
    lower=="CI_l" ~ "Adjusted for interference",
    lower=="CI_l_unadj" ~ "Traditional"
  )) %>% 
  dplyr::select(reservoir, parameter, Psi_hat, lower, lower_value, 
       upper, upper_value, analysis) %>% 
  mutate(pop = "All intervention recipients")


inc_data <- bind_rows(estimates_plot, sens_plot) %>% 
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


plot_data <- inc_data %>% 
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
  
  # # manually imputing the upper bound of of the CI for direct effect combined
  # mutate(
  #   ub = ifelse(reservoir_plot=="C) Combined intervention" & parameter == "Direct effect" &
  #                 Analysis == "Standard", 3.5, ub)
  # ) %>%
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
) ), aes(x = reservoir_plot, y = efficacy+100)) +
  geom_hline(yintercept = 100, size=0.3) +
  geom_point(aes(col = pop),size=3, position = position_dodge(width=0.2)) +
  geom_linerange(aes(col = pop, ymin = efficacy_lb+100, ymax = efficacy_ub+100), 
                 position = position_dodge(width=0.2)) +
  scale_y_continuous(trans= "log10",
                     breaks = c(5,15, 25, 40,  65, 100, 200, 300, 400, 500),
                     labels = c(-95, -85, -75, -60, -35, 0, 100, 200, 300, 400),
                     limits = c(10, 500),
                     expand = c(0, 0)
  ) +
  scale_color_manual("",values = c("black","#b8908d")) +
  ylab("% Effectiveness (95% CI)") + 
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11),
        axis.text.x = element_text(size=10),
        strip.background = element_blank()) 

plot

ggsave(plot, filename = paste0(figure_path, "plot-inc-direct-sens.png"),
       width = 6, height=3.5)







