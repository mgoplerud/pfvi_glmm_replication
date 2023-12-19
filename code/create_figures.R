rm(list=ls())
library(tidyverse)

################################################################################
# Create Figures 1 and 2 describing the results from the simulation study
################################################################################

simulation_output <- readRDS('output/all_simulations.RDS')

limited_plot <- c('Collapse [0]' = 'Fully Factorized',
  'Collapse [FE]' = 'Partially Factorized',
  'CAVI Unfactorized' = 'Unfactorized')

recode_model <- function(x, lvl){
  x <- x %>% filter(model %in% names(lvl)) %>% 
    mutate(model = recode_factor(model, !!!as.list(lvl)))
  x <- x %>% mutate(
    family = recode_factor(family,
                           'linear' = 'Gaussian', 'binomial' = 'Binomial')
  )
  return(x)
}

summary_uqf <- simulation_output$uqf %>% 
  group_by(model, method, levels, family) %>%
  summarize(n_obs = n(), se_uqf = sd(uqf_split_first)/sqrt(n()), 
            acc_RE = mean(acc_RE), acc_FE = mean(acc_FE),
            uqf = mean(uqf_split_sample),
            min_ratio = mean(ratio_split)) %>% ungroup %>%
  mutate(levels = 2^levels)
summary_uqf <- recode_model(summary_uqf, limited_plot)

g_uqf_sim <- ggplot(summary_uqf, 
    aes(x=levels * 2,y=uqf, col = model, 
        group = interaction(model, method))) +
  geom_point() +
  geom_line() +
  facet_grid(. ~ family) + 
  xlab('Number of Parameters (Log-Scale)') + 
  ylab('UQF (Log-Scale)')  +
  scale_x_log10() + scale_y_log10() + theme_bw() +
  theme(legend.position = 'bottom') +
  labs(col = 'Method: ')

ggsave(plot = g_uqf_sim, 
       filename = 'figures/uqf_simulation.pdf',
       width = 8.5, height = 8.5/2)

summary_time <- simulation_output$time %>%
  group_by(levels, family, model) %>%
  summarize(time = mean(time_CAVI_only),
            time_per_iter = mean(time_CAVI_only/iterations)) %>%
  ungroup

summary_time <- recode_model(summary_time, limited_plot)


g_sim_time <- ggplot(summary_time, aes(x=2^levels * 2,
    y =time_per_iter * 60, group = model, col = model)) + 
  geom_point() + geom_line() +
  facet_wrap(~family) + scale_y_log10() +
  scale_x_log10() + theme_bw() +
  theme_bw() + theme(legend.position = 'bottom') +
  xlab('Number of Parameters (Log-Scale)') +
  ylab('Time per Iteration\nin Seconds (Log-Scale)') +
  labs(col = 'Method: ')

ggsave(plot = g_sim_time, 
       filename = 'figures/time_simulation.pdf', 
       width = 8.5, height = 8.5/2)

################################################################################
# Create Table 1 describing the models for the Ghitza-Gelman analysis
################################################################################

pf_models <- dir("output/output_vglmer", pattern='collapse_main.*binomial')

if (length(pf_models) > 0){
  # only run if the models actually exist in output/output_vglmer
  # otherwise, use existing version  
  model_sizes <- lapply(pf_models, FUN=function(i){
    model_i <- readRDS(paste0('output/output_vglmer/', i))
    size_C <- length(unlist(model_i$collapsed$index))
    size_M <- length(unlist(model_i$marginal$index))
    total_size <- nrow(vglmer::format_vglmer(model_i))
    return(data.frame(C = size_C, M = size_M, total = total_size,
                      file = i, stringsAsFactors = FALSE))
  }) %>% bind_rows()
  
  
  model_sizes <- model_sizes %>% 
    mutate(year = str_extract(file, pattern='(?<=y_)[0-9]+')) %>%
    mutate(model = str_extract(file, pattern='(?<=m_)[0-9]+'))
  
  model_sizes <- model_sizes %>% dplyr::select(model, C, M, total) %>% unique
  
  gg_table <- readxl::read_excel('data/gg_table_description.xlsx')
  gg_table <- left_join(
    gg_table %>% mutate(model = as.character(model)), 
    model_sizes, by = c('model'))
  
  gg_table %>% 
    transmute(name, description, C, M) %>% 
    apply(., MARGIN = 1, paste, collapse=' & ') %>%
    paste0(., '\\\\') %>%
    writeLines('figures/gg_table.tex')
  
}

################################################################################
# Create Figure 3 show the UQF for the Ghitza-Gelman models
################################################################################

limited_plot <- c('vglmer_strong' = 'Fully Factorized',
                  'vglmer_collapse_FE' = 'Partially Factorized (Fixed)',
                  'vglmer_collapse_main' = 'Partially Factorized (Fixed + Main)',
                  'vglmer_collapse_all' = 'Unfactorized')

uqf_data <- lapply(dir('output/output_uqf', full.names = TRUE), readRDS)

all_uqf <- bind_rows(lapply(uqf_data, `[[`, 'uqf'))

all_uqf <- all_uqf %>% filter(factor_method %in% names(limited_plot)) %>%
  mutate(method = recode_factor(factor_method, !!! limited_plot))

avg_uqf <- all_uqf %>%
  group_by(model, outcome, method) %>%
  summarize(across(matches('uqf'), mean))

g_uqf <- ggplot(avg_uqf,
      aes(x=factor(model),y=uqf_split_sample, col=method,group=method)) + 
  geom_point() + geom_line() +
  scale_y_log10(breaks = 10^(-4:0), labels = scales::comma, limits = c(10^-4, 1)) + 
  ylab('UQF (Log-Scale)') + xlab('Model') +
  theme_bw() + theme(legend.position = 'bottom') +
  labs(col = 'Method: ')

ggsave(plot = g_uqf, 'figures/gg_uqf.pdf', width = 8.5,
       height = 8.5/2)

################################################################################
# Create Figure 4 calculating the accuracy for the Ghitza-Gelman models
################################################################################

MRP_data <- lapply(dir('output/output_mrp/', full.names = TRUE), readRDS) %>%
  bind_rows()

MRP_data <- MRP_data %>% filter(method %in% names(limited_plot)) %>%
  mutate(method = recode_factor(method, !!! limited_plot))

g_MRP <- ggplot(MRP_data %>% filter(combo == 'stt')) +
  geom_boxplot(aes(x=factor(model),
                   y=accuracy, 
                   group = interaction(model, method),
                   col = method)) +
  xlab('Model') + ylab('Accuracy (Log-Scale)') + theme_bw() +
  labs(col = 'Method: ') +
  scale_y_log10(limits = range(MRP_data$accuracy))
  
g_lp <- ggplot(MRP_data %>% filter(combo == 'stt_eth_inc_age')) +
  geom_boxplot(aes(x=factor(model),
                   y=accuracy, 
                   group = interaction(model, method),
                   col = method)) +
  xlab('Model') + ylab('Accuracy (Log-Scale)') + theme_bw() +
  labs(col = 'Method: ') +
  scale_y_log10(limits = range(MRP_data$accuracy))
  
g_both <- ggpubr::ggarrange(
  g_MRP + ggtitle('(a) State-Level MRP'),
  g_lp + ggtitle('(b) Linear Predictor'),
  common.legend = TRUE, legend = 'bottom'
)

ggsave(plot = g_both, 
       filename = 'figures/gg_acc.pdf', 
       width = 8.5, height = 8.5/2)

MRP_data %>% group_by(method, combo, model) %>%
  summarize(accuracy = mean(accuracy)) %>% ungroup %>%
  pivot_wider(id_cols = c(combo, method), names_from = model,
              values_from = accuracy) %>%
  arrange(combo)

# Additional analysis on computational time

gg_timing <- readRDS('output/vglmer_fit.RDS')
gg_timing <- lapply(gg_timing, `[[`, "timing") %>%
  bind_rows(.id = 'file') %>% 
  separate(file, into = c('year', 'model', 'outcome'), sep = '-') %>%
  filter(outcome == 'binomial')
gg_timing <- gg_timing %>% group_by(model, factor_method) %>%
  summarize(n=n(), iterations = mean(iterations), elbo = min(change_ELBO),
            time = mean(time))
gg_timing <- gg_timing %>% filter(!(factor_method %in% c('glmer', 'hmc')))

limited_plot <- c('collapse_none' = 'Fully Factorized',
                  'collapse_FE' = 'Partially Factorized (Fixed)',
                  'collapse_main' = 'Partially Factorized (Fixed + Main)',
                  'collapse_all' = 'Unfactorized')

gg_timing <- gg_timing %>%
  filter(factor_method %in% names(limited_plot)) %>%
  mutate(factor_method = recode_factor(factor_method, !!! limited_plot))

ggplot(gg_timing,
       aes(x=factor(model),y=time, col = factor_method, group = factor_method)) + 
  geom_point() + geom_line() +
  scale_y_log10()  + ylab('Time (Log-Scale)') + 
    xlab('Model') +
    theme(legend.position = 'bottom') +
    labs(col = 'Method: ')
  

ggplot(gg_timing %>% group_by(model) %>% mutate(ratio = time/time[factor_method == 'Unfactorized']),
       aes(x=factor(model),y=ratio, col = factor_method, group = factor_method)) + 
  geom_point() + geom_line() +
  scale_y_log10()  + ylab('Time Ratio vs. Unfactorized (Log-Scale)') + 
  xlab('Model') +
  theme(legend.position = 'bottom') +
  labs(col = 'Method: ')

print(gg_timing %>% filter(model == 9))