run_sens <- function() {
  # which sensitivity analysis
  args <- commandArgs(trailingOnly = TRUE)
  do_mod <- ifelse(grepl("all|mod", args[1]), T, F)
  do_par <- ifelse(grepl("all|par", args[1]), T, F)
  do_tar <- ifelse(grepl("all|tar", args[1]), T, F)
  do_loo <- ifelse(grepl("all|loo", args[1]), T, F)
  
  # R libraries
  library(tidyverse)
  library(cmdstanr)
  
  # project-specific helper functions
  source("helper/create_stan_data.R")
  source("helper/stan_sampling.R")

  # directory to save files
  save_dir <- 'fitted-models/sens/'
  
  # remove old files
  file.remove(list.files(save_dir, full.names = T))
  
  # load data
  df <- read.csv("data/preprocessed/data.csv")
  
  # default model
  def_model <- "models/active_funding_notr_me.stan"
  
  # compile stan model
  comp_def_model <- cmdstan_model(def_model, cpp_options = list(stan_threads = T)) 
  
  # model comparison
  if (do_mod) {
    # baseline models
    base <- cmdstan_model("models/base.stan", cpp_options = list(stan_threads = T)) 
    base_controls <- cmdstan_model("models/base_controls.stan", cpp_options = list(stan_threads = T)) 
    trend <- cmdstan_model("models/trend.stan", cpp_options = list(stan_threads = T)) 
    main_at <- cmdstan_model("models/active_funding_at_me.stan", cpp_options = list(stan_threads = T)) 
    main_trend <- cmdstan_model("models/active_funding_me.stan", cpp_options = list(stan_threads = T))
    
    # create stan data
    sdat <- create_stan_data(df) 
    sdat_aid_only <- sdat
    sdat_aid_only$D <- 1
    sdat_aid_only$funding <- as.matrix(sdat$funding[,1])
    sdat_dom_only <- sdat
    sdat_dom_only$D <- 1
    sdat_dom_only$funding <- as.matrix(sdat$funding[,2])
    sdat_dom_only$p_a_hyper_diff <- c(-0.5, 0.25)
    
    # model and save
    fit <- default_sampling(base, sdat)
    fit$save_output_files(save_dir, basename = "base")
    fit <- default_sampling(base_controls, sdat)
    fit$save_output_files(save_dir, basename = "base_controls")
    fit <- default_sampling(trend, sdat)
    fit$save_output_files(save_dir, basename = "trend")
    fit <- default_sampling(main_at, sdat)
    fit$save_output_files(save_dir, basename = "main_at")
    fit <- default_sampling(main_trend, sdat)
    fit$save_output_files(save_dir, basename = "main_trend")
    fit <- default_sampling(comp_def_model, sdat_aid_only)
    fit$save_output_files(save_dir, basename = "main_aid_only")
    fit <- default_sampling(comp_def_model, sdat_dom_only)
    fit$save_output_files(save_dir, basename = "main_dom_only")
  }
  
  # parameter sensitivity analysis
  if (do_par) {
    
    # pa_hyper
    mod_starts <- c(6,5,5,7,7)
    grid <- list(
      default = c(2.0, 0.25, 0.5, 0.1),
      very_early = c(1.0, 0.25, 0.5, 0.1),
      early = c(1.5, 0.25, 0.5, 0.1),
      late = c(2.5, 0.25, 0.5, 0.1),
      very_late = c(3.0, 0.25, 0.5, 0.1)
    )
    
    # loop over model starts
    for (g in 1:length(grid)) {
      
      # create stan data
      sdat <- create_stan_data(df, ms = mod_starts[g], pa_hp = grid[[g]], pa_hp_diff = c(0.5, 0.25))
      
      # model
      fit <- default_sampling(comp_def_model, sdat)
      
      # save model
      fit$save_object(file = paste0(save_dir, "par-", names(grid)[g], ".rds"))
    }
  }
  
  # target sensitivity analysis
  if (do_tar) {
    # other than main target
    targets <- c("hiv_new_all", "hiv_new_adult")

    # loop over other targets
    for (t in targets) {

      # create stan data
      sdat <- create_stan_data(df, target = t)

      # model
      fit <- default_sampling(comp_def_model, sdat)

      # save model
      fit$save_object(file = paste0(save_dir, "target-", t, ".rds"))
    }
  }
  
  
  # leave-one-country-out sensitivity analysis
  if (do_loo) {

    # all recipients
    recipients <- unique(df$rcpt)
    
    # loop over recipients
    for (loo in recipients) {
      
      # filter recipient
      df_r <- dplyr::filter(df, rcpt != loo)
      
      # create stan data
      sdat <- create_stan_data(df_r)
      
      # model
      fit <- default_sampling(comp_def_model, sdat)
      
      # save model
      fit$save_object(file = paste0(save_dir, "loo-", loo, ".rds"))
    }
  }
}

# execute
run_sens()