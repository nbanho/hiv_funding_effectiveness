main_model <- function() {
  # arguments
  args <- commandArgs(trailingOnly = TRUE)
  # model file: default model is active_funding_notr_me
  mfile <- paste0("models/", args[1], ".stan")
  # target/outcome: default is hiv_new_adult
  target <- args[2]
  # model start: default is 6 (i.e. five years + year 0)
  model_start <- as.numeric(args[3])
  
  # R libraries
  library(tidyverse)
  library(cmdstanr)
  
  # project-specific helper functions
  source("helper/create_stan_data.R")
  source("helper/stan_sampling.R")

  # directory to save files
  save_dir <- 'fitted-models/main/'
  
  # remove old files
  file.remove(list.files(save_dir, full.names = T))
  
  # load data
  df <- read.csv("data/preprocessed/data.csv")
  
  # compile stan model
  comp_model <- cmdstan_model(mfile, cpp_options = list(stan_threads = T)) 
  
  # create stan data
  sdat <- create_stan_data(df, target = target, ms = model_start) 
  
  # fit model
  fit <- default_sampling(comp_model, sdat)
  
  # save model output
  fit$save_output_files(save_dir, basename = "main")
}


main_model()