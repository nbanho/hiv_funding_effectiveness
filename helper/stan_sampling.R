default_sampling <- function(
  mod, # compiled stan model 
  dat_list # data list
  ) {
  init_fun <- function() {
    list(mu_p_a = runif(1, 1, 3),
         sigma_p_a = runif(1, 0.25, 0.75))
  }
  mod$sample(
    data = dat_list,
    seed = 12345,
    parallel_chains = 4,
    threads_per_chain = 2,
    iter_warmup = 1000,
    iter_sampling = 1000,
    max_treedepth = 15,
    init = init_fun
  )
}