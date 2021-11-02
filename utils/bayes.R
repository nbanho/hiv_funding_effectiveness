require(posterior)
require(cmdstanr)
require(tidybayes)
require(bayesplot)
require(LaplacesDemon)
require(loo)


# Gather draws for posterior draws from read_cmdstan_csv
gather_draws_csv <- function(model, par, mi, ma) {
  model$post_warmup_draws[,,paste0(par, "[", mi:ma, "]")]%>% 
    as_draws_df() %>%
    data.frame() %>%
    reshape2::melt(id.vars = c(".chain", ".iteration", ".draw")) %>%
    mutate(variable = stringi::stri_extract(variable, regex = "\\d+")) 
}

# To support tidybayes gather_draws for cmdstanr
tidy_draws.CmdStanMCMC <- function(model, ...) { return(as_draws_df(model$draws())) }

# LOO for ragged matrices (padded matrices that are originally of unequal length0)
loo.CmdStanCSV <- function(model, N) {
  ll <- paste0("log_lik[", 1:N, "]")
  loo::loo(x = model$post_warmup_draws[,,ll], r_eff = relative_eff(model$post_warmup_draws[,,ll]))
}

# plot influential observations
plot_influential.ragged_mat <- function(model, N, x_axis, facet_g, file, w = 30, h = 30) {
  loos <- loo.CmdStanCSV(model, N)
  
  loo_dat <- data.frame(x_axis =  x_axis,
                        facet_g = facet_g,
                        k = loos$diagnostics$pareto_k) 
  
  loo_pl <- loo_dat %>%
    ggplot(aes(x = x_axis, y = k)) +
    geom_point(shape = 3) +
    facet_wrap(~ facet_g, ncol = floor(sqrt(length(unique(facet_g)))), scales = "free_x") +
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black") +
    geom_hline(aes(yintercept = 0.5), linetype = "dashed", color = "red") +
    geom_hline(aes(yintercept = 0.7), linetype = "dashed", color = "red") +
    scale_y_continuous(limits = c(min(loo_dat$k)-0.05, ifelse(max(loo_dat$k) > 1, max(loo_dat$k), 1)), breaks = c(0, 0.5, 0.7)) +
    scale_x_continuous(breaks = seq(2008, 2018, 4)) +
    labs(y = "Pareto shape k", x = "Date") +
    theme_nature() 
  
  save_plot(loo_pl, file, w = w, h = h)
  return(loo_pl)
}

# Posterior summary
post_summary <- function(fit, pars) {
  fit$post_warmup_draws[,,pars] %>% 
    summarize_draws(mean, median, function(x) quantile2(x, probs = c(0.025, 0.975)), ess_bulk, Rhat) %>%
    mutate_if(is.numeric, round, 3)
}

# Sampler diagnostics
diagnostics <- function(fit) {
  np <- fit$post_warmup_sampler_diagnostics %>% 
    reshape2::melt() %>% 
    set_names(c("Iteration", "Chain", "Parameter", "Value")) %>%
    dplyr::select(Iteration, Parameter, Value, Chain)
  lp <- fit$post_warmup_draws[,,"lp__"] %>% reshape2::melt() %>% 
    dplyr::select(iteration, value, chain) %>% 
    set_names(c("Iteration", "Value", "Chain"))
  return(list(np = np, lp = lp))
}
