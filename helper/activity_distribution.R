lnorm.mu <- function(m, s) {
  log( m ^ 2 / sqrt( s ^ 2 + m ^ 2 ) )
}

lnorm.sigma <- function(m, s) {
  sqrt( log( s ^ 2 / m ^ 2 + 1 ) )
}

p_a_lnorm <- function(
  x, # lag
  max_x, # maximum lag
  m, # mu_p_a
  s # sigma_p_a
  ) {
  log_mu = lnorm.mu(m, s)
  log_sigma = lnorm.sigma(m, s)
  p_tot <- plnorm( max_x + 1, log_mu, log_sigma)
  if (x == 0) {
    return( plnorm( 1, log_mu, log_sigma) / p_tot )
  } else {
    return( (plnorm( x + 1, log_mu, log_sigma) - plnorm( x - 1, log_mu, log_sigma)) / p_tot )
  }
}

p_a_norm <- function(
  x, # lag
  max_x, # maximum lag
  m, # mu_p_a
  s # sigma_p_a
) {
  p_tot <- pnorm( max_x + 1, m, s)
  if (x == 0) {
    return( pnorm( 1, m, s) / p_tot )
  } else {
    return( (pnorm( x + 1, m, s) - pnorm( x - 1, m, s)) / p_tot )
  }
}

p_a.prior <- function(
  method = "dirichlet", # type of prior
  pa = list(lambda = 1e2, w = c(0.05, 0.15, 0.3, 0.3, 0.15, 0.05)), # p_a hyperparameters 
  pa_diff = NULL, # difference to original distribution
  D = 4000, # number of draws
  X = 5 # number of lags including zero (0:X)
  ) {
  
  if (grepl("norm", method)) {
    mu <- rnorm(D, pa[1], pa[2])
    sigma <- rnorm(D, pa[3], pa[4])
    if (!is.null(pa_diff)) {
      mu_diff <- rnorm(D, pa_diff[1], pa_diff[2])
    }
  } else if (method == "dirichlet") {
    #log_mu_k <- lnorm.mu(pa$mu, pa$sigma)
    #log_sigma_k <- lnorm.sigma(pa$mu, pa$sigma)
    w <- pa$w
    k <- rep(pa$lambda, D) #rexp(D, pa$lambda)
  }
  
  
  PA <- matrix(NA, nrow = D, ncol = X+1)
  for (d in 1:D) {
    for (x in 0:X) {
      if (method == "lnorm") {
        PA[d,x+1] <- p_a_lnorm(x, max_x = X, m = mu[d], s = sigma[d])
      } else if (method == "norm") {
        if (!is.null(pa_diff)) {
          PA[d,x+1] <- p_a_norm(x, max_x = X, m = mu[d] + mu_diff[d], s = sigma[d])
        } else {
          PA[d,x+1] <- p_a_norm(x, max_x = X, m = mu[d], s = sigma[d])
        }
      }
    }
    if (method == "dirichlet") {
      PA[d,] <- c(DirichletReg::rdirichlet(1, w * k[d]))
    }
  }
  
  PA <- data.frame(PA) %>%
    mutate(.draw = 1:nrow(.)) %>%
    reshape2::melt(id.vars = ".draw") %>%
    rename(year = variable,
           p_a = value) %>%
    mutate(year = rep(0:X, each = D)) %>%
    dplyr::select(-.draw) %>%
    group_by(year) %>%
    summarize(
      .mean = mean(p_a),
      .lower = quantile(p_a, 0.025),
      .upper = quantile(p_a, 0.975)
    ) %>%
    ungroup() 
  
  
  
  return(PA)
}


p_a.posterior <- function(model_fit, L, d, D) {
  m <- numeric()
  l <- numeric()
  u <- numeric()
  for (i in 1:L) {
    if (D > 1) {
      ps <- post_summary(model_fit, paste0("p_a[",d,",",i,"]"))
    } else {
      ps <- post_summary(model_fit, paste0("p_a[",i,"]"))
    }
    m[i] <- ps$mean
    l[i] <- ps$`q2.5`
    u[i] <- ps$`q97.5`
  }
  p_a_post <- data.frame(
    year = (L-1):0,
    .mean = m,
    .lower = l,
    .upper = u
  )
  return(p_a_post)
}

p_a.posterior_plot <- function(model_fit, L, D, D_names = c("Development assistance", "Domestic funding")) {
  p_a_post_D <- do.call(rbind, lapply(1:D, function(x) p_a.posterior(model_fit, L, d = x, D = D)))
  p_a_post_D$funding_type <- factor(rep(D_names, each = nrow(p_a_post_D) / 2), levels = D_names)
  ggplot(p_a_post_D, aes(x = year)) +
    geom_line(aes(y = .mean, color = funding_type)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = funding_type), alpha = .2) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,NA)) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    labs(x = "t (years)", y = expression(p[RA]*"(t)")) +
    theme_bw2() +
    theme(legend.position = "top", legend.title = element_blank())
}

p_a.prior_vs_posterior <- function(model_fit, L, d = 1, D = 2, method = "norm", pa = c(2.0, 0.25, 0.5, 0.1), pa_diff = c(0.5, 0.25), draws = 4000) {
  p_a_post <- p_a.posterior(model_fit, L, d = d, D = D)
  p_a_prior <- p_a.prior(method = method, pa = pa, pa_diff = pa_diff, X = L-1, D = draws)
  p_a <- rbind(p_a_post,p_a_prior) %>%
    mutate(type = factor(rep(c("Posterior", "Prior"), each = nrow(p_a_post)), levels = c("Posterior", "Prior"))) 
  ggplot(p_a, aes(x = year)) +
    geom_line(aes(y = .mean, color = type)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = type), alpha = .2) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,NA)) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "s (years)", y = "p_a(s)") +
    theme_bw2() +
    theme(legend.title = element_blank(), legend.position = "top")
}
