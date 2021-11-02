create_stan_data <- function(
  dat, # data
  target = "hiv_new_adult", # target age group for new HIV infections
  infectious = "hiv_living_adult", # target age group for people living with HIV
  incidence = "hiv_incidence_adult_mean", # target age group for hiv incidence
  ms = 6, # model start (including year=0)
  pa_method = "norm", # method to use for creating prior activity distribution (lnorm or dirichlet)
  pa_hp = NULL, # hyperparameters for activity distribution p_a
  pa_hp_diff = NULL # hyperparameter for difference between activity distributions p_a_diff
  ) 
{
  
  # set default pa_hp
  if (is.null(pa_hp)) {
    if (pa_method == "norm") {
      pa_hp <- c(2.0, 0.25, 0.5, 0.1)
      pa_hp_diff <- c(0.5, 0.25)
    }
    else if (pa_method == "dirichlet") {
      pa_hp <- rbind(c(rev(c(0.05, 0.15, 0.3, 0.3, 0.15, 0.05))), 
                     c(rev(c(0.10, 0.25, 0.3, 0.2, 0.1, 0.05))))
    }
  } else {
    if (pa_method == "norm") {
      if (length(pa_hp) != 4) stop("PA_HP has wrong dimension.")
    }
    else if (pa_method == "dirichlet") {
      if (ncol(pa_hp) != ms) stop("PA_HP has wrong dimension.")
    }
    
  }
  
  # all recipients
  all_recipients <- unique(dat$rcpt)
  
  # mean, upper and lower bound for HIV
  target_mean <- paste0(target, "_mean")
  target_lower <- paste0(target, "_lower")
  target_upper <- paste0(target, "_upper")
  infectious_mean <- paste0(infectious, "_mean")
  infectious_lower <- paste0(infectious, "_lower")
  infectious_upper <- paste0(infectious, "_upper")
  hiv_new <- dat[[target_mean]]
  hiv_new_lower <- dat[[target_lower]]
  hiv_new_upper <- dat[[target_upper]]
  hiv_living <- dat[[infectious_mean]]
  hiv_living_lower <- dat[[infectious_lower]]
  hiv_living_upper <- dat[[infectious_upper]]
  hiv_inci <- dat[[incidence]]
  dat <- dat %>%
    dplyr::select(-matches("hiv_")) %>%
    mutate(hiv_new = hiv_new,
           hiv_new_lower = hiv_new_lower,
           hiv_new_upper = hiv_new_upper,
           hiv_living = hiv_living,
           hiv_living_lower = hiv_living_lower,
           hiv_living_upper = hiv_living_upper,
           hiv_inci = hiv_inci) %>%
    na.omit() 
  
  # filter countries with not enough observations for model start
  dat <- dat %>%
    group_by(rcpt) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    dplyr::filter(n >= ms-1) %>%
    dplyr::select(-n)
  dat_mod <- dat %>%
    group_by(rcpt) %>%
    arrange(year) %>%
    slice(ms:n()) %>%
    ungroup()
  recipients_mod <- unique(dat$rcpt)
  if (length(recipients_mod) < length(all_recipients)) {
    removed_recipients <- all_recipients[!(all_recipients %in% recipients_mod)]
    removed_recipients <- paste(removed_recipients, collapse = ",")
    warning(sprintf("Recipients removed because #years < ms: %s", removed_recipients))
  }
  
  # filter observations with low incidence
  dat_mod <- dat_mod %>%
    group_by(rcpt) %>%
    mutate(suff_inci = ifelse(any(hiv_inci > 1), T, F)) %>%
    ungroup() %>%
    dplyr::filter(suff_inci) %>%
    dplyr::select(-suff_inci)
  recipients <- unique(dat_mod$rcpt)
  if (length(recipients) < length(recipients_mod)) {
    removed_recipients <- recipients_mod[!(recipients_mod %in% recipients)]
    removed_recipients <- paste(removed_recipients, collapse = ",")
    warning(sprintf("Recipients removed because of incidence: %s", removed_recipients))
    dat <- dat %>%
      dplyr::filter(rcpt %in% recipients)
  }
  
  # number of recipients
  R <- length(recipients)
  
  # year at which to start model = ms-1
  ms <- ms - 1
  y0 <- min(dat$year)
  yT <- max(dat$year)
  Y <- yT - y0 - ms + 1
  
  # initialize vectors for STAN
  N <- numeric()
  cum_N <- numeric()
  cum_N[1] <- 0
  hiv_new <- list()
  hiv_new_lower <- list()
  hiv_new_upper <- list()
  hiv_living <- list()
  hiv_living_lower <- list()
  hiv_living_upper <- list()
  pop <- list()
  std <- list()
  dom <- list()
  year <- list()
  mat <- list()
  edu <- list()
  gdp <- list()
  dens <- list()
  
  # create long vectors for STAN
  for (r in 1:R) {
    dat_r <- dplyr::filter(dat, rcpt == recipients[r])
    N[r] <- nrow(dat_r)
    cum_N[r+1] <- cum_N[r] + N[r]
    hiv_new[[cum_N[r]+1]] <- dat_r$hiv_new
    hiv_new_lower[[cum_N[r]+1]] <- dat_r$hiv_new_lower
    hiv_new_upper[[cum_N[r]+1]] <- dat_r$hiv_new_upper
    hiv_living[[cum_N[r]+1]] <- dat_r$hiv_living
    hiv_living_lower[[cum_N[r]+1]] <- dat_r$hiv_living_lower
    hiv_living_upper[[cum_N[r]+1]] <- dat_r$hiv_living_upper
    pop[[cum_N[r]+1]] <- dat_r$pop
    std[[cum_N[r]+1]] <- dat_r$std
    dom[[cum_N[r]+1]] <- dat_r$dom
    pop[[cum_N[r]+1]] <- dat_r$pop
    year[[cum_N[r]+1]] <- dat_r$year - y0 - ms + 1
    gdp[[cum_N[r]+1]] <- dat_r$fd_gdp
    edu[[cum_N[r]+1]] <- dat_r$fd_edu
    mat[[cum_N[r]+1]] <- dat_r$fd_mat
    dens[[cum_N[r]+1]] <- dat_r$fd_dens
  }
  
  # scale control variables
  controls <- matrix(cbind(unlist(gdp), unlist(edu), unlist(mat), unlist(dens)), ncol = 4)
  controls <- apply(controls, 2, scale)
  
  # combine funding
  funding <- matrix(cbind(unlist(std), unlist(dom)), ncol = 2)
  D <- ncol(funding)
  
  # rcpt identifier
  rcpt_id <- data.frame(
    cid = 1:length(recipients),
    group = recipients
  )
  
  # year identifier
  year_id <- data.frame(
    rid = unique(unlist(year)),
    x = unique(unlist(year)) + y0 + ms - 1
  )
  
  # create data list for STAN
  data_list <- list(
    R = R, Y = Y, N_tot = cum_N[R+1], max_N = max(N), N = N, cum_N = cum_N, model_start = ms, L = ms + 1,
    hiv_new = as.integer(unlist(hiv_new)), hiv_new_lower = unlist(hiv_new_lower), hiv_new_upper = unlist(hiv_new_upper), 
    hiv_living = unlist(hiv_living), hiv_living_lower = unlist(hiv_living_lower), hiv_living_upper = unlist(hiv_living_upper), 
    pop = unlist(pop),
    funding = funding , year = unlist(year), pop = unlist(pop), Z = controls,
    D = D, K = ncol(controls),
    p_a_hyper = pa_hp, p_a_hyper_diff = pa_hp_diff,
    dat_filt = dat, dat_model = dat_mod, rcpt_id = rcpt_id, year_id = year_id
  ) 
}
