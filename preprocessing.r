library(tidyverse)
library(reshape2)
library(stringi)
library(zoo)

start_year = 2003
end_year = 2017

#### Control variables ####

# variables
ctrl_vars <- c("gdp", "mat", "dens", "pop", "edu")
names(ctrl_vars) <- c("GDP per capita (current US$)",
                      "Maternal mortality ratio (modeled estimate, per 100,000 live births)",
                      "Population density (people per sq. km of land area)",
                      "Population, total",
                      "School enrollment, primary (% gross)")

# load data
controls <- read.csv("data/raw/socioeconomic_health_indicators_worldbank.csv") %>%
  rename(country = `Country.Name`,
         iso3 = `Country.Code`,
         series = `Series.Name`) %>%
  dplyr::select(-`Series.Code`) %>%
  dplyr::filter(series %in% names(ctrl_vars)) %>%
  mutate(series = recode(series, !!! ctrl_vars)) %>%
  melt(id.vars = c("country", "iso3", "series")) %>%
  mutate(year = as.numeric(stri_extract(variable, regex = "\\d{4}"))) %>%
  dplyr::filter(year >= 2002 & year <= end_year) %>%
  mutate(value = as.numeric(value)) %>%
  dplyr::select(-variable)

# imputation
impute <- function(x, global_mean_x, min_obs = 3) {
  # impute by mean if more than min_obs
  if (sum(!is.na(x)) < min_obs) {
    y <- global_mean_x
  }
  # impute linearly
  y <- na.approx(x, rule = 2, na.rm = T)
  # impute by previous value
  y <- na.locf(y, fromLast = T, na.rm = T)
  # impute by lagged value
  y <- na.locf(y, fromLast = F, na.rm = T)
  return(y)
}

controls <- controls %>%
  group_by(series, year) %>%
  mutate(mean_value = mean(value, na.rm = T)) %>%
  ungroup() %>%
  group_by(country, iso3, series) %>%
  arrange(year) %>%
  mutate(value = impute(value, mean_value)) %>%
  ungroup() %>%
  dplyr::select(-mean_value)
  
# first differences
controls <- controls %>%
  dcast(country + iso3 + year ~ series) %>%
  group_by(country, iso3) %>%
  arrange(year) %>%
  mutate_at(vars(mat, gdp, dens), list(fd = ~log(. / dplyr::lag(.)))) %>%
  mutate(fd_edu = edu / 100 - dplyr::lag(edu / 100)) %>%
  rename_at(vars( contains( "_fd") ), list( ~paste("fd", gsub("_fd", "", .), sep = "_") )) %>%
  ungroup() %>%
  dplyr::filter(year >= start_year)
  
# save data
write.csv(controls, "data/preprocessed/control_variables.csv", row.names = F)

  
#### HIV ####

# convert censored data to numeric
hiv_to_numeric <- function(x) {
  x <- gsub("<", "", x) # caveat: consider it to be the maximum value
  x <- gsub(" ", "", x)
  x <- as.numeric(x)
  return(x)
}

# wrapper function to read HIV data
read.hiv <- function(file, col_name) {
  read.csv(file, stringsAsFactors = F) %>%
    dplyr::select(!matches("Footnote")) %>%
    reshape2::melt("Country") %>%
    mutate(year = as.integer(stringi::stri_extract(variable, regex = "\\d{4}"))) %>%
    mutate(variable = ifelse(grepl("lower", variable), "lower", ifelse(grepl("upper", variable), "upper", "mean"))) %>%
    dplyr::filter(year > 2002 & year <= end_year) %>%
    reshape2::dcast(Country + year ~ variable) %>%
    set_names(c("country", "year", paste0(col_name, "_", tail(colnames(.), 3)))) %>%
    mutate_at(vars(matches("hiv_")), .funs = list(cens = ~ifelse(grepl("<", .), hiv_to_numeric(.), NA))) %>%
    mutate_at(vars(matches("hiv_")), hiv_to_numeric)
}

# read and join data
hiv <- read.hiv("data/raw/hiv_unaids/hiv_new_infections_all.csv", "hiv_new_all") %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_new_infections_0-14.csv", "hiv_new_children")) %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_new_infections_10-19.csv", "hiv_new_adolescent")) %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_new_infections_15-24.csv", "hiv_new_young")) %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_new_infections_15-49.csv", "hiv_new_mid")) %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_new_infections_15+.csv", "hiv_new_adult")) %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_new_infections_50+.csv", "hiv_new_old")) %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_living_all.csv", "hiv_living_all")) %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_living_15-24.csv", "hiv_living_young")) %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_living_15-49.csv", "hiv_living_mid")) %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_living_15+.csv", "hiv_living_adult")) %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_living_50+.csv", "hiv_living_old")) %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_incidence_all.csv", "hiv_incidence_all")) %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_incidence_15-24.csv", "hiv_incidence_young")) %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_incidence_15-49.csv", "hiv_incidence_mid")) %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_incidence_15+.csv", "hiv_incidence_adult")) %>%
  left_join(read.hiv("data/raw/hiv_unaids/hiv_incidence_50+.csv", "hiv_incidence_old")) %>%
  mutate(country = ifelse(country == "Côte d'Ivoire", "Cote d'Ivoire", country)) %>%
  mutate(country = ifelse(country == "Democratic Republic of the Congo", "Congo, Dem. Rep.", country)) %>%
  mutate(country = ifelse(country == "Congo", "Congo, Rep.", country)) %>%
  mutate(country = ifelse(country == "Gambia", "Gambia, The", country)) %>%
  mutate(country = ifelse(country == "United Republic of Tanzania", "Tanzania", country))

# consider censored data as if value = censored value
hiv <- hiv %>%
  dplyr::select(-contains("_cens"))

# filter Sub-Saharan Africa (only Seychelles will be missing)
hiv <- hiv %>%
  dplyr::filter(country %in% unique(controls$country))

# save data
write.csv(hiv, "data/preprocessed/hiv.csv", row.names = F)


#### ODA ####

# load data
oda <- read.csv("data/raw/oda_ihme.csv")  %>%
  dplyr::select(recipient_country, recipient_isocode, year, wb_regioncode, matches("hiv_")) %>%
  rename(country = recipient_country, iso3 = recipient_isocode) %>%
  dplyr::filter(wb_regioncode == "SSA") %>%
  dplyr::select(-wb_regioncode) %>%
  dplyr::filter(year >= start_year & year <= end_year)

# summarize data
oda <- oda %>%
  melt(id.vars = c("country", "iso3", "year")) %>%
  mutate(value = as.numeric(value)) %>%
  group_by(country, iso3, year, variable) %>%
  summarize(value = 1e3 * sum(value, na.rm = T)) %>%
  ungroup() %>%
  mutate(variable = as.character(variable)) %>%
  mutate(variable = ifelse(variable == "hiv_dah_20", "total", variable)) %>%
  mutate(variable = gsub("hiv_", "", variable)) %>%
  mutate(variable = gsub("_dah_20", "", variable)) %>%
  dcast(country + iso3 + year ~ variable) %>%
  dplyr::select(country, iso3, year, total, everything()) %>%
  dplyr::filter(year >= start_year & year <= end_year)

# save data
write.csv(oda, "data/preprocessed/oda_funding.csv", row.names = F)


#### Domestic funding ####

# load data
dom <- read.csv("data/raw/domestic_funding_ihme.csv") %>%
  dplyr::select(location_name, iso3, year, ghes_total_mean, ghes_total_lower, ghes_total_upper) %>%
  dplyr::filter(iso3 %in% unique(controls$iso3)) %>%
  rename(country = location_name) %>%
  dplyr::filter(year >= start_year & year <= end_year) %>%
  mutate_at(vars(matches("ghes")), ~ . * 1e3)

# translation of value (dom in 2019 USD and oda is in 2020 USD)
# based on https://www.in2013dollars.com/us/inflation/2019?endYear=2020&amount=1
dom <- dom %>%
  mutate_at(vars(matches("ghes")), ~ . * 1.01) 

# save data
write.csv(dom, "data/preprocessed/domestic_funding.csv")



#### Population projections ####

# load data
ppop30 <- read_csv("data/raw/population_projections_un.csv") %>%
  mutate(country = ifelse(country == "Côte d'Ivoire", "Cote d'Ivoire",
                          ifelse(country == "Congo", "Congo, Rep.",
                                 ifelse(country == "United Republic of Tanzania", "Tanzania", 
                                        ifelse(country == "Democratic Republic of the Congo", "Congo, Dem. Rep.", 
                                               ifelse(country == "Gambia", "Gambia, The", country)))))) %>%
  mutate(medium_estimate_2030 = as.numeric(medium_estimate_2030 * 1e3)) %>%
  dplyr::filter(country %in% unique(controls$country))

# save data
write.csv(ppop30, "data/preprocessed/population_projections_2030.csv", row.names = F)



#### Join data ####

# join data
df <- hiv %>%
  dplyr::select(country, year, matches("adult"), matches("all")) %>%
  dplyr::filter(!is.na(hiv_new_adult_mean)) %>%
  left_join(controls, 
            by = c("country", "year")) %>%
  left_join(ppop30 %>% rename(pop30_est = medium_estimate_2030), 
            by = c("country")) %>%
  left_join(oda %>% dplyr::select(iso3, year, total) %>% rename(std = total), 
            by = c("iso3", "year")) %>%
  left_join(dom %>% dplyr::select(iso3, year, ghes_total_mean) %>% rename(dom = ghes_total_mean), 
            by = c("iso3", "year")) %>%
  rename(rcpt = iso3) %>%
  dplyr::select(country, rcpt, year, matches("hiv"), std, dom, everything()) %>%
  mutate(country = ifelse(rcpt == "COG", "Republic of the Congo", 
                          ifelse(rcpt == "COD", "Democratic Republic of the Congo",
                                 ifelse(rcpt == "GMB", "Gambia", country))))

# save data
write.csv(df, "data/preprocessed/data.csv", row.names = F)
