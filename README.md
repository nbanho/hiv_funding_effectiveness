# Estimating the effect of development assistance and domestic funding on new infections with HIV in Sub-Saharan Africa

This repository contains code and data to reproduce the results from the paper "Estimating the effect of development assistance and domestic funding on new infections with HIV in Sub-Saharan Africa". 

## Data 

A preprocessed data file is provided in `data/preprocessed/data.csv`. The file includes the following variables:
- rcpt = ISO-3 country code
- year = year
- std = official development assistance (ODA) funding for sexually transmitted diseases (USD)
- dom = domestic funding for HIV (USD)
- edu = education
- gdp = gross domestic product (USD per capita)
- pop = population
- dens = log population density
- mat = maternal mortality per capita
- fd_* = first difference of control variables
- hiv_new_* = new infections with HIV by age group
- hiv_living_* = people living with HIV by age group
- hiv_incidence_* = incidence with HIV by age group (new infections per 1,000 population)


## Model

Models were developed in the probabilistic programming language STAN and provided in the folder `models`. The model used for the main analysis is `active_funding_notr_me.stan`. 

- To run any model, you need to have `CmdStan` installed together with `CmdStanR` in order execute STAN models from R. 
- To run the main model from Terminal in Linux, direct into the main directory of the repository and run the following command: `Rscript main_model.R active_funding_notr_me hiv_new_adult 6`. 
- Analogously, to run the full sensitivity analysis, use the following command: `Rscript sensitivity_analysis.R all`

By default, results from the main analysis will be saved into `fitted-models/main` and results from the sensitivity analysis will be saved into `fitted-models/sens`.


## Analysis

Descriptive analysis, analysis of the model estimation results, and the sensitivity analysis are documented and conducted in the R Markdown file `analysis.Rmd`. 

- Descriptive results can be reproduced. 
- Results from the main analysis can also be reproduced with the uploaded model estimation results in `fitted-models/main`. 
- To reproduce the results from the sensitivity analysis, first obtain the model outputs from the full sensitivity analysis, which could not all be uploaded into the repository to conserve space.

By default, results will be saved into the folder `results`.

