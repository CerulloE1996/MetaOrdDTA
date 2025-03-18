**MetaOrdDTA** is an R package for the meta-analysis (MA) and network-meta-analysis (NMA) of medical (e.g., for diagnostic and screening) tests at all thresholds using ordinal regression-based models. 

Note that, eventhough **MetaOrdDTA** was originally made to analyse ordianl tests (such as queestionnaires), it can still be used for continuous tests (such as biomarkers) - using either semi-parameteric meta-analysis model proposed by Jones et al (Jones et al, 2019) - which was designed for continuous tests - or the ordinal regression based models (based on: Cerullo et al, 2022;Xu et al, 2013)

**Some facts about MetaOrdDTA:**
- All models are coded in **Stan**, using the cmdstanr R package.
- MCMC summary estimates of model parameters (including "raw" parameters - and also generated qantities such as sensitivity and specificy at each test threshold) can either be estimated using the cmdstanr ```summary```

**Coming soon:**
- **Covariates**: Inclusion of one (or more) covariates to conduct (network or single) meta-regression of test accuracy.
- 




An R package for the meta-analysis and network-meta-analysis of diagnostic and screening tests, optimised for ordinal tests such as questionnaires (with possibly missing threshold data is some or all studies). However, it can also be used for binary and continuous tests (e.g. using the Jones et al model). 

**To install MetaOrdDTA:**

First, install the cmdstanr R package (if not already installed) by running the following:
```
##
## ----------------- install cmdstanr first:
## Install the cmdstanr "outer" R package:
##
remotes::install_github("stan-dev/cmdstanr", force = TRUE)
#### Load cmdstanr R outer package:
require(cmdstanr)
#### Then, install the latest version of CmdStan:
install_cmdstan(cores = parallel::detectCores(),
                overwrite = TRUE,
                cpp_options = list("STAN_MODEL_LDFLAGS" = "-shared",   "CXXFLAGS" = "-fPIC"))
```

Then, you can install MetaOrdDTA by running the following R code:
```
remotes::install_github("https://github.com/CerulloE1996/MetaOrdDTA", force = TRUE)
```



**References:**

**Stan:**
Stan Development Team, 2024, Stan Reference Manual, version 2.36.0. https://mc-stan.org 

**cmdstanr** R package (to run Stan / NUTS-HMC):
Gabry J, Češnovar R, Johnson A, Bronder S (2024). cmdstanr: R Interface to 'CmdStan'. R package version 0.8.1. https://github.com/stan-dev/cmdstanr.

**BayesMVP R package:**
https://github.com/CerulloE1996/BayesMVP

**Cerullo et al, 2024:**
https://onlinelibrary.wiley.com/doi/10.1002/jrsm.1567
https://arxiv.org/abs/2103.06858 [may be more up to date in the future]

**Nyaga et al, 2016:**
https://journals.sagepub.com/doi/10.1177/0962280216669182

**Jones et al, 2019:**
https://pubmed.ncbi.nlm.nih.gov/31571244/

**Xu et al, 2013:**
https://pubmed.ncbi.nlm.nih.gov/23212851/

**Rutter & Gatsonis, 2001:**
https://onlinelibrary.wiley.com/doi/10.1002/sim.942

**SNAPER-HMC (Sountsov & Hoffman, 2022) :**
https://arxiv.org/abs/2110.11576v1

