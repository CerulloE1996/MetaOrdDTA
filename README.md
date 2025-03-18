**MetaOrdDTA** is an R package for the meta-analysis (MA) and network-meta-analysis (NMA) of medical (e.g., for diagnostic and screening) tests at all thresholds using ordinal regression-based models. 

Note that, eventhough **MetaOrdDTA** was originally made to analyse ordianl tests (such as queestionnaires), it can still be used for continuous tests (such as biomarkers) - using either semi-parameteric meta-analysis model proposed by Jones et al (Jones et al, 2019) - which was designed for continuous tests - or the ordinal regression based models. The latter may give 

Some facts about MetaOrdDTA:
- All models are coded in **Stan**, using the cmdstanr R package.
- MCMC summary estimates of model parameters (including "raw" parameters - and also generated qantities such as sensitivity and specificy at each test threshold)

**Coming soon:**
- **Covariates**: Inclusion of one (or more) covariates to conduct (network or single) meta-regression of test accuracy.
- 




An R package for the meta-analysis and network-meta-analysis of diagnostic and screening tests, optimised for ordinal tests such as questionnaires (with possibly missing threshold data is some or all studies). However, it can also be used for binary and continuous tests (e.g. using the Jones et al model). 



      ## path to Stan model
      file <- file.path(pkg_dir, "inst/stan_models/LC_MVP_bin_w_mnl_cpp_grad_v1.stan") 
      ## path to the C++ .hpp header file
      path_to_cpp_user_header <- file.path(pkg_dir, "src/lp_grad_fn_for_Stan.hpp") 
      ## compile model together with the C++ functions
      mod <- cmdstan_model(file,  force_recompile = TRUE, user_header = path_to_cpp_user_header) 




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

