**MetaOrdDTA** is an R package for the meta-analysis (MA) and network-meta-analysis (NMA) of medical (e.g., for diagnostic and screening) tests at all thresholds (including unreported or "missing" thresholds) using ordinal regression-based models. 


Unlike most methods, it allows users to input ALL data in the analysis and still produce summary measures of sensitivty (Se) and specificity (Sp) at every test threshold - rather than "partitioning" the data at each threshold (which does not make full use of the data - unless every study reports data at every single test threshold - which is unfortunately still very rare).




# 1. Models:
MetaOrdDTA can fit the following models:

- Ordinal regression-based models (based on our previous work [Cerullo et al, 2022]). At present, users can fit 2 broad types of ordinal-regression -based MA/NMA models using MetaOrdDTA:
  - THe first (and currently the default - but this may change soon depending on our findings) is the **"Xu"** parameterisation (which you can select by specifying: ```model_parameterisation = "Xu"```). This is based on the work of Xu et al, 2013 (but extended to any number of thresholds and extended to meta-analysis, as well as other differences - which we will detail in the upcoming paper). This model does not model a scale parameter at all (assumes scale = 1 in both diseased and non-diseased groups) and instead allows both location parameters to be completely "free" (i.e., freely estimated) in both groups, and models the between-study correlation **explicitly** using a bivariate model - much like the commonly-used "bivariate model" proposed by Reitsma et al (Reitsma et al, 2005).
  - The second option is the **"R&G" parameterisation** (which you can select by specifying: ```model_parameterisation = R&G```). This model is based on another widely used model for the meta-analysis of test accuracy - typically referred to as the "HSROC" model (Rutter & Gatsonis, 2001). Our model implemented in MetaOrdDTA extends this model to the "multiple threshold" case, and so is a "true" ordinal regression model - unlike the original HSROC model (Rutter & Gatsonis, 2001). Furthermore, just like the standard HSROC model, this model does not model the between-study correlation explicitly - instead it models this correlation **implicitly** by using **shared** location and scale parameters in both groups (but not identical in both groups - namely we assume that the location in the diseased group is $+\beta$ and $-\beta$ in the diseased and non-diseased groups, respectively.
  - Both the "Xu" and "R&G" -based models can be fit using either fixed-effects or random-effects thresholds/cutpoints (please see **section 3** below for more details on this).
- MetaOrdDTA can also fit the "multiple thresholds" model proposed by **Jones et al (Jones et al, 2019)**. This can be done by specifying ```cts = TRUE``` and/or ``` model_parameterisation = "Jones"```. This model is more appropriate for continuous data as it makes more restrictive assumptions about the data compared to the models proposed/mentioned above (see below/upcoming paper for more details).
- All of the aformentioned models can be run either for meta-analysis (MA; i.e. the analysis of a **single index test**) or for **network-meta-analysis (NMA)**. The NMA versions of these models are based on the NMA model proposed by Nyaga et al (Nyaga et al, 2016).




# 2. More about MetaOrdDTA:
Unlike the more "standard", commonly-used methods of meta-analysis of test accuracy - such as the "bivariate" model (Reitsma et al, 2005) or the "HSROC" model (Rutter & Gatsonis, 2001) - MetaOrdDTA allows users to include * all * studies in the analysis, regardless of the test threshold reported ** AND ** **simultanously** estimate test accuracy (sensitivity [Se] and specificity [Sp]) at **every test threshold** - even if a given study does not report accuaracy at any given test threshold (i.e. at "missing" thresholds). 

Furthermore, it uses a true ordinal regression-based model - based on a previous model we proposed (Cerullo et al, 2022). Unlike other proposed "multiple threshold" (or "missing threshold") methods - this method does not make strong assumptions about the data, such as assuming it is continuous (parameteric or "semi-parameteric"). This means that MetaOrdDTA is likely to produce better summary estimates compared to other propsoed "multiple threshold" DTA methods (and our preliminary simulations based on real-world data support this - see below for more details). 

# Why should you use MetaOrdDTA for analysing test accuracy data at multiple thresholds? (based on preliminary findings):
Note that we are currently writing a paper which contains a real-life-based simulation study which suggests that for ordinal tests (even ones which have 25+ ordinal categories) MetaOrdDTA may produce notably better summary sensitivty and specificity estimates compared to other recently proposed methods (see **table 1** below for more details), due to making much less restrictive (and more realistic) assumptions about the data (i.e., modelling it as it really is - ordinal).

The preliminary findings suggest that these differences in summary Se and Sp may sometimes be quite substantial. To be more specific, thus far our (preliminary) findings suggest that our proposed model (based on Cerullo et al, 2022; Jones et al, 2019; Xu et al, 2013, and Rutter & Gatsonis, 2001) may obtain mean absolute differences in summary Se of **over 5%** and a maximum absolute difference of **around 8% or more** on our "real-life" simulated data - which is based on simulated (but realistic) "real-world-based" data we obtained on a PHQ-9 screening test for depression (which is scored from 0-27 and hence has 28 ordered categories and 27 thresholds).  

The table below shows some **priliminary** results (mentioned above) obtained from our simulation study:

# (PRELIMINARY) Model Performance on simulated data (based on real-life PHQ-9 data w/ 28 ordinal categories)

| Model | Mean Absolute Difference |  | Max. Absolute Difference |  |
| --- | :---: | :---: | :---: | :---: |
|  | **Se (%)** | **Sp (%)** | **Se (%)** | **Sp (%)** |
| Jones (log-probit)          | 4.3 | 3.0 | 11.4 | 8.8 |
| Jones (Box-Cox)             | 5.6 | 1.5 | 11.7 | 4.0 |
| Cerullo (fixed cutpoints)   | 1.2 | 0.8 | 3.9  | 1.4 |
| Cerullo (random cutpoints)  | 1.7 | 0.6 | 3.8  | 1.4 |




# 3. Some other facts about MetaOrdDTA:

- All of the ordinal models (both Xu-based [Xu et al, 2013] and R&G-based [Rutter & Gatsonis, 2001] models) make use of one of the following "induced Dirichlet" cutpoint models, depending on whether you use a fixed-effect (aka "complete-pooling") or a random-effect (aka "partial-pooling") between study model on the cutpoint parameters (for more information on the induced-Dirichlet model please refer to Michael Betancourt's "ordinal regression" case study (Betancourt et al, 2019). More specifically:
  - If you decide to use a **fixed-effect cutpoint model** (i.e. you specify ```random_thresholds = FALSE``` - which is the default as it currently standads), then MetaOrdDTA will use induced Dircihlet model as a ** prior ** model. This prior allows users to put priors directly on the "induced" ordinal probabilities (which are transformed parameters - not the "raw" cutpoint parameters). This makes prior modelling much more intuitive, and also allows one to set a "flat" prior on the cutpoints (which is otherwise virtually impossible to do). Furthermore, it also allows users to much more easily input domain-specific prior knowledge into the model.
  - On the other hand, if you decide to use a **random-effect cutpoint model** (i.e. you directly model the between-study heterogenity in the cutpoints - by specifying ```random_thresholds = TRUE```), then MetaOrdDTA will use the induced-Dirichlet model as part of the between-study model itself (i.e. it will form part of the likelihood itself - not "only" the prior). 
  - Note that the default prior for the ordinal models in MetaBayesDTA is a "flat" dirichlet (i.e., a prior of ${1,...,1}$ for the "alpha" concentration parameters of the Dirichlet distribution). 
- All models are coded in **Stan**, using the cmdstanr R package.
- MCMC summary estimates of all model parameters (including "raw" parameters - and also generated qantities such as sensitivity and specificy at each test threshold) can either be estimated using the cmdstanr ```summary$()``` method, or using the more efficient (sometimes by over a factor of 10) ```BayesMVP::generate_summary_tibble()``` function.




# 4. Coming soon:
- **Covariates**: Inclusion of one (or more) covariates to conduct (network or single) meta-regression of test accuracy.
- **Faster estimation using BayesMVP:** I am currently working on implementing a second MCMC sampling algorithm (based on SPANER-HMC; see Sountsov & Hoffman, 2022) into MetaOrdDTA - which is the same algorithm my other R package (BayesMVP; see: https://github.com/CerulloE1996/BayesMVP) uses. 

Note that, eventhough **MetaOrdDTA** was originally made to analyse ordianl tests (such as queestionnaires), it can still be used for continuous tests (such as biomarkers) - using either semi-parameteric meta-analysis model proposed by Jones et al (Jones et al, 2019) - which was designed for continuous tests - or the ordinal regression based models (based on: Cerullo et al, 2022; Xu et al, 2013).




# 5. How to install the MetaOrdDTA R package:

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






# 6. References:

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

**Reitsma et al, 2005:**
https://www.jclinepi.com/article/S0895-4356(05)00162-9/abstract

**SNAPER-HMC (Sountsov & Hoffman, 2022) :**
https://arxiv.org/abs/2110.11576v1

