# Spillover effects of reactive, focal malaria interventions

## Directory structure
**`0-config.R` :** configuration file that sets data directories, sources base functions, and loads required libraries, as well as the check_washout functions for matching

**`0-base-functions` :** folder containing R scripts with general functions used across the analysis
* `0-base-functions.R`: contains miscellaneous custom functions
* `1-helper-functions-covar`: contains helper functions generating individual level covariates for analaysis data
* `2-namibia-data-prep-combine-dist.R` : contain a helper function that combine the dist variables
* `3-prescreen-function.R` : containing pre-screening function using a likelihood ratio test.
* `4-hTMLE.R` : the main script containg the implementation of hierachical TMLE based on Laura Balzer's paper.
* `5-wrapper_hTMLE.R` : a wrapper function consist 1) pre-screening, 2) adaptively select Q model, and 3) do hTMLE.
* `6-incidence-wrapper-functions.R` : wrapper functions for incidence analyses 
* `7-glm-functions.R` : functions for unadjusted GLM models
* `8-adaptively-choose-Qdata.R` : function for adaptively selecting whether Q model is fit on individual or cohort level data
* `9-tmle-prevalence-function.R` : wrapper functions for prevalence analyses with TMLE
* `10-wrapper_hTMLE_direct.R` : wrapper functions for sensitivity direct effects analyses with TMLE


**`1-data-cleaning` :** folder containing data processing scripts  

* `1-namibia` : folder containing data processing scripts for Namibia  
  * `1-clean-namibia-baseline.R` : clean baseline data from geographical reconaissance survey and create a dataset with one row for each person in the study 
  * `2-clean-namibia-intervention.R` : clean dataset for intervention recipients 
  * `3-clean-namibia-index.R` : clean dataset for index cases
  * `4a-match-index-int.R` : match index cases that triggered interventions to individual imputed dataset 
  * `4b-match-int.R` : match intervention individuals to individual imputed dataset 
  * `4c-match-index-nonint.R` : match index cases that did not trigger interventions to individual imputed dataset 
  * `4d-cohort-subsets.R` : create matched cohort subsets
  * `4e-cohort-classification.R` : classify cohorts using different radii for spillover zones
  * `5-clean-namibia-xs.R` : clean cross-sectional survey dataset for prevalence analyses
  * `6a-data-prep.R` : generate the analysis data set with arm-specific observation period 
  * `6b-data-prep-short.R` : generate the analysis data set using short observation period
  * `6c-data-prep-long.R` : generate the analysis data set using long observation period
  * `7a-data-prep-short-sens-obs.R` : generate the analysis data set using short observation period (sensitivity analysis with different short observation period)
  * `7b-data-prep-long-sens-obs.R` : generate the analysis data set using long observation period (sensitivity analysis with different short observation period)
  * `7c-data-prep-short-sens-spzone-2km.R` : generate the analysis data set using short observation period (sensitivity analysis with 2km spillover zone radius)
  * `7d-data-prep-long-sens-spzone-2km.R` : generate the analysis data set using long observation period (sensitivity analysis with 2km spillover zone radius)
  * `7e-data-prep-short-sens-spzone-3km.R` : generate the analysis data set using short observation period (sensitivity analysis with 3km spillover zone radius)
  * `7f-data-prep-long-sens-spzone-3km.R` : generate the analysis data set using long observation period (sensitivity analysis with 3km spillover zone radius)
  * `8-data-prep-sens-no-overlap.R` : generate the analysis data set with no overlapping target areas / spillover zones
    
**`2-analysis` :** folder containing data analysis scripts 
* `0-namibia-incidence-learners.R` : define learners for SuperLearner analyses
* `1-preprocess-inc` : folder with scripts that pre-process data for incidence analysis  
* `2-incidence-tmle` : folder with scripts for primary analyses using hTMLE for incidence data
* `2b-incidence-tmle-adaptQ` : folder with scripts that select the Q level for hTMLE analyses
* `3-incidence-tmle-sens` : folder with scripts for sensitivity analyses using hTMLE for incidence data
* `4-incidence-tmle-EM` : folder with scripts for effect modification analyses using hTMLE for incidence data
* `5-prevalence` : folder with scripts for secondary analyses using prevalence data 
* `6-inc-glm-indiv-unadj.R` : unadjusted estimates of effects on incidence using GLM
* `7-namibia-btw-clus-spillover-prev.R` : assess between cluster spillover effects using prevalence data
* `7-namibia-btw-clus-spillover.R` : assess between cluster spillover effects using incidence data
* `8-dist-to-HF.R` : assess whether incidence and prevalence vary by distance to healthcare facility
* `9-cost-effectiveness.R` : cost-effectiveness analysis
* `10-percent-treated.R` : calculate percentage of individuals treated inside and outside target area


**`3-figure-scripts` :** folder containing figure generation scripts  

**`4-table-scripts` :** folder containing table generation scripts  

**`5-figures` :** folder containing figures

**`6-tables` :** folder containing tables

