# GrayWhaleCalfProduction
This repository contains R code to extract and analyze gray whale calf production data that are collected at Piedras Blancas, CA, annually. 
Necessary functions are stored in GrayWhaleCalfProduction_fcns.R. Some of the functions in the file may be obsolete and unused. The analytical
part is conducted via JAGS (Just Another Gibbs Sampler, which can be downloaded from https://sourceforge.net/projects/mcmc-jags/files/. 

It is assumed that this repository is downloaded to your local computer via RStudio and a project created. Please create three empty folders
within the project folder and name them "data," "figures," and "RData." 

In 2023, new data extraction code was developed due to some inconsistencies that were found in the previous version (Collating and Formatting Passdown.R).
The comparison of the two extraction methods is provided in a notebook file named "compare_data_extraction.Rmd." In the past report, I did not 
use the new extraction method. I appended the new estimate for the 2023 season to the 2022 report. For the 2023 report, I did use the new extraction but
the table was replaced with the old version. These reports can be found in Report_calf_production_2022.Rmd and Report_calf_production_2023.Rmd.

The order of the data processing and analysis should be:
1. Raw data should be saved in the "data/All data" directory. Folder names are C_C YYYY, where YYYY indicates a four-digit year. 
2. Extract Excel Data.Rmd
3. Run calf_production_jags.R
4. Edit Report_calf_production_2023.Rmd to create a report.

In Stewart and Weller (2021), the mean ($\lambda$) of the number of whales passing by the survey area is assumed to be the same within a week. 
There is no assumption about the mean except it is bounded between 0 and 40 (i.e., $ \lambda \sim UNIF(0,40)$). The weekly mean assumption is a bit arbitrary. 
In general, the number of mother-calf pairs increases at the beginning of each survey season and reaches a peak, then decreases thereafter. 
So, it may make sense to make an assumption about means that change daily. 

It was also assumed that the true number of whales in the survey area per 3-hr period was an independent Poisson random deviate (i.e., $n_{true} \sim POI(\lambda)$).  
I think that there should be some auto-correlations between two consecutive observation periods; they are arbitrarily separated by changes in observers. 
Each shift is conducted by two observers, which lasts up to 3 hrs (1.5 hr x 2 observers). It's likely that the number of whales passing by 
the survey area would be similar between two consecutive shifts than those that are further apart. 

The observed number of whales per shift is assumed to be a binomial random deviate (i.e., $n_{obs} \sim BIN(n_{true}, p_{obs})$), where $p_{obs}$ was 
estimated through a calibration study ($p_{obs} \sim N(0.889, 0.06375^2)$). Looking at the data, there are many zeros. So, it may be better to use 
some other distributions (e.g., zero-inflated Poisson, zero-inflated negative binomial, Tweedie) to model these data. It's not easy to get zero observations 
unless there are no whales to observe when the detection probability is close to 0.9. Furthermore, there probably are other factors affecting the observations. 
Perhaps, the sighting conditions may be used to model "actual" sighting probability (i.e., $ p_{obs}$ = f(sighting condition) and $p_{obs}$ reaches 0.889 as 
the sighting condition becomes ideal).
