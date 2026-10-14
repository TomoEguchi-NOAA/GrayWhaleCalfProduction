# GrayWhaleCalfProduction

This repository contains R code to extract and analyze gray whale calf production data that are collected at Piedras Blancas, CA, annually. 


Necessary functions are stored in GrayWhaleCalfProduction_fcns.R. Some of the functions in the file may be obsolete and unused. The analytical part was conducted using JAGS (Just Another Gibbs Sampler, which can be downloaded from https://sourceforge.net/projects/mcmc-jags/files/). The model was developed by Josh Stewart. In 2026, an updated analytical approach was introduced due to the previous model's failure when counts declined and were spread over several weeks. The new approach uses the Negative Binomial model instead of the Poisson-Binomial and is coded in Stan (https://mc-stan.org/). The updated model was developed by Tomo Eguchi.

It is assumed that this repository is downloaded to your local computer and a project is created in RStudio. Please create three empty folders within the project folder and name them "data," "figures," and "RData". They are case sensitive. 

In 2023, new data extraction code was developed due to inconsistencies found in the previous version (Collating and Formatting Passdown.R). The comparison of the two extraction methods is provided in an RStudio notebook file named "compare_data_extraction.Rmd." The new extraction method was first used for the 2023 season. For the 2023 report, even though the new extraction method was used, the table was replaced with the results from the old extraction version. These reports can be found in Report_calf_production_2022.Rmd and Report_calf_production_2023.Rmd. The new extraction method was reported in Lang et al. (2024). 

The order of the data processing and analysis should be: 

1. Raw data should be saved in the "data/All data" directory. Folder names are C_C YYYY, where YYYY indicates a four-digit year. 
2. Edit list.sheet.names.inshore, list.col.types.inshore, and list.col.names.inshore in GrayWhaleCalfProduction_fcns.R (if you pull the most recent version, this should be done) 
3. Run Extract Excel Data.Rmd
4. Run calf_production_hierarchical_stan.R. This script can run fixed_year and shifting_peak models. See Eguchi et al. (2026) for details. 
5. Edit 'gray whale calf production tech memo 2026 v4.qmd' to create a new report. This is a Quarto document. 
