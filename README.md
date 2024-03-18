# GrayWhaleCalfProduction
This repository contains R code to extract and analyze gray whale calf production data that are collected at Piedras Blancas, CA, annually. 
Necessary functions are stored in GrayWhaleCalfProduction_fcns.R. Some of the functions in the file may be obsolete and unused. The analytical
part is conducted via JAGS (Just Another Gibbs Sampler, which can be downloaded from https://sourceforge.net/projects/mcmc-jags/files/. 

It is assumed that this repository is donwloaded to your local computer via RStudio and a project created. Please create three empty folders
within the project folder and name them "data," "figures," and "RData." 

In 2023, new data extraction code was developed due to some inconsistencies that were found in the previous version (Collating and Formatting Passdown.R).
The comparison of the two extraction methods is provided in a notebook file named "compare_data_extraction.Rmd." In the past report, I did not 
use the new extraction method. I appended the new estimate for the 2023 season to the 2022 report. For the 2023 report, I did use the new extraction but
the table was replaced with the old version. These reports can be found in Report_calf_production_2022.Rmd and Report_calf_production_2023.Rmd.

The order of the data processing and analysis should be:
1. Raw data should be saved in the "data/All data" directory. Foloder names are C_C YYYY, where YYYY indicates a four-digit year. 
2. Extract Excel Data.Rmd
3. Run calf_production_jags.R
4. Edit Report_calf_production_2023.Rmd to create a report.

