# Age-structured model for COVID-19 transmisison
Renata Retkute and Christopher Aidan Gilligan

Department of Plant Sciences, University of Cambridge, Downing Street, Cambridge
CB2 3EA, UK

This is an adaptation of age-structured model from Minter and Retkute 2019, Epidemics 29, 100368, "Approximate Bayesian Computation for infectious disease modelling" https://doi.org/10.1016/j.epidem.2019.100368

Parameter are estimated using the ABC-SMC algorithm with the tolerance calculated as the median from the previous generation.

Outbreak data is the COVID-19 epidemic in Hubei province in China as curated @: https://github.com/jriou/covid_adjusted_cfr

Figure 1. Comparing observed outbreak data with model simulations: number of daily cases (left) and percentage of cases in different age groups (right).
![](ABC_COVID19_results.png)


