# BIOSTAT 682 Final Project
This repository contains part of the code and training/test data that can illustrate how the Bayesian analysis was done for this project. The raw data was downloaded from [Kaggle](https://www.kaggle.com/andrewmvd/heart-failure-clinical-data).

* `jags_binary.R` contains the code for best subset selection based on DIC for logistic regression.
* `Bayesian_Lasso_binary.R` contains the code for lasso logistic regression in Bayesian way.
* `jags_weibull.R` contains the code for best subset selection based on DIC for survival analysis based on Weibull model.
* `jags_binary.slurm` is a script used to submit job arrays. It uses `jags_binary.R`.
* `jags_weibull.slurm` is a script used to submit job arrays. It uses `jags_weibull.R`.