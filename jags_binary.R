# logistic regression with slightly informative prior

library(R2jags)
library(foreach)
library(doParallel)

# set up parallel computing
cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset='2'))
cl <- makeCluster(cores)
registerDoParallel(cl)

rm(list=ls())

# get number of covariates
ncovar <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
ndim <- ncovar+1


heart_train <- read.csv("heart_train.csv")

train_full <- heart_train[,c(1:11, 13)]

# training covariates and binary outcomes
train_x <- as.matrix(heart_train[, c(1:11)])
y = heart_train[, 13]

covarnames <- colnames(train_x)

# all possible combinations of covariates given covariate numbers
covarchoice <- combn(11, ncovar)

# run all the models parallelly
models_summary <- foreach(i=1:ncol(covarchoice),
                          .packages = 'R2jags',
                          .inorder=FALSE) %dopar%{
                            
                            design_mat <- cbind(1, train_x[, covarchoice[, i]])
                            covars <- covarnames[covarchoice[,i]]
                            mean_vec <- rep(0, ndim)
                            precision_mat <- diag(1e-2, ndim)
                            
                            heart_data <- list('design_mat', 'y', 'mean_vec', 'precision_mat')
                            
                            heart_model <- function(){
                              for (j in 1:240){
                                logit[j] <- design_mat[j,] %*% beta
                                prob[j] <- exp(logit[j]) / (1 + exp(logit[j]))
                                y[j] ~ dbern(prob[j])
                              }
                              beta ~ dmnorm(mean_vec, precision_mat)
                            }
                            
                            save.parms <- c("beta")
                            
                            heart.out <- jags(data = heart_data, parameters.to.save = save.parms, model.file = heart_model,
                                              n.chains = 3, n.iter = 10000, n.burnin = 5000, n.thin = 1)
                            
                            beta_estimate <- heart.out$BUGSoutput$sims.list$beta
                            DIC <- heart.out$BUGSoutput$DIC
                            
                            list(beta = beta_estimate, DIC = DIC, covariates = covars)
                            
                          }


filename <- paste0("logistic_model_", ncovar, ".RDS")
saveRDS(models_summary, filename)




  

