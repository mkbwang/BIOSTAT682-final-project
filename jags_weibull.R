
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
ncovar <- 1
ndim <- ncovar+1

mean_vec <- rep(0, ndim)
precision_mat <- diag(ndim) * 1e-2

heart <- read.csv("heart_train.csv")

# set up the follow up time and censorship information
follow_time <- heart$time
indicator <- 1-heart$DEATH_EVENT
censor_time <- follow_time
censor_time[indicator == 0] <- 360 
follow_time[indicator == 1] <- NA

train_x <- heart[, c(1:11)]

covarlist <- colnames(train_x)

covarchoice <- combn(11, ncovar)


models_summary <- foreach(i=1:ncol(covarchoice),
        .packages = 'R2jags',
        .inorder=FALSE) %dopar% {
          
          covariates <- cbind(1, as.matrix(train_x[, covarchoice[, i]]))
          covarnames <- covarlist[covarchoice[, i]]
          
          heart_init <- list(list("beta" = rep(0, ndim), "v" = 1),
                             list("beta" = rep(0, ndim), "v" = 1),
                             list("beta" = rep(0, ndim), "v" = 1))
          
          heart_data <- list("covariates", "mean_vec", "precision_mat",
                             "follow_time", "indicator", "censor_time")
          
          heart_surv_model <- function(){
            for (i in 1:240){
              lambda[i] <- (exp(-covariates[i,] %*% beta))^v
              follow_time[i] ~ dweib(v, lambda[i])
              indicator[i] ~ dinterval(follow_time[i], censor_time[i])
            }
            v ~ dunif(0,5)
            sigma <- 1/v
            beta ~ dmnorm(mean_vec, precision_mat)
          }
          
          save.parms <- c("beta", "sigma")
          
          output <- tryCatch({
            heart.out <- jags(data = heart_data, inits = heart_init, parameters.to.save = save.parms, model.file = heart_surv_model,
                              n.chains = 3, n.iter = 10000, n.burnin = 5000, n.thin = 1)
            
            beta_estimate <- heart.out$BUGSoutput$sims.list$beta
            sigma_estimate <- heart.out$BUGSoutput$sims.list$sigma
            DIC <- heart.out$BUGSoutput$DIC
            
            return(list(covariates = covarnames, beta = beta_estimate, sigma = sigma_estimate, DIC=DIC))
          }, 
          error = function(cond){
            return(list(failureid = i))
          },
          warning = function(cond){
            return(list(warningid = i))
          })
          output
        }


filename <- paste0("weibull_model_", ncovar, ".RDS")
saveRDS(models_summary, filename)


