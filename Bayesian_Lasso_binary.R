
# EBglmnet trial for Bayesian Lasso

library(EBglmnet)

heart <- read.csv("heart_train.csv")

heart_x <- as.matrix(heart[, c(1:11)])
heart_y <- heart[,13]

best_regularization_param <- cv.EBglmnet(heart_x, heart_y, family = "binomial",
                                         prior = "lasso", nfolds=3) 

covariates <- colnames(heart_x)

lasso_output <- EBglmnet(heart_x, heart_y, family="binomial",
                         prior = "lasso",
                         hyperparameters = 0.045)

selected_varnames <- covariates[lasso_output$fit[,1]]

final_result <- cbind(selected_varnames, lasso_output$fit[,c(3:6)])

saveRDS(list(intercept = lasso_output$Intercept, covariates = final_result), "Bayes_Lasso_Logistic.RDS")

