library(far)
library(glmnet)
library(xtable)

#############################################################
# task 1

# TODO by hand:
# 1. calculate lambda that minimize MSE in ridge regression for beta
# 2. calculate the bias, the variance and MSE of above optimal estimator 
# 3. Find critical value and calculate the power of the statistical
#    test based on the ridge regression and controlling FWER at level 0.1
#    REMINDER: controlling FWER at given level is equavilent to controlling
#              probability of type I error at the same level


final_table <- matrix(nrow=1, ncol=9)
k <- 20
rp1 <- 200
optimal_lambdas <- c()
ptm <- proc.time()
j <- 1

b <- c(rep(3.5, k), rep(0, 950-k))
b_hat_ridge <- matrix(nrow=950, ncol=rp1)
b_hat_ols <- matrix(nrow=950, ncol=rp1)
power_ridge_vector <- c()
power_ols_vector <- c()
for(i in 1:rp1){
  X <- matrix(rnorm(1000 * 950, 0, sqrt(1/1000)), nrow=1000)
  eps <- rnorm(1000)
  Y <- X %*% b + eps
  
  ridge_obj <- cv.glmnet(X, Y, standarize=F, intercept=F, alpha=0)
  ols_obj <- lm(Y ~ X-1)
  
  optimal_lambdas[i] <- ridge_obj$lambda.min
  
  
  b_hat_ridge[,i] <- coef(ridge_obj, s='lambda.min')[2:951]
  b_hat_ols[,i] <- summary(ols_obj)$coefficients[,1]
  
  # TODO: add power computation
  gamma_ridge <- optimal_lambdas[i]
  TD_ridge <- sum(abs(b_hat_ridge[1:k,i]) > 1/(1+gamma_ridge) * qnorm(1-0.1/2))
  power_ridge_vector[i] <- TD_ridge/k
  
  gamma_ols <- 0
  TD_ols <- sum(abs(b_hat_ols[1:k,i]) > 1/(1+gamma_ols) * qnorm(1-0.1/2))
  power_ols_vector[i] <- TD_ols/k
  
  cat('Iteration: ', i, 'for k = ', k, '\n')
  cat('Time elapsed: ', (proc.time()-ptm)[3], '\n')
}
mean_b_hat_ridge <- apply(b_hat_ridge, 1, mean)
mean_b_hat_ols <- apply(b_hat_ols, 1, mean)
var_b_hat_ridge <- apply(b_hat_ridge, 1, var)
var_b_hat_ols <- apply(b_hat_ols, 1, var)
bias_ridge <- sum((mean_b_hat_ridge - b)^2)
bias_ols <- sum((mean_b_hat_ols - b)^2)
var_ridge <- sum(var_b_hat_ridge)
var_ols <- sum(var_b_hat_ols)
mse_ridge <- var_ridge + bias_ridge
mse_ols <- var_ols + bias_ols
power_ridge <- mean(power_ridge_vector)
power_ols <- mean(power_ols_vector)

final_table[j,] <- c(
  bias_ols,
  bias_ridge,
  var_ols,
  var_ridge,
  mse_ols,
  mse_ridge,
  power_ols,
  power_ridge,
  mean(optimal_lambdas)
)


final_table <- as.data.frame(
  final_table,
  row.names = c('k=20', 'k=100', 'k=200')
)
colnames(final_table) <- c(
  'Bias_OLS',
  'Bias_ridge',
  'Variance_OLS',
  'Variance_ridge',
  'MSE_OLS',
  'MSE_ridge',
  'Power_OLS',
  'Power_ridge',
  'Mean_optimal_lambda'
)

file_name <- paste('r2result_task1_k', k, '.csv', sep='')
write.csv(final_table, file=file_name)

#############################################################
# task 2


























