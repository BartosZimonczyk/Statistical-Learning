library(far)
library(glmnet)
library(xtable)

set.seed(42)
#############################################################
# task 1

# TODO by hand:
# 1. calculate lambda that minimize MSE in ridge regression for beta
# 2. calculate the bias, the variance and MSE of above optimal estimator 
# 3. Find critical value and calculate the power of the statistical
#    test based on the ridge regression and controlling FWER at level 0.1


final_table <- matrix(nrow=1, ncol=9)
k <- 20
rp1 <- 200
optimal_lambdas <- c()
ptm <- proc.time()
j <- 1

cat('Computing orthonormal matrix...\n')
X <- matrix(rnorm(950*1000), nrow=1000, ncol=950)
X <- orthonormalization(X, FALSE, TRUE)
b <- c(rep(3.5, k), rep(0, 950-k))
b_hat_ridge <- matrix(nrow=rp1, ncol=950)
b_hat_ols <- matrix(nrow=rp1, ncol=950)
bias_ridge <- c()
bias_ols <- c()
mse_ridge_iter <- c()
mse_ols_iter <- c()
power_ridge_vector <- c()
power_ols_vector <- c()
for(i in 1:rp1){
  eps <- rnorm(1000)
  Y <- X %*% b + eps
  
  ridge_obj <- cv.glmnet(X, Y, standarize=F, intercept=F, alpha=0)
  ols_obj <- lm(Y ~ X-1)
  
  optimal_lambdas[i] <- ridge_obj$lambda.min
  
  b_hat_ridge[i,] <- coef(ridge_obj, s='lambda.min')[2:951]
  b_hat_ols[i,] <- summary(ols_obj)$coefficients[,1]
  
  bias_ridge[i] <- sum((b_hat_ridge[i,] - b)^2)
  bias_ols[i] <- sum((b_hat_ols[i,] - b)^2)
  
  # TODO: add power computation
  gamma_ridge <- optimal_lambdas[i]
  TD_ridge <- sum(abs(b_hat_ridge[i,1:k]) > 1/(1+gamma_ridge) * qnorm(1-0.1/2/950))
  power_ridge_vector[i] <- TD_ridge/k
  
  gamma_ols <- 0
  TD_ols <- sum(abs(b_hat_ols[i,1:k]) > 1/(1+gamma_ols) * qnorm(1-0.1/2/950))
  power_ols_vector[i] <- TD_ols/k
  
  mse_ridge_iter[i] <- sum((Y - X%*%b_hat_ridge[i,])^2)
  mse_ols_iter[i] <- sum((Y - X%*%b_hat_ols[i,])^2)
  
  cat('Iteration: ', i, 'for k = ', k, '\n')
  cat('Time elapsed: ', (proc.time()-ptm)[3], '\n')
}
var_ridge <- sum(apply(b_hat_ridge, 2, var))
var_ols <- sum(apply(b_hat_ols, 2, var))
bias_ridge <- mean(bias_ridge)
bias_ols <- mean(bias_ols)
mse_ridge <- mean(mse_ridge_iter)
mse_ols <- mean(mse_ols_iter)
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
# task 1 power computation


p <- seq(0, 1.5, by=0.01)
    
temp <- function(x){
  1 - pnorm(qnorm(1-0.1/2/950) - 3.5/(1+x)) + pnorm(-qnorm(1-0.1/2/950) - 3.5/(1+x))
}

plot(p, temp(p), type='l', lwd=2, ylab='Power', xlab=expression(gamma))

temp(p)













