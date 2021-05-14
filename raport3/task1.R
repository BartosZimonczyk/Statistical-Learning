library(glmnet)
library(xtable)
library(mvtnorm)
library(crayon)
library(SLOPE)

set.seed(42)

X <- matrix(rnorm(500*450, 0, 1/sqrt(500)), nrow=500, ncol=450)

# s <- min(eigen(sigma)$values)
# s <- min(2*s,1)
# sseq <- c(rep(s,p))
# V <- 2*diag(sseq)-diag(sseq)%*%solve(sigma)%*%diag(sseq)
# mu<-X-X%*%solve(sigma)%*%diag(sseq)

# knockoffs setup, used in all knockoff scenarios
s <- 1
sseq <- rep(s, 450)
V <- 2 * diag(sseq) - diag(sseq) %*% diag(sseq)
mu_knock <- X - X%*%diag(sseq)
Xn <- mu_knock + rmvnorm(500, rep(0, 450), V) / sqrt(500)
cX <- cbind(X, Xn)

# basic setup to perform simulations
ks <- c(5, 20, 50)
task <- 'task1'
rp1 <- 100

# TODO: for iii) - vii) estimate FDR and power  
# TODO: estimate MSE of b and mu for all apart iii)

final_table_fdr <- matrix(nrow=3, ncol=6)
final_table_power <- matrix(nrow=3, ncol=6)
final_table_mse_b <- matrix(nrow=3, ncol=7)
final_table_mse_mu <- matrix(nrow=3, ncol=7)

compute_fraction <- function(t, w){
  (1 + sum(w <= -t))/sum(w >= t)
}

# my own way to estimate threshold statistic
# it looks promising, but I am not sure
# moreover prof. Bogdan code is working as well,
# so I will stick to it

# power_general <- function(w, b, b_hat, k, q=0.2){
#   # ordered_w <- sort(abs(w), index.return=T)
#   fractions <- sapply(abs(w), compute_fraction, w=w)
#   TD <- 0
#   FD <- 0
#   FDR <- 0
#   if(length(which(fractions < q)) > 0){
#     treshold <- abs(w)[min(which(fractions < q))]
#     selected_variables <- which(w >= treshold)
#     if(length(selected_variables) > 0){
#       TD <- sum(selected_variables < k+1)
#       FD <- sum(selected_variables > k)
#       FDR <- FD/(TD+FD)
#     }
#   }
#   # power, FDR
#   return(c(TD/k, FDR))
# }

# totaly based on prof. Bogdan idea to calculate power and FDR
power_general <- function(u, b, b_hat, k, q=0.2){
  result <- sort(abs(u), decreasing = T, index.return = T)
  fd <- cumsum(u[result$ix] < 0)
  nd <- cumsum(u[result$ix] > 0)

  fdr <- (fd + 1)/nd

  fdp <- 0
  tp <- 0
  u1 <- which(fdr < q)
  if(length(u1) > 0){
    indopt <- max(u1)
    a1 <- result$ix[1:indopt]
    a2 <- which(u > 0)
    a3 <- intersect(a1, a2)
    new_b_hat <- rep(0, 450)
    new_b_hat[a3] <- b_hat[a3]
    tp <- sum(abs(b[a3]) > 0)
    fd <- length(a3) - tp
    if(length(a3) > 0){
      fdp <- fd/(fd+tp)
    }

  }
  return(c(tp/k, fdp))
}

for(j in 1:3){
    # creating new beta for each value of k
    k <- ks[j]
    b <- c(rep(10, k), rep(0, 450-k))
    
    # a) power and FDR
    # columns:
    # 1. knockoffs ridge
    # 2. knockoffs LASSO
    # 3. adaptive LASSO 1
    # 4. adaptive LASSO 2
    # 5. adaptive SLOPE
    # 6. adaptive Bayesian SLOPE
    fdr_matrix <- matrix(nrow=rp1, ncol=6)
    power_matrix <- matrix(nrow=rp1, ncol=6)
    
    # b) MSE of b and mu
    # columns:
    # 1. least squeares
    # 2. ridge
    # 3. LASSO
    # 4. adaptive LASSO 1
    # 5. adaptive LASSO 2
    # 6. adaptive SLOPE
    # 7. adaptive Bayesian SLOPE
    mse_b_matrix <- matrix(nrow=rp1, ncol=7)
    mse_mu_matrix <- matrix(nrow=rp1, ncol=7)
    cat(bold(blue('\nCurrently at k =', k, '\n')))
    for(i in 1:rp1){
        cat('\nCurrent iteration:', bold(i), '\n')
        eps <- 2*rnorm(500)
        Y <- X %*% b + eps
        
        # i)
        # TODO: least squares regression
        cat('OLS...\t\t\t\t')
        start_time_ols <- proc.time()[3]
        
        obj_ols <- lm(Y~X-1)
        b_hat_ols <- summary(obj_ols)$coefficients[,1]
        
        mse_b_matrix[i, 1] <- sum((b_hat_ols - b)^2)
        mse_mu_matrix[i, 1] <- sum((X %*% (b_hat_ols - b))^2)
        
        elapsed <- proc.time()[3] - start_time_ols
        cat(green('DONE\t'), round(elapsed, 3), '\n')
        
        # ii)
        # TODO: ridge with tuning parameter selected by CV
        cat('Ridge...\t\t\t')
        start_time_ridge <- proc.time()[3]
        
        obj_ridge <- cv.glmnet(X, Y, standarize=F, intercept=F, alpha=0)
        b_hat_ridge <- coef(obj_ridge, s='lambda.min')[2:451]
        
        mse_b_matrix[i, 2] <- sum((b_hat_ridge - b)^2)
        mse_mu_matrix[i, 2] <- sum((X %*% (b_hat_ridge - b))^2)
        
        elapsed <- proc.time()[3] - start_time_ridge
        cat(green('DONE\t'), round(elapsed, 3), '\n')
        
        # TODO: LASSO with tuning parameters selected by CV
        cat('LASSO...\t\t\t')
        start_time_lasso <- proc.time()[3]
        
        obj_lasso <- cv.glmnet(X, Y, standarize=F, intercept=F)
        b_hat_lasso <- coef(obj_lasso, s='lambda.min')[2:451]
        
        mse_b_matrix[i, 3] <- sum((b_hat_lasso - b)^2)
        mse_mu_matrix[i, 3] <- sum((X %*% (b_hat_lasso - b))^2)
        
        elapsed <- proc.time()[3] - start_time_lasso
        cat(green('DONE\t'), round(elapsed, 3), '\n')
        
        # iii)
        # TODO: knockoffs with ridge with FDR at level 0.2
        cat('Ridge knockoffs...\t\t')
        start_time_rk <- proc.time()[3]
        
        obj_ridge_knockoffs <- cv.glmnet(cX, Y, standarize=F, intercept=F, alpha=0)
        b_hat_ridge_knockoffs <- coef(obj_ridge_knockoffs, s='lambda.min')
        bhrk <- b_hat_ridge_knockoffs
        w_ridge <- abs(bhrk[2:451]) - abs(bhrk[452:901])
        
        res_stats_ridge <- power_general(w_ridge, b, b_hat_ridge, k, q=0.2)
        
        fdr_matrix[i, 1] <- res_stats_ridge[2]
        power_matrix[i, 1] <- res_stats_ridge[1]
        
        elapsed <- proc.time()[3] - start_time_rk
        cat(green('DONE\t'), round(elapsed, 3), '\n')
        
        # TODO: knockoffs with LASSO with FDR at level 0.2
        cat('LASSO knockoffs...\t\t')
        start_time_lk <- proc.time()[3]
        
        obj_lasso_knockoffs <- cv.glmnet(cX, Y, standarize=F, intercept=F)
        b_hat_lasso_knockoffs <- coef(obj_lasso_knockoffs, s='lambda.min')
        bhlk <- b_hat_lasso_knockoffs
        w_lasso <- abs(bhlk[2:451]) - abs(bhlk[452:901])
        
        res_stats_lasso <- power_general(w_lasso, b, b_hat_lasso, k, q=0.2)
        
        fdr_matrix[i, 2] <- res_stats_lasso[2]
        power_matrix[i, 2] <- res_stats_lasso[1]
        
        elapsed <- proc.time()[3] - start_time_lk
        cat(green('DONE\t'), round(elapsed, 3), '\n')
        
        # iv)
        # TODO: Adaptive LASSO I
        cat('Adaptive LASSO I...\t\t')
        start_time_al1 <- proc.time()[3]
        
        # for this case we will reuse LASSO setup
        # to reduce computation time
        ind_selected_lasso <- which(abs(b_hat_lasso) > 0)
        no_of_v_selected <- length(ind_selected_lasso)
        b_hat_selected_lasso <- b_hat_lasso[ind_selected_lasso]
        weights_alasso_1 <- 1/abs(b_hat_selected_lasso)
        X_alasso_1 <- X[, ind_selected_lasso]
        p_alasso_1 <- length(ind_selected_lasso)
        
        obj_alasso_1 <- cv.glmnet(X_alasso_1, Y, penalty.factor=weights_alasso_1, standarize=F, intercept=F)
        b_hat_alasso_1 <- rep(0, 450)
        b_hat_alasso_1[ind_selected_lasso] <- coef(obj_alasso_1, s='lambda.min')[2:(p_alasso_1+1)]
        ind_selected_alasso_1 <- which(abs(b_hat_alasso_1) > 0)
        
        fdr_matrix[i, 3] <- sum(ind_selected_alasso_1 > k)/max(1, length(ind_selected_alasso_1))
        power_matrix[i, 3] <- sum(ind_selected_alasso_1 <= k)/k
        mse_b_matrix[i, 4] <- sum((b_hat_alasso_1 - b)^2)
        mse_mu_matrix[i, 4] <- sum((X %*% (b_hat_alasso_1 - b))^2)
        
        elapsed <- proc.time()[3] - start_time_al1
        cat(green('DONE\t'), round(elapsed, 3), '\n')
        
        # v)
        # TODO: Adaptive LASSO II
        cat('Adaptive LASSO II...\t\t')
        start_time_al2 <- proc.time()[3]
        
        # for this case we will reuse LASSO and AdaLASSO1 setup
        # to reduce computation time
        rss_alasso_2 <- sum((Y - X%*%b_hat_lasso)^2)
        sigma_alasso_2 <- sqrt(rss_alasso_2/(500 - no_of_v_selected))
        weights_alasso_2 <- sigma_alasso_2/abs(b_hat_selected_lasso)
        X_alasso_2 <- X_alasso_1
        p_alasso_2 <- p_alasso_1
        
        obj_alasso_2 <- glmnet(
          X_alasso_2, Y, standardize=F, intercept=F, lambda=qnorm(1-0.2/2/450)*sigma_alasso_2/500,
          penalty.factor=weights_alasso_2
        )
        b_hat_alasso_2 <- rep(0, 450)
        b_hat_alasso_2[ind_selected_lasso] <- coef(obj_alasso_2)[2:(p_alasso_2+1)]
        ind_selected_alasso_2 <- which(abs(b_hat_alasso_2) > 0)
        
        fdr_matrix[i, 4] <- sum(ind_selected_alasso_2 > k)/max(1, length(ind_selected_alasso_2))
        power_matrix[i, 4] <- sum(ind_selected_alasso_2 <= k)/k
        mse_b_matrix[i, 5] <- sum((b_hat_alasso_2 - b)^2)
        mse_mu_matrix[i, 5] <- sum((X %*% (b_hat_alasso_2 - b))^2)
        
        elapsed <- proc.time()[3] - start_time_al2
        cat(green('DONE\t'), round(elapsed, 3), '\n')
        
        # vi)
        # TODO: Adaptive SLOPE, keep FDR at level 0.2
        cat('Adaptive SLOPE...\t\t')
        start_time_aslope <- proc.time()[3]
        
        # for this case we will reuse LASSO and AdaLASSOs setups
        # to reduce computations
        
        weights_aslope <- abs(b_hat_lasso)/sigma_alasso_2
        X_aslope <- sweep(X, 2, weights_aslope, '*')
        
        obj_aslope <- SLOPE(
          X_aslope, Y, q=0.2, lambda='bh', solver='admm',
          max_passes=100, scale = 'none', alpha = sigma_alasso_2/500
        )
        
        b_hat_aslope <- rep(0, 450)
        b_hat_aslope[ind_selected_lasso] <- coef(obj_aslope)[2:(p_alasso_1+1)] * weights_aslope
        ind_selected_aslope <- which(abs(b_hat_aslope) > 0)
        
        fdr_matrix[i, 5] <- sum(ind_selected_aslope > k)/max(1, length(ind_selected_aslope))
        power_matrix[i, 5] <- sum(ind_selected_aslope <= k)/k
        mse_b_matrix[i, 6] <- sum((b_hat_aslope - b)^2)
        mse_mu_matrix[i, 6] <- sum((X %*% (b_hat_aslope - b))^2)
        
        elapsed <- proc.time()[3] - start_time_aslope
        cat(green('DONE\t'), round(elapsed, 3), '\n')
        
        # vii)
        # For extra points +5, tempting
        # TODO: Adaptive Bayesian SLOPE with FDR at level 0.2
        cat('Adaptive Bayesian SLOPE...\t')
        start_time_abslope <- proc.time()[3]
        
        fdr_matrix[i, 6] <- NA_real_
        power_matrix[i, 6] <- NA_real_
        mse_b_matrix[i, 7] <- NA_real_
        mse_mu_matrix[i, 7] <- NA_real_
        
        elapsed <- proc.time()[3] - start_time_abslope
        cat(green('DONE\t'), round(elapsed, 3), '\n')
    }
    final_table_fdr[j, ] <- apply(fdr_matrix, 2, mean)
    final_table_power[j, ] <- apply(power_matrix, 2, mean)
    final_table_mse_b[j, ] <- apply(mse_b_matrix, 2, mean)
    final_table_mse_mu[j, ] <- apply(mse_mu_matrix, 2, mean)
}

df_fdr <- as.data.frame(
  final_table_fdr,
  row.names=c('k=5', 'k=20', 'k=50')
)

colnames(df_fdr) <- c(
  'ridge_knockoffs',
  'LASSO_knockoffs',
  'AdaLASSO1',
  'AdaLASSO2',
  'AdaSLOPE',
  'AdaBayesSLOPE'
)

df_power <- as.data.frame(
  final_table_power,
  row.names=c('k=5', 'k=20', 'k=50')
)

colnames(df_power) <- c(
  'ridge_knockoffs',
  'LASSO_knockoffs',
  'AdaLASSO1',
  'AdaLASSO2',
  'AdaSLOPE',
  'AdaBayesSLOPE'
)

df_mse_b <- as.data.frame(
  final_table_mse_b,
  row.names=c('k=5', 'k=20', 'k=50')
)

colnames(df_mse_b) <- c(
  'LS',
  'ridge',
  'LASSO',
  'AdaLASSO1',
  'AdaLASSO2',
  'AdaSLOPE',
  'AdaBayesSLOPE'
)

df_mse_mu <- as.data.frame(
  final_table_mse_mu,
  row.names=c('k=5', 'k=20', 'k=50')
)

colnames(df_mse_mu) <- c(
  'LS',
  'ridge',
  'LASSO',
  'AdaLASSO1',
  'AdaLASSO2',
  'AdaSLOPE',
  'AdaBayesSLOPE'
)

write.csv(df_fdr, file=paste(task, 'FDR.csv', sep=''))
write.csv(df_power, file=paste(task, 'power.csv', sep=''))
write.csv(df_mse_b, file=paste(task, 'MSE_b.csv', sep=''))
write.csv(df_mse_mu, file=paste(task, 'MSE_mu.csv', sep=''))

cat(bgGreen('Results wrtitten to .csv files!'), '\n')

##################
# Test section















