library(MASS)
library(bigstep)
library(xtable)


# ex.1
set.seed(2137)

n <- 1000
p <- c(5, 10, 20, 100, 500, 950)
X <- matrix(rnorm(n*p[6], 0, sqrt(1/n)), nrow=n, ncol=p[6])
b <- c(rep(3, p[1]), rep(0, p[6]-p[1]))
eps <- rnorm(n)
Y <- X %*% b + eps

b_hat <- function(X, Y){solve(t(X) %*% X) %*% t(X) %*% Y}


my_model <- function(p, X, Y, b){
  X <- X[,1:p]
  b <- b[1:p]
  true_signal <- c(rep(T, 5), rep(F, p-5))
  
  b_h <- b_hat(X, Y)
  Y_h <- X %*% b_h
  M <- X %*% solve(t(X) %*% X) %*% t(X)
  RSS <- sum((Y-Y_h)^2)
  PE_1 <- sum((X%*%b - Y_h)^2) + n
  PE_2 <- RSS + 2*sum(diag(M))
  PE_3 <- RSS + 2*(RSS/(n-p))*sum(diag(M))
  PE_4 <- sum(((Y-Y_h)/(1-diag(M)))^2)

  FP_known <- sum(abs(b_h) > sqrt(2) & !true_signal)
  FN_known <- sum(!(abs(b_h) > sqrt(2)) & true_signal)
  FP_unknown <- sum(abs(b_h) > sqrt(RSS/(n-p)*2) & !true_signal)
  FN_unknown <- sum(!(abs(b_h) > sqrt(RSS/(n-p)*2)) & true_signal)
  AIC_known <- -n*log(sqrt(2*pi)) - RSS/2 - p
  AIC_unknown <- -n*log(sqrt(2*pi)) - n/2 * log(RSS) - p
  c(RSS, PE_1, PE_2, PE_3, PE_4, FP_known, FN_known, FP_unknown, FN_unknown, AIC_known, AIC_unknown)
}

final <- sapply(p, my_model, X, Y, b)
final <- as.data.frame(final, row.names=c('RSS',
                                 'PE1',
                                 'PE2',
                                 'PE3',
                                 'PE4',
                                 'FP_known',
                                 'FN_known',
                                 'FP_unknown',
                                 'FN_unknown',
                                 'AIC_known',
                                 'AIC_unknown'))
colnames(final) <- c('p=5', 'p=10', 'p=20', 'p=100', 'p=500', 'p=950')
final

set.seed(42)
n <- 1000
p <- c(5, 10, 20, 100, 500, 950)
results <- array(dim=c(100, 11, 6))
for(i in 1:100){
  X <- matrix(rnorm(n*p[6], 0, sqrt(1/n)), nrow=n, ncol=p[6])
  b <- c(rep(3, p[1]), rep(0, p[6]-p[1]))
  eps <- rnorm(n)
  Y <- X %*% b + eps
  results[i,,] <- sapply(p, my_model, X, Y, b)
  cat(i, '\n')
}

p_str <- c('p=5', 'p=10', 'p=20', 'p=100', 'p=500', 'p=950')

boxplot(results[,3,] - results[,2,],
        use.cols = T,
        names=c('p=5', 'p=10', 'p=20', 'p=100', 'p=500', 'p=950'),
        ylab='PE - PE1',
        main='Boxplots of difference between PE and it\'s first estimator')
boxplot(results[,4,] - results[,2,],
        use.cols = T,
        names=c('p=5', 'p=10', 'p=20', 'p=100', 'p=500', 'p=950'),
        ylab='PE - PE2',
        main='Boxplots of difference between PE and it\'s second estimator')
boxplot(results[,5,] - results[,2,],
        use.cols = T,
        names=c('p=5', 'p=10', 'p=20', 'p=100', 'p=500', 'p=950'),
        ylab='PE - PE3',
        main='Boxplots of difference between PE and it\'s third estimator')


xt <- xtable(final, digits=rep(3, 7))
print(xt, file="l1z1_1.txt")

# ex.2

library(bigstep)

ps <- c(20, 100, 500, 950)

ric <- function(loglik, k, p){2*k*log(p) - 2*loglik}

m_aic_ff <- c()
m_bic_ff <- c()
m_ric_ff <- c()
m_mbic_ff <- c()
m_mbic2_ff <- c()

set.seed(42)

# additional info: stepwise procedure is to slow to repeat it
# 100 times for each criterion, for each value of p

2*(1-pnorm(sqrt(log(1000))))

results_aic <- array(dim = c(100, 6, 4))
results_bic <- array(dim = c(100, 6, 4))
results_ric <- array(dim = c(100, 6, 4))
results_mbic <- array(dim = c(100, 6, 4))
results_mbic2 <- array(dim = c(100, 6, 4))

compute_stat <- function(big_model, b){
  s_big_model <- summary(big_model)
  b_hat <- rep(0, length(b))
  b_hat[as.numeric(big_model$model)] <- s_big_model$coefficients[-1, 1]
  true_signal <- c(rep(T, 5), rep(F, length(b)-5))
  FD <- sum(as.numeric(big_model$model) > 5)
  TD <- sum(as.numeric(big_model$model) <= 5)
  SEy <- mean((big_model$y - big_model$Xm %*% s_big_model$coefficients[,1])^2)
  SEb <- mean((b - b_hat)^2)
  SEex <- mean((big_model$X %*% b - big_model$Xm %*% s_big_model$coefficients[,1])^2)
  pow <- TD/5
  c(FD, TD, pow, SEy, SEb, SEex)
}
  
ps <- c(20, 100, 500, 950)

for(i in 1:100){
  n <- 1000
  p <- c(5, 10, 20, 100, 500, 950)
  X <- matrix(rnorm(n*p[6], 0, sqrt(1/n)), nrow=n, ncol=p[6])
  b <- c(rep(3, p[1]), rep(0, p[6]-p[1]))
  eps <- rnorm(n)
  Y <- X %*% b + eps
  
  cat('Iteration:', i, '\n')
  for(j in 1:4){
    cat('\tp=', ps[j], '\n')
    d <- prepare_data(Y, X[,1:ps[j]], verbose = F)
    
    cat('\t\tAIC\n')
    m_aic_ff <- fast_forward(d, crit=aic)
    results_aic[i,,j] <- compute_stat(m_aic_ff, b[1:ps[j]])
    
    cat('\t\tBIC\n')
    m_bic_ff <- fast_forward(d, crit=bic)
    results_bic[i,,j] <- compute_stat(m_bic_ff, b[1:ps[j]])
    
    cat('\t\tRIC\n')
    m_ric_ff <- fast_forward(d, crit=ric)
    results_ric[i,,j] <- compute_stat(m_ric_ff, b[1:ps[j]])
    
    cat('\t\tmBIC\n')
    m_mbic_ff <- fast_forward(d, crit=mbic)
    results_mbic[i,,j] <- compute_stat(m_mbic_ff, b[1:ps[j]])
    
    cat('\t\tmBIC2\n')
    m_mbic2_ff <- fast_forward(d, crit=mbic2)
    results_mbic2[i,,j] <- compute_stat(m_mbic2_ff, b[1:ps[j]])
  }
}


# there are models for which some criterions choose 0 variables, in this case we define
# FDR equal to 0

FDR <- function(result){
  fdr_ps <- apply(result, 2, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2])))
  fdr_ps
}

mean_FDR_aic <- apply(apply(results_aic, 1, FDR), 1, mean)
mean_FDR_bic <- apply(apply(results_bic, 1, FDR), 1, mean)
mean_FDR_ric <- apply(apply(results_ric, 1, FDR), 1, mean)
mean_FDR_mbic <- apply(apply(results_mbic, 1, FDR), 1, mean)
mean_FDR_mbic2 <- apply(apply(results_mbic2, 1, FDR), 1, mean)

mean_FDR_aic
mean_FDR_bic
mean_FDR_ric
mean_FDR_mbic
mean_FDR_mbic2

mean_power_aic <- apply(results_aic, c(2,3), mean)[3,]
mean_power_bic <- apply(results_bic, c(2,3), mean)[3,]
mean_power_ric <- apply(results_ric, c(2,3), mean)[3,]
mean_power_mbic <- apply(results_mbic, c(2,3), mean)[3,]
mean_power_mbic2 <- apply(results_mbic2, c(2,3), mean)[3,]

mean_power_aic
mean_power_bic
mean_power_ric
mean_power_mbic
mean_power_mbic2

mean_MSEy_aic <- apply(results_aic, c(2,3), mean)[4,]
mean_MSEy_bic <- apply(results_bic, c(2,3), mean)[4,]
mean_MSEy_ric <- apply(results_ric, c(2,3), mean)[4,]
mean_MSEy_mbic <- apply(results_mbic, c(2,3), mean)[4,]
mean_MSEy_mbic2 <- apply(results_mbic2, c(2,3), mean)[4,]

mean_SE_aic
mean_SE_bic
mean_SE_ric
mean_SE_mbic
mean_SE_mbic2

mean_FD_aic <- apply(results_aic, c(2,3), mean)[1,]
mean_FD_bic <- apply(results_bic, c(2,3), mean)[1,]
mean_FD_ric <- apply(results_ric, c(2,3), mean)[1,]
mean_FD_mbic <- apply(results_mbic, c(2,3), mean)[1,]
mean_FD_mbic2 <- apply(results_mbic2, c(2,3), mean)[1,]

mean_TD_aic <- apply(results_aic, c(2,3), mean)[2,]
mean_TD_bic <- apply(results_bic, c(2,3), mean)[2,]
mean_TD_ric <- apply(results_ric, c(2,3), mean)[2,]
mean_TD_mbic <- apply(results_mbic, c(2,3), mean)[2,]
mean_TD_mbic2 <- apply(results_mbic2, c(2,3), mean)[2,]

mean_MSEb_aic <- apply(results_aic, c(2,3), mean)[5,]
mean_MSEb_bic <- apply(results_bic, c(2,3), mean)[5,]
mean_MSEb_ric <- apply(results_ric, c(2,3), mean)[5,]
mean_MSEb_mbic <- apply(results_mbic, c(2,3), mean)[5,]
mean_MSEb_mbic2 <- apply(results_mbic2, c(2,3), mean)[5,]

mean_MSEex_aic <- apply(results_aic, c(2,3), mean)[6,]
mean_MSEex_bic <- apply(results_bic, c(2,3), mean)[6,]
mean_MSEex_ric <- apply(results_ric, c(2,3), mean)[6,]
mean_MSEex_mbic <- apply(results_mbic, c(2,3), mean)[6,]
mean_MSEex_mbic2 <- apply(results_mbic2, c(2,3), mean)[6,]


library(xtable)

xt <- xtable(rbind(mean_FD_aic,
             mean_FD_bic,
             mean_FD_ric,
             mean_FD_mbic,
             mean_FD_mbic2), digits = rep(3,5))

print(xt, file="l1z2_1.txt")

xt <- xtable(rbind(mean_TD_aic,
             mean_TD_bic,
             mean_TD_ric,
             mean_TD_mbic,
             mean_TD_mbic2), digits = rep(3,5))

print(xt, file="l1z2_2.txt")

xt <- xtable(rbind(mean_FDR_aic,
                mean_FDR_bic,
                mean_FDR_ric,
                mean_FDR_mbic,
                mean_FDR_mbic2), digits = rep(3,5))

print(xt, file="l1z2_3.txt")

xt <- xtable(rbind(mean_power_aic,
                mean_power_bic,
                mean_power_ric,
                mean_power_mbic,
                mean_power_mbic2), digits = rep(3,5))

print(xt, file="l1z2_4.txt")

xt <- xtable(rbind(mean_MSEy_aic,
                mean_MSEy_bic,
                mean_MSEy_ric,
                mean_MSEy_mbic,
                mean_MSEy_mbic2), digits = rep(3,5))

print(xt, file="l1z2_5.txt")

xt <- xtable(rbind(mean_MSEb_aic,
             mean_MSEb_bic,
             mean_MSEb_ric,
             mean_MSEb_mbic,
             mean_MSEb_mbic2), digits = rep(3,5))

print(xt, file="l1z2_6.txt")

xt <- xtable(rbind(mean_MSEex_aic,
                   mean_MSEex_bic,
                   mean_MSEex_ric,
                   mean_MSEex_mbic,
                   mean_MSEex_mbic2), digits = rep(3,5))

print(xt, file="l1z2_7.txt")


# ex.3

library(bigstep)

ric <- function(loglik, k, n, p){mbic(loglik, k, n, p, const=exp(log(n)/2))}

set.seed(42)

# additional info: stepwise procedure is to slow to repeat it
# 100 times for each criterion, for each value of p

results_ric <- array(dim = c(100, 5))
results_mbic <- array(dim = c(100, 5))
results_mbic2 <- array(dim = c(100, 5))

compute_stat <- function(big_model, b){
  s_big_model <- summary(big_model)
  b_hat <- rep(0, length(b))
  b_hat[as.numeric(big_model$model)] <- s_big_model$coefficients[-1, 1]
  true_signal <- c(rep(T, 50), rep(F, length(b)-50))
  FD <- sum(as.numeric(big_model$model) > 50)
  TD <- sum(as.numeric(big_model$model) <= 50)
  SEy <- mean((big_model$y - big_model$Xm %*% s_big_model$coefficients[,1])^2)
  SEb <- mean((b - b_hat)^2)
  pow <- TD/50
  c(FD, TD, pow, SEy, SEb)
}


for(i in 1:100){
  n <- 1000
  X <- matrix(rnorm(n*950, 0, sqrt(1/n)), nrow=n, ncol=950)
  b <- c(rep(3, 50), rep(0, 900))
  eps <- rnorm(n)
  Y <- X %*% b + eps
  
  cat('Iteration:', i, '\n')
  d <- prepare_data(Y, X, verbose = F)
  
  cat('\t\tRIC\n')
  m_ric_ff <- fast_forward(d, crit=ric)
  results_ric[i,] <- compute_stat(m_ric_ff, b)
  
  cat('\t\tmBIC\n')
  m_mbic_ff <- fast_forward(d, crit=mbic)
  results_mbic[i,] <- compute_stat(m_mbic_ff, b)
  
  cat('\t\tmBIC2\n')
  m_mbic2_ff <- fast_forward(d, crit=mbic2)
  results_mbic2[i,] <- compute_stat(m_mbic2_ff, b)
}


mean_FDR_ric <- mean(apply(results_ric, 1, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2]))))
mean_FDR_mbic <- mean(apply(results_mbic, 1, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2]))))
mean_FDR_mbic2 <- mean(apply(results_mbic2, 1, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2]))))

mean_FDR_ric
mean_FDR_mbic
mean_FDR_mbic2

mean_ric <- apply(results_ric, 2, mean)
mean_mbic <- apply(results_mbic, 2, mean)
mean_mbic2 <- apply(results_mbic2, 2, mean)

mean_ric
mean_mbic
mean_mbic2

df <- data.frame(cbind(c(mean_FDR_ric,
                         mean_FDR_mbic,
                         mean_FDR_mbic2),
                       rbind(mean_ric,
                             mean_mbic,
                             mean_mbic2)))
colnames(df) <- c('FDR', 'FD', 'TD', 'power', 'MSEy', 'MSEb')
xt <- xtable(df, digits=rep(3,7))

print(xt, file="l1z3_1.txt")


# ex.4

rbic <- function(loglik, k, n, p){mbic(loglik, k, n, p)}
rbic2 <- function(loglik, k, n, p){mbic2(loglik, k, n, p)}

compute_stat <- function(big_model, b, Y){
  s_big_model <- summary(big_model)
  
  linear_model <- lm(Y ~ big_model$Xm - 1)
  robust_model_huber <- rlm(Y ~ big_model$Xm - 1, psi=psi.huber)
  robust_model_bisq <- rlm(Y ~ big_model$Xm - 1, psi=psi.bisquare)
  
  b_hat_ls <- rep(0, length(b))
  b_hat_ls[as.numeric(big_model$model)] <- summary(linear_model)$coefficients[-1, 1]
  b_hat_rrh <- rep(0, length(b))
  b_hat_rrh[as.numeric(big_model$model)] <- summary(robust_model_huber)$coefficients[-1, 1]
  b_hat_rrb <- rep(0, length(b))
  b_hat_rrb[as.numeric(big_model$model)] <- summary(robust_model_bisq)$coefficients[-1, 1]
  
  
  true_signal <- c(rep(T, 30), rep(F, length(b)-30))
  FD <- sum(as.numeric(big_model$model) > 30)
  TD <- sum(as.numeric(big_model$model) <= 30)
  SEy_ls <- mean((Y - big_model$Xm %*% summary(linear_model)$coefficients[,1])^2)
  SEy_rrh <- mean((Y - big_model$Xm %*% summary(robust_model_huber)$coefficients[, 1])^2)
  SEy_rrb <- mean((Y - big_model$Xm %*% summary(robust_model_bisq)$coefficients[, 1])^2)
  SEb_ls <- mean((b - b_hat_ls)^2)
  SEb_rrh <- mean((b - b_hat_rrh)^2)
  SEb_rrb <- mean((b - b_hat_rrb)^2)
  pow <- TD/30
  c(FD, TD, pow, SEy_ls, SEy_rrh, SEy_rrb, SEb_ls, SEb_rrh, SEb_rrb)
}

#exponential errors
set.seed(42)

results_mbic <- array(dim = c(100, 9))
results_mbic2 <- array(dim = c(100, 9))
results_rbic <- array(dim = c(100, 9))
results_rbic2 <- array(dim = c(100, 9))

for(i in 1:100){
  n <- 1000
  X <- matrix(rnorm(n*950, 0, sqrt(1/n)), nrow=n, ncol=950)
  b <- c(rep(10, 30), rep(0, 920))
  eps <- rexp(n) - 1
  Y <- X %*% b + eps
  rY <- rank(Y)
  
  d <- prepare_data(Y, X, verbose = F)
  rd <- prepare_data(rY, X, verbose = F)
  
  cat('Iteration:', i, '\n\n')
  cat('mBIC\n')
  mbic_ff <- fast_forward(d, crit=mbic)
  results_mbic[i,] <- compute_stat(mbic_ff, b, Y)
  
  cat('mBIC2\n')
  mbic2_ff <- fast_forward(d, crit=mbic2)
  results_mbic2[i,] <- compute_stat(mbic2_ff, b, Y)
  
  cat('rBIC\n')
  rbic_ff <- fast_forward(rd, crit=rbic)
  results_rbic[i,] <- compute_stat(rbic_ff, b, Y)
  
  cat('rBIC2\n')
  rbic2_ff <- fast_forward(rd, crit=rbic2)
  results_rbic2[i,] <- compute_stat(rbic2_ff, b, Y)
}

mean_FDR_mbic <- mean(apply(results_mbic, 1, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2]))))
mean_FDR_mbic2 <- mean(apply(results_mbic2, 1, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2]))))
mean_FDR_rbic <- mean(apply(results_rbic, 1, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2]))))
mean_FDR_rbic2 <- mean(apply(results_rbic2, 1, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2]))))

mean_mbic <- apply(results_mbic, 2, mean)
mean_mbic2 <- apply(results_mbic2, 2, mean)
mean_rbic <- apply(results_rbic, 2, mean)
mean_rbic2 <- apply(results_rbic2, 2, mean)

df <- data.frame(cbind(c(mean_FDR_mbic,
                         mean_FDR_mbic2,
                         mean_FDR_rbic,
                         mean_FDR_rbic2),
                       rbind(mean_mbic[1:6],
                             mean_mbic2[1:6],
                             mean_rbic[1:6],
                             mean_rbic2[1:6])))
colnames(df) <- c('FDR', 'FD', 'TD', 'power', 'MSEy_ls', 'MSEy_rrh', 'MSEy_rrb')
xt <- xtable(df, digits=rep(3,8))

print(xt, file="l1z4_1.txt")

df <- data.frame(rbind(mean_mbic[7:9],
                       mean_mbic2[7:9],
                       mean_rbic[7:9],
                       mean_rbic2[7:9]))
colnames(df) <- c('MSEb_ls', 'MSEb_rrh', 'MSEb_rrb')
xt <- xtable(df, digits=rep(3,4))
print(xt, file="l1z4_2.txt")


# cauchy distributed errors

#cauchy errors
set.seed(42)

compute_stat <- function(big_model, b, Y){
  s_big_model <- summary(big_model)
  
  linear_model <- lm(Y ~ big_model$Xm - 1)
  robust_model_huber <- rlm(Y ~ big_model$Xm - 1, psi=psi.huber)
  robust_model_bisq <- rlm(Y ~ big_model$Xm - 1, psi=psi.bisquare)
  
  b_hat_ls <- rep(0, length(b))
  b_hat_ls[as.numeric(big_model$model)] <- summary(linear_model)$coefficients[-1, 1]
  b_hat_rrh <- rep(0, length(b))
  b_hat_rrh[as.numeric(big_model$model)] <- summary(robust_model_huber)$coefficients[-1, 1]
  b_hat_rrb <- rep(0, length(b))
  b_hat_rrb[as.numeric(big_model$model)] <- summary(robust_model_bisq)$coefficients[-1, 1]
  
  
  true_signal <- c(rep(T, 30), rep(F, length(b)-30))
  FD <- sum(as.numeric(big_model$model) > 30)
  TD <- sum(as.numeric(big_model$model) <= 30)
  SEy_ls <- mean((Y - big_model$Xm %*% summary(linear_model)$coefficients[,1])^2)
  SEy_rrh <- mean((Y - big_model$Xm %*% summary(robust_model_huber)$coefficients[, 1])^2)
  SEy_rrb <- mean((Y - big_model$Xm %*% summary(robust_model_bisq)$coefficients[, 1])^2)
  SEb_ls <- mean((b - b_hat_ls)^2)
  SEb_rrh <- mean((b - b_hat_rrh)^2)
  SEb_rrb <- mean((b - b_hat_rrb)^2)
  pow <- TD/30
  c(FD, TD, pow, SEy_ls, SEy_rrh, SEy_rrb, SEb_ls, SEb_rrh, SEb_rrb)
}

results_mbic <- array(dim = c(100, 9))
results_mbic2 <- array(dim = c(100, 9))
results_rbic <- array(dim = c(100, 9))
results_rbic2 <- array(dim = c(100, 9))

for(i in 1:100){
  n <- 1000
  X <- matrix(rnorm(n*950, 0, sqrt(1/n)), nrow=n, ncol=950)
  b <- c(rep(10, 30), rep(0, 920))
  eps <- rcauchy(n)
  Y <- X %*% b + eps
  rY <- rank(Y)
  
  d <- prepare_data(Y, X, verbose = F)
  rd <- prepare_data(rY, X, verbose = F)
  
  cat('Iteration:', i, '\n\n')
  cat('mBIC\n')
  mbic_ff <- fast_forward(d, crit=mbic)
  results_mbic[i,] <- compute_stat(mbic_ff, b, Y)
  
  cat('mBIC2\n')
  mbic2_ff <- fast_forward(d, crit=mbic2)
  results_mbic2[i,] <- compute_stat(mbic2_ff, b, Y)
  
  cat('rBIC\n')
  rbic_ff <- fast_forward(rd, crit=rbic)
  results_rbic[i,] <- compute_stat(rbic_ff, b, Y)
  
  cat('rBIC2\n')
  rbic2_ff <- fast_forward(rd, crit=rbic2)
  results_rbic2[i,] <- compute_stat(rbic2_ff, b, Y)
}

mean_FDR_mbic <- mean(apply(results_mbic, 1, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2]))))
mean_FDR_mbic2 <- mean(apply(results_mbic2, 1, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2]))))
mean_FDR_rbic <- mean(apply(results_rbic, 1, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2]))))
mean_FDR_rbic2 <- mean(apply(results_rbic2, 1, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2]))))

mean_mbic <- apply(results_mbic, 2, mean)
mean_mbic2 <- apply(results_mbic2, 2, mean)
mean_rbic <- apply(results_rbic, 2, mean)
mean_rbic2 <- apply(results_rbic2, 2, mean)

df <- data.frame(cbind(c(mean_FDR_mbic,
                         mean_FDR_mbic2,
                         mean_FDR_rbic,
                         mean_FDR_rbic2),
                       rbind(mean_mbic[1:6],
                             mean_mbic2[1:6],
                             mean_rbic[1:6],
                             mean_rbic2[1:6])))
colnames(df) <- c('FDR', 'FD', 'TD', 'power', 'MSEy_ls', 'MSEy_rrh', 'MSEy_rrb')
xt <- xtable(df, digits=rep(3,8))

print(xt, file="l1z4_3.txt")

df <- data.frame(rbind(mean_mbic[7:9],
                       mean_mbic2[7:9],
                       mean_rbic[7:9],
                       mean_rbic2[7:9]))
colnames(df) <- c('MSEb_ls', 'MSEb_rrh', 'MSEb_rrb')
xt <- xtable(df, digits=rep(3,4))

print(xt, file="l1z4_4.txt")

rbic2_ff

# ex.5

set.seed(42)
sigmoid <- function(x){1/(1+exp(-x))}

compute_stat <- function(big_model, b, Y){
  s_big_model <- summary(big_model)
  b_hat <- rep(0, length(b))
  b_hat[as.numeric(big_model$model)] <- s_big_model$coefficients[-1, 1]
  true_signal <- c(rep(T, 30), rep(F, length(b)-30))
  FD <- sum(as.numeric(big_model$model) > 30)
  TD <- sum(as.numeric(big_model$model) <= 30)
  SEy <- mean((Y - ifelse(sigmoid(big_model$Xm %*% s_big_model$coefficients[,1]) > 0.5, 1, 0))^2)
  pow <- TD/30
  SEb <- mean((b - b_hat)^2)
  c(FD, TD, pow, SEy, SEb)
}

results_aic <- array(dim = c(100, 5))
results_bic <- array(dim = c(100, 5))
results_mbic <- array(dim = c(100, 5))
results_mbic2 <- array(dim = c(100, 5))

for(i in 1:100){
  n <- 1000
  X <- matrix(rnorm(n*950, 0, sqrt(1/n)), nrow=n, ncol=950)
  b <- c(rep(10, 30), rep(0, 920))
  p <- sigmoid(X %*% b)
  Y <- rbinom(n, 1, prob=p)
  
  d <- prepare_data(Y, X, verbose = F, type='logistic')
  
  cat('\nIteration:', i, '\n\n')
  
  cat('AIC\n')
  aic_ff <- fast_forward(rd, crit=aic)
  results_aic[i,] <- compute_stat(aic_ff, b, Y)
  
  cat('BIC\n')
  bic_ff <- fast_forward(rd, crit=bic)
  results_bic[i,] <- compute_stat(bic_ff, b, Y)
  
  cat('mBIC\n')
  mbic_ff <- fast_forward(d, crit=mbic)
  results_mbic[i,] <- compute_stat(mbic_ff, b, Y)
  
  cat('mBIC2\n')
  mbic2_ff <- fast_forward(d, crit=mbic2)
  results_mbic2[i,] <- compute_stat(mbic2_ff, b, Y)
}

mean_FDR_aic <- mean(apply(results_aic, 1, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2]))))
mean_FDR_bic <- mean(apply(results_bic, 1, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2]))))
mean_FDR_mbic <- mean(apply(results_mbic, 1, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2]))))
mean_FDR_mbic2 <- mean(apply(results_mbic2, 1, function(x) ifelse(x[1]+x[2] == 0, 0, x[1]/(x[1] + x[2]))))

mean_FDR_aic
mean_FDR_bic
mean_FDR_mbic
mean_FDR_mbic2

mean_aic <- apply(results_aic, 2, mean)
mean_bic <- apply(results_bic, 2, mean)
mean_mbic <- apply(results_mbic, 2, mean)
mean_mbic2 <- apply(results_mbic2, 2, mean)

mean_aic
mean_bic
mean_mbic
mean_mbic2

df <- data.frame(cbind(c(mean_FDR_aic,
                         mean_FDR_bic,
                         mean_FDR_mbic,
                         mean_FDR_mbic2),
                       rbind(mean_aic,
                             mean_bic,
                             mean_mbic,
                             mean_mbic2)))
colnames(df) <- c('FDR', 'FD', 'TD', 'power', 'MSEy', 'MSEb')
xt <- xtable(df, digits=rep(3,7))

print(xt, file="l1z5_1.txt")
  








