library(pesel)

all_p <- c(100, 500, 1000, 5000)

n <- 50
k <- 5
for(p in all_p){
    cat('Currently working with p=', p, '\n')
    dims_normal <- c()
    dims_exp <- c()
    dims_cauchy <- c()
    for(i in 1:100){
        cat('Iteration: ', i, '\n')
        F_matrix <- matrix(rnorm(n * k), nrow=n, ncol=k)
        W_matrix <- matrix(rnorm(p * k), nrow=p, ncol=k)
        W_matrix <- t(t(W_matrix) * 5:1)
        
        E_matrix_normal <- matrix(rnorm(n * p, 0, 10), nrow=n, ncol=p)
        E_matrix_exp <- matrix(rexp(n * p, 1/10), nrow=n, ncol=p)
        E_matrix_cauchy <- matrix(rcauchy(n * p), nrow=n, ncol=p)
        
        X <- F_matrix %*% t(W_matrix)
        X_normal <- X + E_matrix_normal
        X_exp <- X + E_matrix_exp
        X_cauchy <- X + E_matrix_cauchy
        
        dims_normal[i] <- pesel(X_normal)$nPCs
        dims_exp[i] <- pesel(X_exp)$nPCs
        dims_cauchy[i] <- pesel(X_cauchy)$nPCs
    }
    png(paste('task2_hist_normal_', p, '.png', sep=''))
    hist(dims_normal, xlab='Number of PCs selected by PESEL', ylab='Frequency',
         main=paste('For normal errors and p =', p), breaks=0:6)
    dev.off()
    
    png(paste('task2_hist_exp_', p, '.png', sep=''))
    hist(dims_exp, xlab='Number of PCs selected by PESEL', ylab='Frequency',
         main=paste('For exp errors and p =', p), breaks=0:6)
    dev.off()
    
    png(paste('task2_hist_cauchy_', p, '.png', sep=''))
    hist(dims_cauchy, xlab='Number of PCs selected by PESEL', ylab='Frequency',
         main=paste('For cauchy errors and p =', p), breaks=0:6)
    dev.off()
}