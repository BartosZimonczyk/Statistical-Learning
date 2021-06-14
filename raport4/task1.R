set.seed(42)
all_n <- c(500, 2000)

results <- matrix(nrow=8, ncol=2)
for(n in all_n){
    cat('Doing n =', n, '...\n')
    corr_matrix<- matrix(ncol=100, nrow=6)
    sigma_mse <- matrix(ncol=100, nrow=2)
    for(i in 1:100){
        cat('Iteration: ', i, '\n')
        F_matrix <- matrix(rnorm(n*3), nrow=n, ncol=3)
        W_matrix <- matrix(rnorm(100 * 3), nrow=100, ncol=3)
        W_matrix <- t(t(W_matrix) * 3:1)
        E_matrix <- matrix(rnorm(500 * 100, 0, 1), nrow=n, ncol=100)
        X <- F_matrix %*% t(W_matrix) + E_matrix
        
        # a)
        # eigenvalue decomposition of XTX
        e_out <- eigen(t(X) %*% X)
        e_values <- e_out$values
        e_vectors <- e_out$vectors
        e_v_sorted_idx <- sort(abs(e_values), decreasing = T, index.return = T)$ix
        weights_matrix <- e_vectors[,e_v_sorted_idx]
        PC_eigen <- X %*% weights_matrix
        corr_matrix[1:3, i] <- c(
            abs(cor(F_matrix[,1], PC_eigen[,1])),
            abs(cor(F_matrix[,2], PC_eigen[,2])),
            abs(cor(F_matrix[,3], PC_eigen[,3]))
        )
        sigma_mse[1, i] <- (1/(n * (100-3)) * sum(e_values[4:100]) - 1) ^ 2
        
        
        # b)
        # SVD of X
        
        X.svd <- svd(X)
        U <- X.svd$u
        D <- diag(X.svd$d)
        PC_svd <- U %*% D
        corr_matrix[4:6, i] <- c(
            abs(cor(F_matrix[,1], PC_svd[,1])),
            abs(cor(F_matrix[,2], PC_svd[,2])),
            abs(cor(F_matrix[,3], PC_svd[,3]))
        )
        sigma_mse[2, i] <- (1/(n * (100-3)) * sum(diag(D)[4:100]^2) - 1) ^ 2
    
    }
    corrs <- apply(corr_matrix, 1, mean)
    mses <- apply(sigma_mse, 1, mean)
    j <- ifelse(n == 500, 1, 2)
    results[,j] <- c(corrs, mses)
    
    png(paste('task1_boxplot_eigen_', n, '.png', sep=''))
    boxplot(t(corr_matrix[1:3, ]), xlab='Column of F matrix', ylab='Correlation',
            main='Eigenvalue decomposition')
    dev.off()
    
    png(paste('task1_boxplot_SVD_', n, '.png', sep=''))
    boxplot(t(corr_matrix[4:6, ]), xlab='Column of F matrix', ylab='Correlation',
            main='SVD')
    dev.off()
    
    png(paste('task1_boxplot_mse_', n, '.png', sep=''))
    boxplot(t(sigma_mse), xlab='Method', ylab='MSE of sigma',
            main='MSE of sigma for both methods')
    dev.off()
}

results.df <- as.data.frame(results, row.names=c(
    'corr_eigen_1',
    'corr_eigen_2',
    'corr_eigen_3',
    'corr_SVD_1',
    'corr_SVD_2',
    'corr_SVD_3',
    'mse_eigen',
    'mse_SVD'
))
colnames(results.df) <- c('n=500', 'n=2000')

write.csv(results.df, file='task1_table.csv')







