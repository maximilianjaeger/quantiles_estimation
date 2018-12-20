library(Biobase)
library(RANN) # nearest neighbour search algorithm - nn2
library(Hmisc) # calculation of weighted quantiles - wtd.quantile
library(QuantifQuantile)
library(parallel) # parallel processing - use of more than one kernel


QuantifQuantile.max <- function (
    X, 
    Y, 
    x, 
    alpha  = c(0.05, 0.25, 0.5, 0.75, 0.95), 
    testN  = seq(5, 70, by = 5), 
    p      = 2, 
    B      = 50, 
    tildeB = 20, 
    same_N = TRUE,
    type   = 7, 
    ncores = 1
    ) {
    
    # determines how much real and CPU time (in seconds) the currently running R process has already taken
    total_time <- proc.time() 
    
    # control for infeasible input data
    input.control(X, x, Y, testN, alpha, B, tildeB, p, same_N)
    
    if (is.vector(X)) {
        d <- 1
        n <- length(X)
        X <- matrix(X, nrow = d)
        x <- matrix(x, nrow = d)
    } else {
        n <- ncol(X)
        d <- nrow(X)
    }
    
    m        <- length(x) / d
    hatISE_N <- array(0, dim = c(length(alpha), length(testN))) 
    hatq_N   <- array(0, dim = c(length(alpha), m, length(testN))) 
    primeX   <- array(X[, sample(c(1:n), n * (B + tildeB), replace = T)], dim = c(d, n, B + tildeB))
    
    calc_hatq_N <- function(N) {
        # This function provides ng optimal quantization grids for X, with N fixed. 
        hatX <- choice.grid(X, N, ng = (B + tildeB))$opti_grid 
        
        # indeces and distances of the 3 grids on which X is projected
        amount.neighbours <- 1
        tmp               <- projection(hatX, X, amount.neighbours, n, B + tildeB)
        proj.X.index      <- tmp[[1]]
        proj.X.weight     <- tmp[[2]]
        
        # estimation of q_alpha(x) for N fixed 
        # save the B+tildeB estimation of q_alpha(x)
        Hatq <- array(0, dim = c(m, length(alpha), B + tildeB))
        
        # save by Voronoi cell
        Hatq_cell <- array(0, dim = c(N, length(alpha), B + tildeB))
        
        # index of the grid on which x is projected
        projection.x.index <- projection(hatX, x, 1, m, B + tildeB)[[1]]
        
        # calculation of the conditional quantile for each cell
        # Since any point of a cell is projected on the center of this cell, the corresponding conditional quantiles are equal
        for (k in 1:(B + tildeB)) {
            for (j in 1:N) {
                # find the data points that belong to the Voronoi cell and their weights
                cell.elements <- NULL
                Y.weight      <- NULL
                
                for (l in 1:amount.neighbours) {
                    tmp           <- which(proj.X.index[, k, l] == j)
                    cell.elements <- c(cell.elements, tmp)
                    Y.weight      <- c(Y.weight, proj.X.weight[tmp, k, l])
                }
                
                Y.weight      <- Y.weight     [!duplicated(cell.elements)]
                cell.elements <- cell.elements[!duplicated(cell.elements)]
                
                # calculate the quantiles for the Voronoi cell
                if (length(cell.elements) > 0) {
                    Hatq_cell[j, , k] <- wtd.quantile(Y[cell.elements], weights = Y.weight, probs = alpha)
                } else {
                    Hatq_cell[j, , k] <- NA 
                }
            }
            Hatq[, , k] = Hatq_cell[projection.x.index[, k, 1], , k]
        }
        
        # the final estimation is the mean of the B estimations
        hatq <- array(0, dim = c(m, length(alpha)))
        
        for (k in 1:B) {
            hatq <- hatq + Hatq[, , k]
        }
        
        hatq <- hatq / B
        
        # the last tilde B are used to estimate the ISE
        HATq   <- array(rep(hatq, tildeB), dim = c(m, length(alpha), tildeB))
        hatISE <- (HATq - Hatq[, , c((1 + B):(B + tildeB)), drop = FALSE]) ^ 2
        hatISE <- apply(hatISE, 2, sum) / (m * tildeB)
        
        print(N)
        list(hatq = hatq, hatISE = hatISE)
    }
    
    parallel_hatq_hatISE <- mclapply(testN, calc_hatq_N, mc.cores = ncores, mc.set.seed = F)
    
    for (i in 1:length(testN)) {
        hatq_N[, , i] <- t(parallel_hatq_hatISE[[i]]$hatq) 
        hatISE_N[, i] <- parallel_hatq_hatISE[[i]]$hatISE
    }
    
    if (same_N) {
        # choice of optimal N
        hatISEmean_N <- apply(hatISE_N, 2, mean)
        i_opt        <- which.min(hatISEmean_N)
        
        # optimal value for N chosen as minimizing the sum of hatISE for the different alpha's
        N_opt <- testN[i_opt]
        
        # table of the associated estimated conditional quantiles
        hatq_opt <- hatq_N[, , i_opt, drop = F]
        hatq_opt <- matrix(hatq_opt, nrow = length(alpha))
    } else {
        # choice of optimal N
        i_opt <- apply(hatISE_N, 1, which.min)
        
        # optimal value for N chosen as minimizing the sum of hatISE for the different alpha's
        N_opt <- testN[i_opt]
        
        # table of the associated estimated conditional quantiles
        hatq_opt <- array(0, dim = c(length(alpha), length(x)/d))
        
        for (i in 1:length(alpha)) {
            hatq_opt[i, ] <- hatq_N[i, , i_opt[i]]
        }
    }
    
    if (length(N_opt) == 1) {
        if (N_opt == min(testN)) {
            warning("N_opt is on the left boundary of testN")
        }
        if (N_opt == max(testN)) {
            warning("N_opt is on the right boundary of testN")
        }
    } else {
        if (any(N_opt == min(testN))) {
            warning("N_opt is on the left boundary of testN for at least one value of alpha")
        }
        if (any(N_opt == max(testN))) {
            warning("N_opt is on the right boundary of testN for at least one value of alpha")
        }
    }
    
    rownames(hatISE_N) <- alpha
    colnames(hatISE_N) <- testN
    
    output <- list(
        hatq_opt      = hatq_opt, 
        fitted.values = hatq_opt, 
        N_opt         = N_opt, 
        hatISE_N      = hatISE_N, 
        hatq_N        = hatq_N, 
        X             = X, 
        Y             = Y, 
        x             = x, 
        alpha         = alpha, 
        testN         = testN, 
        B             = B, 
        tildeB        = tildeB, 
        type          = type, 
        total_time    = proc.time() - total_time
        )
    
    class(output) <- "QuantifQuantile.max"
    
    return(output)
}

#' This function controls for infeasible input data.
input.control <- function(X, x, Y, testN, alpha, B, tildeB, p, same_N) {
    if (!is.numeric(X)) {
        stop("Y must be numeric")
    }
    if (!is.numeric(x)) {
        stop("x must be numeric")
    }
    if (!is.vector(Y)) {
        stop("Y must be a vector")
    }
    if (!all(floor(testN) == testN & testN > 0)) {
        stop("testN must have entire positive entries")
    }
    if (!all(alpha > 0 & alpha < 1)) {
        stop("alpha must be strictly between 0 and 1")
    }
    if ((!(floor(B) == B)) | (B <= 0)) {
        stop("B must be a positive entire")
    }
    if ((!(floor(tildeB) == tildeB)) | (tildeB <= 0)) {
        stop("tildeB must be a positive entire")
    }
    if (p < 1) {
        stop("p must be at least 1")
    } 
    if (!is.logical(same_N)) {
        stop("same_N must be logical")
    }
    if (is.vector(X)) {
        d <- 1
        n <- length(X)
        X <- matrix(X, nrow = d)
        
        if (!is.vector(x)) {
            stop("x must have same dimension as X")
        }
        x <- matrix(x, nrow = d)
    } else {
        if (!is.matrix(X)) {
            stop("X must be a matrix with d rows")
        }
        n <- ncol(X)
        d <- nrow(X)
        
        if (!is.matrix(x)) {
            stop("x must be a matrix with d rows")
        }
        if (nrow(x) != d) {
            stop("x must be a matrix with d rows")
        }
    }
}


#' This function returns the indeces and weights of the nearest neighbours of @param X in a grid @param hatX.
#' The weights are calculated by the method of inverse distance weighting.
projection <- function(hatX, X, amount.neighbours, rows, cols) {
    proj.index    <- array(0, dim = c(rows, cols, amount.neighbours)) 
    proj.distance <- array(0, dim = c(rows, cols, amount.neighbours))
    proj.weight   <- array(1, dim = c(rows, cols, amount.neighbours))
    
    for (k in 1:cols) {
        # search for the nearest neighbours of X in hatX
        tmp <- nn2(t(hatX[, , k]), t(X), amount.neighbours)
        
        # safe the indeces and distances of the neighbours
        for (l in 1:amount.neighbours) {
            proj.index   [, k, l] <- tmp[[1]][, l]
            proj.distance[, k, l] <- tmp[[2]][, l]
        }
    }
    
    # inverse distance weighting
    pow <- 0.5
    
    for (j in 1:rows) {
        for (k in 1:cols) {
            full.inverse.dist <- 0
            
            if (proj.distance[j, k, l] <= 0.00001) {
                proj.weight[j, k, l] <- 1
            } else {
                for (l in 1:amount.neighbours) {
                    full.inverse.dist <- full.inverse.dist + 1 / (proj.distance[j, k, l] ^ pow)
                }
                
                for (l in 1:amount.neighbours) {
                    proj.weight[j, k, l] <- 1 / (full.inverse.dist * proj.distance[j, k, l] ^ pow)
                }
            }
        }
    }
    return(list(proj.index, proj.weight))
}

