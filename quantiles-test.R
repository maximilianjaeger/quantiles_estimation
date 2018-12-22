library(geozoo) # create multidimensional objects (grid)
library(scales) # rescale numerical values
library(ggplot2)


# This method is used for comparing the estimation of the conditional quantiles via kNN 
# and the quantization algorithm.
quantiles.test <- function(variables) {
       
       # test set for both estimation models
       n      <- 500
       d      <- 2
       data.X <- matrix(runif(n * d, -2, 2), nrow = d)
       data.Y <- apply(data.X ^ 2, 2, sum) + rnorm(n)
       x      <- create.grid(10, d, -2, 2)
       alpha  <- c(0.05, 0.25, 0.5, 0.75, 0.95)
       testN  <- c(10, 15, 30, 60, 100, 150)
       B      <- 20
       tildeB <- 15
       
       # calculation of the quantiles for every point in x using kNN
       amount.neighbors.optimum <- kNN.quantiles.crossValidation(data.X, data.Y, alpha)
       quantiles.knn            <- array(NA, dim = c(length(alpha), nrow(x)))
       
       for (i in 1:nrow(x)) {
              quantiles.knn[, i] <- kNN.quantiles(data.X, data.Y, t(x[i, ]), alpha, amount.neighbors.optimum)
       }
       
       # calculation of the quantiles for every point in x using quantization
       quantiles.quantif <- quantif.quantiles.max(
              X      = data.X,
              Y      = data.Y,
              x      = t(x),
              alpha  = alpha,
              B      = B,
              tildeB = tildeB,
              testN  = testN,
              ncores = 3
       )$hatq_opt
       
       
       # evalution of the models
       evaluation.mean    <- apply(quantiles.knn - quantiles.quantif, 1, mean)
       evaluation.sd      <- apply(quantiles.knn - quantiles.quantif, 1, sd)
       evaluation.summary <- apply((quantiles.knn - quantiles.quantif) ^ 2, 1, summary)
       evaluation.table   <- rbind(evaluation.mean, evaluation.sd, evaluation.summary)
       evaluation.table   <- cbind(evaluation.table, apply(evaluation.table, 1, mean))
       
       colnames(evaluation.table) <- c(alpha, "Mean")
       rownames(evaluation.table) <- c(
              "Mean Error", 
              "SD", 
              "Min. Qu.",
              "1st Qu.", 
              "Median Qu.", 
              "Mean Qu.", 
              "3rd Qu.", 
              "Max. Qu."
       )
       
       return(round(evaluation.table, 2))
}


# This method creates a p-dimensional grid with n equispaced points each dimension between 0 and 1.
# The method then rescales each value to a specified range.
create.grid <- function(points.amount, dimension, range.min, range.max) {
       grid <- cube.solid.grid(p = dimension, n = points.amount)$points
       grid <- rescale(grid, to = c(range.min, range.max))
       
       return(grid)
}

