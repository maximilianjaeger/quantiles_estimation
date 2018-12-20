library("RANN")
library("caTools")
library("caret")


# This is a method to estimate the conditional quantiles of a given point
# using only the values of its k nearest neighbors in a given data set.
kNN.quantiles <- function(data.X, data.Y, x, alpha = c(0.05, 0.25, 0.5, 0.75, 0.95), amount.neighbors) {

       # find the nearest neighbors of x in data.X
       neighbors <- nn2(t(data.X), query = x, k = amount.neighbors)[[1]]
       
       # calculate the quantiles for x using the nearest neighbors
       x.quantiles <- quantile(data.Y[neighbors], probs = alpha, na.rm = TRUE)
       
       return(x.quantiles)
}


# This is a method for finding the optimal amount of nearest neighbors for estimating conditional quantiles.
# Therefore the given data set is split into a training and a test set for cross validation of the estimations.
kNN.quantiles.crossValidation <- function(data.X, data.Y, alpha = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
       # set.seed(00)
       
       # Data is split into different training and test sets.
       folds.amount <- 2
       folds        <- createFolds(data.Y, k = folds.amount)
       
       # Set containing different amount of neigbors that is to be optimized.
       neighbors <- c(1, 2,3,4, 5, 10, 20, 50)
       
       # array containing the results of the check function (cost function of quantile regression)
       quantiles.error <- array(NA, dim = c(length(alpha), length(neighbors), length(data.Y), folds.amount))
       
       for (f in 1:folds.amount) {
              data.X.train <- data.X[, -folds[[f]]]
              data.X.test  <- data.X[,  folds[[f]]]
              data.Y.train <- data.Y[-folds[[f]]]
              data.Y.test  <- data.Y[ folds[[f]]]
              
              # array containing the estimated quantiles for different alpha, neighbors and test points
              quantiles <- array(NA, dim = c(length(alpha), length(neighbors), length(data.Y.test)))
              
              for (i in 1:length(neighbors)) {
                     for (j in 1:length(data.Y.test)) {
                            quantiles[, i, j] <- kNN.quantiles(data.X.train, data.Y.train, t(data.X.test[, j]), alpha, neighbors[i])[[2]]
                            
                            for (a in 1:length(alpha)) {
                                   quantiles.error[a, i, j, f] <- check.function(data.Y.test[j], quantiles[a, i, j], alpha[a])
                            }
                     }
              }
       }
       
       # mean of errors over all test points
       k.optimum <- apply(quantiles.error, c(1, 2, 4), mean, na.rm = TRUE)
       
       # mean of errors over all folds
       k.optimum <- apply(k.optimum, c(1, 2), mean, na.rm = TRUE)
       
       # mean of errors over all alpha
       k.optimum <- apply(k.optimum, c(2), mean, na.rm = TRUE)
       
       # return index of the minimum error / optimal amount of neighbors
       return(neighbors[which.min(k.optimum)])
}


# Cost function of quantile regression (check function)
check.function <- function(value.Y, quantile, alpha) {
       if (value.Y < quantile) {
              return((alpha - 1) * (value.Y - quantile))
       } else {
              return(alpha * abs(value.Y - quantile))
       }
}


# This is a method to declare data sets for testing of the functions above.
kNN.quantiles.test <- function() {
       # set.seed(1)
       
       n      <- 500
       d      <- 2
       data.X <- matrix(runif(n * d, -2, 2), nrow = d)
       data.Y <- apply(data.X, 2, sum) + 2*rnorm(n)
       x      <- array(0, dim = c(1, d))
       
       amount.neighbors.optimum <- kNN.quantiles.crossValidation(data.X = data.X, data.Y = data.Y)
       results                  <- kNN.quantiles(data.X, data.Y, x, amount.neighbors = amount.neighbors.optimum)
       
       return(list(k = amount.neighbors.optimum, q = results))
}

