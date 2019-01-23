library(RANN)
library(caTools)
library(caret)
library(geozoo) # create multidimensional objects (grid)
library(scales) # rescale numerical values


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

       # Data is split into different training and test sets.
       folds.amount <- 2
       folds        <- createFolds(data.Y, k = folds.amount)
       
       # Set containing different amount of neigbors that is to be optimized.
       neighbors <- c(1, 2, 3, 4, 5, 10, 20)
       
       # array containing the results of the check function (cost function of quantile regression)
       quantiles.error <- array(NA, dim = c(length(alpha), length(neighbors), length(data.Y), folds.amount))
       
       for (f in 1:folds.amount) {
              
              # check for dimensionality issues with the data and split data into train and test sets
              if (length(data.X) == length(data.Y)) {
                     data.X.train <- data.X[-folds[[f]]]
                     data.X.test  <- data.X[ folds[[f]]]
              } else {
                     data.X.train <- data.X[, -folds[[f]]]
                     data.X.test  <- data.X[,  folds[[f]]]
              }
              data.Y.train <- data.Y[-folds[[f]]]
              data.Y.test  <- data.Y[ folds[[f]]]
              
              # array containing the estimated quantiles for different alpha, neighbors and test points
              quantiles <- array(NA, dim = c(length(alpha), length(neighbors), length(data.Y.test)))
              
              for (i in 1:length(neighbors)) {
                     for (j in 1:length(data.Y.test)) {
                            
                            # check for dimensionality issues with the data
                            if (length(data.X) == length(data.Y)) {
                                   tmp1 <- t(data.X.train)
                                   tmp2 <- t(data.X.test[j])
                            } else {
                                   tmp1 <- data.X.train
                                   tmp2 <- t(data.X.test[, j])
                            }
                            
                            # estimate the quantiles
                            quantiles[, i, j] <- kNN.quantiles(
                                   data.X           = tmp1, 
                                   data.Y           = data.Y.train, 
                                   x                = tmp2, 
                                   alpha            = alpha, 
                                   amount.neighbors = neighbors[i]
                            )[[2]]
                            
                            # calculate the errors
                            for (a in 1:length(alpha)) {
                                   quantiles.error[a, i, j, f] <- check.function(
                                          value.Y  = data.Y.test[j], 
                                          quantile = quantiles[a, i, j], 
                                          alpha    = alpha[a]
                                   )
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
              return(alpha * (value.Y - quantile))
       }
}


# This is a method to declare data sets for testing of the functions above.
kNN.quantiles.test <- function() {
       
       # simulation of data
       n         <- 500
       dimension <- 1
       interval  <- c(-5, 5)
       data.X    <- matrix(runif(n * dimension, interval[1], interval[2]), nrow = dimension)
       data.Y    <- apply(data.X ^ 2, 2, sum) + rnorm(n)
       test.x    <- create.grid(n / 10, dimension, interval[1], interval[2])
       alpha     <- c(0.05, 0.25, 0.5, 0.75, 0.95)
       
       # calculation of the quantiles for every point in x using kNN
       amount.neighbors.optimum <- kNN.quantiles.crossValidation(data.X, data.Y, alpha)
       quantiles.knn            <- array(NA, dim = c(length(alpha), nrow(test.x)))
       
       for (i in 1:nrow(test.x)) {
              quantiles.knn[, i] <- kNN.quantiles(
                     data.X           = data.X, 
                     data.Y           = data.Y, 
                     x                = t(test.x[i, ]), 
                     alpha            = alpha, 
                     amount.neighbors = amount.neighbors.optimum
              )
       }
       
       if (dimension == 1) {
              colours <- c(
                     "orange", 
                     "dark green", 
                     "dark blue", 
                     "dark red", 
                     "violet"
              )
              
              plot(data.X, data.Y, col = "light gray", bty = "n")
              
              for (i in 1:length(alpha)) {
                     df <- rbind(t(test.x), quantiles.knn[i, ])
                     df <- df[, order(df[1, ])]
                     
                     lines(df[1, ], df[2, ], col = colours[i])   
              }
              
              legend(
                     "topleft",
                     inset     = 0.02,
                     cex       = 0.75,
                     ncol      = length(alpha),
                     legend    = alpha,
                     lty       = array(1, dim = length(alpha)),
                     col       = colours,
                     box.lty   = 0
              )
       }
}


# This method creates a p-dimensional grid with n equispaced points each dimension between 0 and 1.
# The method then rescales each value to a specified range.
create.grid <- function(points.amount, dimension, range.min, range.max) {
       grid <- cube.solid.grid(p = dimension, n = points.amount)$points
       grid <- rescale(grid, to = c(range.min, range.max))
       
       return(grid)
}

