library(geozoo) # create multidimensional objects (grid)
library(scales) # rescale numerical values


# This method is used for comparing the estimations of the conditional quantiles with different 
# algorithms to the theoretical quantiles of simulated data.
# Note that the results can be printed for onedimensional x.
quantiles.test <- function(variables) {
       
       # test set for both estimation models
       n         <- 1000
       dimension <- 1
       interval  <- c(-5, 5)
       data.X    <- matrix(runif(n * dimension, interval[1], interval[2]), nrow = dimension)
       data.Y    <- apply(data.X ^ 2, 2, sum) + rnorm(n)
       test.x    <- create.grid(50, dimension, interval[1], interval[2])
       alpha     <- c(0.05, 0.25, 0.5, 0.75, 0.95)
       testN     <- c(10, 15, 30, 60, 100, 150)
       B         <- 20
       tildeB    <- 15
       
       # calculation of the true quantiles
       quantiles <- array(NA, dim = c(length(alpha), nrow(test.x)))
       
       for (i in 1:length(alpha)) {
              if (dimension > 1) {
                     quantiles[i, ] <- apply(test.x ^ 2, 1, sum) + qnorm(alpha[i])
              } else {
                     quantiles[i, ] <- test.x ^ 2 + qnorm(alpha[i])
              }
       }
       
       # calculation of the quantiles for every point in x using kNN
       amount.neighbors.optimum <- kNN.quantiles.crossValidation(data.X, data.Y, alpha)
       quantiles.knn            <- array(NA, dim = c(length(alpha), nrow(test.x)))
       
       for (i in 1:nrow(test.x)) {
              quantiles.knn[, i] <- kNN.quantiles(data.X, data.Y, t(test.x[i, ]), alpha, amount.neighbors.optimum)
       }
       
       # calculation of the quantiles for every point in x using quantization
       # Note that (at the moment) this method only works with onedimensional x.
       if (dimension > 1) {
              quantiles.quantif <- quantif.quantiles.max(
                     X      = data.X,
                     Y      = data.Y,
                     x      = t(test.x),
                     alpha  = alpha,
                     B      = B,
                     tildeB = tildeB,
                     testN  = testN,
                     ncores = 3
              )$hatq_opt
       }
       
       # calculation of the quantiles for every point in x using XGBoost
       quantiles.xgboost <- xgboost.quantiles(
              data.X = t(data.X),
              data.Y = data.Y,
              test.x = test.x
       )
       
       # draw the different results
       if (dimension == 1) {
              quantiles.plot(data.X, data.Y, test.x, "Theoretical Quantiles", alpha, quantiles)
              quantiles.plot(data.X, data.Y, test.x, "KNN Quantiles", alpha, quantiles.knn)
              quantiles.plot(data.X, data.Y, test.x, "XGBoost Quantiles", alpha, quantiles.xgboost)
       }
       
       # evalution of the models
       if (dimension == 1) {
              models        <- list(quantiles.knn,   quantiles.xgboost  )
              names(models) <- list("KNN Quantiles", "XGBoost Quantiles")
       } else {
              models        <- list(quantiles.knn,   quantiles.xgboost  , quantiles.quantif  )
              names(models) <- list("KNN Quantiles", "XGBoost Quantiles", "Quantif Quantiles")
       }
       
       tables <- vector(mode = "list", length = length(models))
       
       for (i in 1:length(models)) {
              evaluation.mean    <- apply( quantiles - models[[i]]     , 1, mean   )
              evaluation.sd      <- apply( quantiles - models[[i]]     , 1, sd     )
              evaluation.summary <- apply((quantiles - models[[i]]) ^ 2, 1, summary)
              evaluation.table   <- rbind(evaluation.mean, evaluation.sd, evaluation.summary)
              evaluation.table   <- cbind(evaluation.table, apply(evaluation.table, 1, mean))
              
              colnames(evaluation.table) <- c(alpha, "Mean")
              rownames(evaluation.table) <- c(
                     "Mean Error",
                     "SD",
                     "Min.   Qu.",
                     "1st    Qu.",
                     "Median Qu.",
                     "Mean   Qu.",
                     "3rd    Qu.",
                     "Max.   Qu."
              )
              
              tables[[i]]        <- round(evaluation.table, 2)
              names(tables)[[i]] <- paste("Differences between theoretical Quantiles and ", names(models)[[i]])
       }

       return(tables)
}


# This is a method that plots the quantile estimations for onedimensional x.
quantiles.plot <- function(data.X, data.Y, test.x, title, alpha, quantiles) {
       
       plot(data.X, data.Y, col = "light gray", bty = "n")
       
       colours <- c(
              "orange", 
              "dark green", 
              "dark blue", 
              "dark red", 
              "violet"
       )
       
       legend(
              "top",
              title   = title,
              legend  = alpha,
              lty     = array(1, dim = length(alpha)),
              col     = colours,
              box.lty = 0
       )
       
       for (i in 1:length(alpha)) {
              dataframe <- rbind(t(test.x), quantiles[i, ])
              dataframe <- dataframe[, order(dataframe[1, ])]
              
              lines(dataframe[1, ], dataframe[2, ], col = colours[i])   
       }
}


# This method creates a p-dimensional grid with n equispaced points each dimension between 0 and 1.
# The method then rescales each value to a specified range.
create.grid <- function(points.amount, dimension, range.min, range.max) {
       grid <- cube.solid.grid(p = dimension, n = points.amount)$points
       grid <- rescale(grid, to = c(range.min, range.max))
       
       return(grid)
}

