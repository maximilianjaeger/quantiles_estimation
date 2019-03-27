library(geozoo) # create multidimensional objects (grid)
library(scales) # rescale numerical values
library(splines)
library(quantreg)


# This is a method to test the following method for estimating conditional quantiles using Smoothed Splines.
# A data set is simulated and passed on to the qss.quantiles() function.
qss.quantiles.test <- function() {
       n         <- 500
       dimension <- 1
       interval  <- c(-2, 2)
       data.X    <- matrix(runif(n * dimension, interval[1], interval[2]), nrow = dimension)
       data.Y    <- apply(data.X ^ 2, 2, sum) + rnorm(n)
       test.x    <- create.grid2(10, dimension, data.X)
       test.x    <- matrix(sort(test.x))
       alpha     <- c(0.05, 0.25, 0.5, 0.75, 0.95)
       
       quantiles <- qss.quantiles(data.X, data.Y, test.x, alpha, draw = TRUE)
}


# This is a method to estimate conditional quantiles using Smoothed Splines. 
# It requires training data and test points for which the alpha-quantiles will be returned.
qss.quantiles <- function(data.X, data.Y, test.x, alpha = c(0.05, 0.25, 0.5, 0.75, 0.95), draw = FALSE) {

       quantiles <- array(NA, dim = c(length(alpha), nrow(test.x)))

       for (i in 1:length(alpha)) {
              x          <- data.X[1, ]
              lambda     <- seq(0, 20, by = 0.2)
              models     <- vector(mode = "list", length = length(lambda))
              models.aic <- array(NA, dim = length(lambda))
              
              for (j in 1:length(lambda)) {
                     models    [[j]] <- rqss(data.Y ~ qss(x, lambda = lambda[j]), tau = alpha[i])
                     models.aic [j]  <- AIC(models[[j]])
              }
              
              model.qss      <- models[[which.min(abs(models.aic))]]
              quantiles[i, ] <- predict(model.qss, newdata = data.frame(x = test.x))
       }
       
       if (draw) {
              quantiles.plot(data.X, data.Y, test.x, "Quantile Smoothing Splines", alpha, quantiles)
       }
       
       return(quantiles)
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
# The method then rescales each value to a specified range within the convex hull of the fitting data.
create.grid2 <- function(points.amount, dimension, data.X) {
       range.min <- min(data.X)
       range.max <- max(data.X)
       
       # if (dimension > 1) {
       #        for (i in 1:dimension) {
       #               if (range.min < min(data.X[i, ])) {
       #                      range.min <- min(data.X[i, ])
       #               }
       #               
       #               if (range.max > max(data.X[i, ])) {
       #                      range.max <- max(data.X[i, ])
       #               }
       #        }
       # }
       
       grid <- cube.solid.grid(p = dimension, n = points.amount)$points
       grid <- rescale(grid, to = c(range.min, range.max))
       
       return(grid)
}


