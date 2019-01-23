require(xgboost)
require(purrr) # partial function

library(geozoo) # create multidimensional objects (grid)
library(scales) # rescale numerical values


# This is a method to test the following method for estimating conditional quantiles using XGBoost.
# A data set is simulated and passed on to the xgboost.quantiles() function.
xgboost.quantiles.test <- function() {
       
       # simulation of test data
       n         <- 1000
       dimension <- 1
       interval  <- c(-5, 5)
       data.X    <- matrix(runif(n * dimension, interval[1], interval[2]), nrow = dimension)
       test.x    <- create.grid(floor(sqrt(n)), dimension, interval[1], interval[2])
       data.Y    <- apply(data.X ^ 2, 2, sum) + rnorm(n)
       
       # simulation of f(x) = x * sin(x)
       # data.Y <- data.X * sin(data.X)
       # dy     <- 1.5 + 1.0 * runif(length(data.Y), 0, 1)
       # noise  <- rnorm(500, 0, dy)
       # data.Y <- data.Y + noise
      
       if (dimension == 1) {
              draw <- TRUE
       } else {
              draw <- FALSE
       }
       
       xgboost.quantiles(
              data.X = t(data.X),
              data.Y = data.Y,
              test.x = test.x,
              draw   = draw
       )
}

# This is a method to estimate conditional quantiles using XGBoost. 
# It requires training data and test points for which the alpha-quantiles will be returned.
xgboost.quantiles <- function(data.X, data.Y, test.x, alpha = c(0.05, 0.25, 0.5, 0.75, 0.95), draw = FALSE) {
       
       dtrain <- xgb.DMatrix(data = as.matrix(data.X), label = as.matrix(data.Y))
       
       params.obj <- list(
              # controls the smoothing of the step function
              # the higher the delta, the steeper the connection is between the two steps.
              delta = 1,
              
              # controls the size of the interval where the cost-function is similar to
              # the usual quantile regression cost function
              threshold = 5.5,
              
              # controls how drastically you want to remix predictions within the leafs. 
              # The higher rho, the more the gradient is randomized.
              rho = 3.7
       )
       
       predictions <- array(0, dim = c(length(alpha), nrow(test.x)))
       
       for (i in 1:length(alpha)) {
              
              # training of the model
              model <- xgb.train(
                     
                     # Maximum depth of a tree. 
                     # Increasing this value will make the model more complex and more likely to overfit.
                     # range: [0,inf]
                     max_depth = 8,
                     
                     # L1 regularization term on weights. 
                     # Increasing this value will make model more conservative.
                     # default = 1
                     alpha = 1,
                     
                     # Minimum loss reduction required to make a further partition on a leaf node. 
                     # The larger gamma is, the more conservative the algorithm will be.
                     # range: [0,inf]
                     gamma = 0.5,
                     
                     # training data
                     data = dtrain,
                     
                     # amount of iterations
                     nrounds = 100, 
                     
                     # objective function
                     obj = partial(
                            
                            # Choose between the smoothed or unsmoothed variant of the quantile loss function.
                            # Note that the smoothed variant requires smoothing parameters.
                            loss.quantile.unsmoothed, 
                            
                            # Note that it may be useful to change the parameters corresponding to alpha.
                            alpha     = alpha[i]
                            # delta     = params.obj$delta,
                            # threshold = params.obj$threshold,
                            # rho       = params.obj$rho
                     ),
                     
                     # evaluation function
                     feval = partial(loss.quantile.score, alpha = alpha)
              )
              
              # predict using the trained model
              predictions[i, ] <- predict(model, as.matrix(test.x))
       }
       
       # Note that drawing only works with onedimensional x.
       if (draw) {
              colours <- c(
                     "orange", 
                     "dark green", 
                     "dark blue", 
                     "dark red", 
                     "violet"
              )
              
              plot(data.X, data.Y, col = "light gray", bty = "n")
              
              for (i in 1:length(alpha)) {
                     df <- rbind(t(test.x), predictions[i, ])
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
       
       return(predictions)
}


# This method is a customized evaluation function for the XGBoost package.
# It requires the true and predicted values of a response variable 
# and returns the metric name and the result.
loss.quantile.score <- function(preds, dtrain, alpha) {
       labels <- getinfo(dtrain, "label")
       score  <- array(0, dim = length(preds))
       
       for (i in 1:length(preds)) {
              if (labels[i] < preds[i]) {
                     score[i] <- (alpha - 1) * (labels[i] - preds[i])
              } else {
                     score[i] <- alpha * (labels[i] - preds[i])
              }
       }
      
       score <- sum(score)
       
       return(list(metric = "error", value = score))
}


# This method requires the true and predicted values of a response variable,
# the wanted quantiles alpha and certain smoothing parameters found by grid search 
# to return smoothed variants of both gradient and hessian of the quantile loss function.
loss.quantile.smoothed <- function(
              preds, 
              dtrain,
              alpha     = 0.5,
              delta     = 1.0, 
              threshold = 5.0, 
              rho       = 3.5
       ) {
       
       labels   <- getinfo(dtrain, "label")
       x        <- labels - preds
       gradient <- array(0, dim = length(x))
       hessian  <- array(0, dim = length(x))
       
       for (i in 1:length(x)) {
              # force a split by adding randomization
              rho <- (2 * sample(0:1, 1) - 1.0) * rho
              
              # smooth the gradient and the hessian of the quantile loss function
              if (abs(x[i]) > threshold) {
                     gradient[i] <- rho
                     hessian [i] <- 1
              } else if (x[i] < (alpha - 1) * delta) {
                     gradient[i] <- 1 - alpha
                     hessian [i] <- 0
              } else if (x[i] < alpha * delta) {
                     gradient[i] <- -x[i] / delta
                     hessian [i] <- 1.0 / delta
              } else {
                     gradient[i] <- -alpha
                     hessian [i] <- 0
              }
       }
       
       # plot smoothed variants of both gradient and hessian
       # plot(x, gradient, type = "p")
       # plot(x, hessian , type = "p")
       
       return(list(grad = gradient, hess = hessian))
}


# This method requires the true and predicted values of a response variable,
# the wanted quantiles alpha and certain smoothing parameters found by grid search 
# to return both gradient and hessian of the quantile loss function.
loss.quantile.unsmoothed <- function(preds, dtrain, alpha = 0.5) {
       
       labels   <- getinfo(dtrain, "label")
       x        <- labels - preds
       n        <- 5
       gradient <- array(0, dim = length(x))
       hessian  <- array(0, dim = length(x))
       
       for (i in 1:length(x)) {
              if (x[i] < (-1.0) / n) {
                     gradient[i] <- 1 - alpha 
                     hessian [i] <- 0
              } else if (x[i] <= 0) {
                     gradient[i] <- (-1) * n / 2 * x[i] 
                     hessian [i] <- n ^ 2 * x[i] + n
              } else if (x[i] <= 1.0 / n) {
                     gradient[i] <- (-1) * n / 2 * x[i] 
                     hessian [i] <- (-1) * n ^ 2 * x[i] + n
              } else {
                     gradient[i] <- (-1) * alpha
                     hessian [i] <- 0
              }
       } 
       
       return(list(grad = gradient, hess = hessian))
}


# This method creates a p-dimensional grid with n equispaced points each dimension between 0 and 1.
# The method then rescales each value to a specified range.
create.grid <- function(points.amount, dimension, range.min, range.max) {
       grid <- cube.solid.grid(p = dimension, n = points.amount)$points
       grid <- rescale(grid, to = c(range.min, range.max))
       
       return(grid)
}

