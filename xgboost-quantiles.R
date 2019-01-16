require(xgboost)
require(purrr) # partial function


# This is a method to test the following method for estimating conditional quantiles using XGBoost.
# A data set is simulated and passed on to the xgboost.quantiles() function.
xgboost.quantiles.test <- function() {
       x     <- sort(runif(10000, 0, 10))
       y     <- x * sin(x)
       dy    <- 1.5 + 1.0 * runif(length(y), 0, 1)
       noise <- rnorm(10000, 0, dy)
       y     <- y + noise
       
       xgboost.quantiles(x, y)
}


# This is a method for estimating conditional quantiles using XGBoost given data points (x, y).
# Note that x is required to be onedimensional.
xgboost.quantiles <- function(x, y) {
       
       # TODO: Train and test data have to be splitted/different.
       dtrain <- xgb.DMatrix(as.matrix(X), label = as.matrix(y))
       dtest  <- xgb.DMatrix(as.matrix(X), label = as.matrix(y))
       
       # observation of train and test errors
       watchlist <- list(train = dtrain, eval = dtest)
       
       params.xgb <- list(
              # Maximum depth of a tree. 
              # Increasing this value will make the model more complex and more likely to overfit.
              # range: [0,inf]
              max_depth = 3,
              
              # L1 regularization term on weights. 
              # Increasing this value will make model more conservative.
              # default = 1
              alpha = 5.0,
              
              # Minimum loss reduction required to make a further partition on a leaf node. 
              # The larger gamma is, the more conservative the algorithm will be.
              # range: [0,inf]
              gamma = 0.5
       )
       
       params.obj <- list(
              # determine the alpha-quantiles
              alpha = c(0.05, 0.5, 0.95),
              
              # controls the smoothing of the step function
              # the higher the delta, the steeper the connection is between the two steps.
              delta = 1,
              
              # controls the size of the interval where the cost-function is similar to
              # the usual quantile regression cost function
              threshold = c(6, 5.5, 5),
              
              # controls how drastically you want to remix predictions within the leafs. 
              # The higher rho, the more the gradient is randomized.
              rho  = c(3.2, 3.7, 4.2),
              
              # colours for plotting the predictions
              colour = c("orange", "dark green", "dark blue")
       )
       
       plot (X, y         , col = "light gray")
       lines(X, X * sin(X), col = "dark red"  )
       
       for (i in 1:length(params.obj$alpha)) {
              
              # training of the model
              model.boosted <- xgb.train(
                     params    = params.xgb,
                     data      = dtrain,
                     nrounds   = 50, 
                     watchlist = watchlist,
                     obj       = partial(
                            loss.quantile, 
                            alpha         = params.obj$alpha    [i], 
                            threshold     = params.obj$threshold[i], 
                            rho           = params.obj$rho      [i]
                     ),
                     feval = partial(loss.quantile.score, alpha = params.obj$alpha[i])
              )
              
              # perform the prediction
              predictions <- predict(model.boosted, dtest)
              
              lines(X, predictions, col = params.obj$colour[i])
       }
}


# This method is a customized evaluation function for the XGBoost package.
# It requires the true and predicted values of a response variable 
# and returns the metric name and the result.
loss.quantile.score <- function(preds, dtrain, alpha) {
       labels <- getinfo(dtrain, "label")
       score  <- array(0, dim = length(preds))
       
       for (i in 1:length(preds)) {
              if (labels[i] < preds[i]) {
                     score[i] <- (alpha - 1) * abs(labels[i] - preds[i])
              } else {
                     score[i] <- alpha * abs(labels[i] - preds[i])
              }
       }
      
       score <- 1.0 / sum(score)
       
       return(list(metric = "error", value = score))
}


# This method requires the true and predicted values of a response variable,
# the wanted quantiles alpha and certain smoothing parameters found by grid search 
# to return smoothed variants of both gradient and hessian of the quantile loss function.
loss.quantile <- function(
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
              gradient[i] <- loss.quantiles.gradient(x[i], alpha, delta, threshold, rho)
              hessian [i] <- loss.quantiles.hessian (x[i], alpha, delta, threshold)
       }
       
       # plot smoothed variants of both gradient and hessian
       # plot(x, gradient, type = "p")
       # plot(x, hessian , type = "p")
       
       return(list(grad = gradient, hess = hessian))
}


# This method provides a smoothed variant of the gradient of the quantile loss function.
loss.quantile.gradient <- function(x, alpha, delta, threshold, rho) {
       
       # smooth the gradient of the quantile loss function
       if ((x >= (alpha - 1) * delta) & (x < alpha * delta)) {
              gradient <- (-1) * x / delta
       } else if (x < (alpha - 1) * delta) {
              gradient <- 1 - alpha
       } else if (x >= alpha * delta) {
              gradient <- (-1) * alpha
       }
       
       # force a split by adding randomization
       rho <- (2 * sample(0:1, 1) - 1.0) * rho
       
       if (abs(x) >= threshold) {
              gradient <- rho
       }
       return(gradient)
}


# This method provides a smoothed variant of the hessian of the quantile loss function.
loss.quantile.hessian <- function(x, alpha, delta, threshold) {
       
       # smooth the hessian of the quantile loss function
       hessian <- ((x >= (alpha - 1) * delta) & (x < alpha * delta)) / delta
       
       if (abs(x) >= threshold) {
              hessian <- 1
       }
       return(hessian)
}

