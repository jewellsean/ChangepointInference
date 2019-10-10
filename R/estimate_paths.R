# Created by: jewellsean
# Created on: 2018-09-26

## implement crops algo
## https://www.tandfonline.com/doi/abs/10.1080/10618600.2015.1116445

# CROPS algorithm
# input: A dataset y1:n = (y1,y2,...,yn);
# Minimum and maximum values of the penalty,
# lambda_min and lambda_max;
# output: The details of optimal segmentations

## get fpop stats
## input: fpop object
## output: dataframe with cols (lambda, m(lambda), Q_m(lambda)(y1:n))
fpop_stats <- function(fit) {
  df <- data.frame(lambda = fit$tuning_parameter,
                   changepoints_n = length(fit$change_pts),
                   cost = fit$cost[length(fit$cost)] - length(fit$change_pts) * fit$tuning_parameter)
  return(df)
}

update_path_stats <- function(df, new_fit) {
  return(rbind(df, fpop_stats(new_fit)))
}

get_num_changepts <- function(lambda, df) {
  df$changepoints_n[df$lambda == lambda][1]
}

get_cost <- function(lambda, df) {
  df$cost[df$lambda == lambda][1]
}


#' Estimate changepoints using
#' automatic tuning parameter selection within
#' a range of values [lambda_min, lambda_max]
#'
#' @param dat observed data
#' @param lambda_min minimum lamba value
#' @param lambda_max maximum lamba value
#' @param max_iters maximum number of tuning parameters attempted
#'
#' @return Returns a list with elements:
#' @return \code{path_stats} a dataframe with summary statistics (number of changepoints, tuning parameters, cost)
#' @return \code{path_fits} a list containing the fitted model for each tuning parameter
#' @return \code{approximate_path} a boolean indicating whether an early stopping criterion condition occurred
#'
#' @details
#' Uses the CROPS algorithm (https://www.tandfonline.com/doi/abs/10.1080/10618600.2015.1116445) to efficiently 
#' choose values of the tuning parameter between [lambda_min, lambda_max].
#' 
#' @examples
#' 
#' set.seed(1)
#' mu <- rep(c(1, 3), each = 100)
#' dat <- mu + rnorm(length(mu))
#' 
#' fits <- estimate_paths(dat, lambda_min = 2, lambda_max = 3)
#'  
#' @export
estimate_paths <- function(dat, lambda_min = 1e-2, lambda_max = 1e1, max_iters = 10) {
  lambdas_used <- c(lambda_min, lambda_max)
  path_fits <- list()
  path_stats <- NULL
  approximate_path <- FALSE
  
  ## input validity
  stopifnot(lambda_min > 0)
  stopifnot(lambda_max > 0)
  stopifnot(max_iters > 4)
  stopifnot(is.numeric(dat))
  
  ## 1. Run CPD for penalty values lambda_min and lambda_max
  path_fits[[1]] <- changepoint_estimates(dat, "L0", lambda_min)
  path_stats <- update_path_stats(path_stats, path_fits[[1]])
  path_fits[[2]] <- changepoint_estimates(dat, "L0", lambda_max)
  path_stats <- update_path_stats(path_stats, path_fits[[2]])
  
  n_fits <- 2
  insert_ind <- 2
  ## 2. Set lambda_star = {[lambda_min, lambda_max]}
  lambda_star <- list(lambdas_used)
  
  while (length(lambda_star) > 0 &&
         n_fits <= max_iters) {
    # 3. Choose an element of lambda_star; denote this element as [lambda_0,lambda_1];
    # here always take the first element of list
    ind <- sample(1:length(lambda_star), size = 1)
    current_interal <- lambda_star[[ind]]
    
    if (get_num_changepts(current_interal[1], path_stats) >
        get_num_changepts(current_interal[2], path_stats) + 1) {
      lambda_int <- (get_cost(current_interal[2], path_stats) -
                       get_cost(current_interal[1], path_stats)) / (
                         get_num_changepts(current_interal[1], path_stats) -
                           get_num_changepts(current_interal[2], path_stats)
                       )
      
      n_fits <- n_fits + 1
      test_fit <- changepoint_estimates(dat, "L0", lambda_int)  
      
      if (!(length(test_fit$change_pts) %in% path_stats$changepoints_n)) {
        insert_ind <- insert_ind + 1
        path_fits[[insert_ind]] <- test_fit
        path_stats <- update_path_stats(path_stats, path_fits[[insert_ind]])
        
        if (get_num_changepts(lambda_int, path_stats) !=
            get_num_changepts(current_interal[2], path_stats)) {
          ## Set lambda_star = {lambda_star,[lambda_0,lambda_int),[lambda_int,lambda_1]}.;
          n_intervals <- length(lambda_star)
          lambda_star[[n_intervals + 1]] <- c(current_interal[1], lambda_int)
          lambda_star[[n_intervals + 2]] <- c(lambda_int, current_interal[2])
        }
      }
      lambda_star[[ind]] <- NULL
      }
      
      
    
    if (n_fits == max_iters) {
      warning(paste0("Full search path terminated early since maximum number of iterations (", max_iters,
                     ") reached. Rerun with larger 'max_iter' parameter for full path."))
      approximate_path = TRUE
    }
  }
  print(n_fits)
  out <- list(path_stats = path_stats, path_fits = path_fits, approximate_path = approximate_path)
  class(out) <- "estimated_paths"
  
  return(out)
}
