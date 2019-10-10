#' Estimate changepoints via L0 and binary segmentation
#'
#' @param dat observed data
#' @param type specifies changepoint algorithm: L0 ("L0") or binary segmentation
#'   ("BS");
#' @param tuning_parameter L0 segmentation: tuning parameter \eqn{\lambda};
#'   Binary segmentation: number of segments K.
#' @param functional_pruning_out boolean specifying whether, in the case of L0
#'   segmentation, cost functions are returned. Defaults to FALSE.
#'
#'
#' @return For L0 segmentation, returns a list with elements:
#' @return \code{estimated_means} estimated means
#' @return \code{change_pts} the set of changepoints
#' @return \code{cost} the cost at each time point
#' @return \code{n_intervals} the number of piecewise quadratics used at each
#'   point
#' @return \code{piecewise_square_losses} dataframe of optimal cost functions
#'   Cost_s*(mu) for s = 1,..., T.
#' @return For binary segmentation, returns a list with elements
#' @return \code{change_pts} the set of changepoints
#' @return \code{orders} estimated changepoints in the order that they entered
#'   the model
#' @return \code{signs} the signs of the change in mean due to each changepoint
#'   (in the same order as orders)
#'
#' @details
#'
#' Consider the changepoint model \deqn{Y_t = \mu_t + \epsilon_t} where
#' \eqn{\epsilon_t iid N(0, \sigma^2)} for \eqn{ t=1,\ldots,T} and assume that
#' \eqn{\mu_1,\ldots,\mu_T} is  piecewise constant, in the sense that
#' \eqn{\mu_{\tau_j+1}=\mu_{\tau_j + 2 } = \ldots = \mu_{\tau_{j+1}}}, and
#' \eqn{\mu_{\tau_{j+1}} != \mu_{\tau_{j+1}+1}}, for \eqn{j=0,\ldots,K-1}, where
#' \eqn{0 = \tau_0 < \tau_1 < \ldots < \tau_K < \tau_{K+1} = T}, and where
#' \eqn{\tau_1,\ldots,\tau_K} represent the true changepoints.
#'
#' Estimation:
#'
#' This function estimates changepoints via L0 or binary segmentation based on
#' noisy observations \eqn{y_t, t = \ldots, T}. In the case of L0 segmentation,
#' changepoints are estimated by solving the optimization problem
#'
#' \deqn{minimize_{\mu_1,...,\mu_T} 0.5 \sum_{t=1}^T ( y_t - \mu_t )^2 + \lambda
#' \sum_{t=2}^T 1_[\mu_t != \mu_{t-1}]}
#'
#' for the global optimum. Our estimates for the changepoints correspond to the
#' breakpoints of \eqn{\hat{\mu}_1, \ldots, \hat{\mu}_T}. See Rigaill, G. (2015)
#' and Maidstone, R. et al. (2017) for more information and for discussion of
#' cost functions Cost_s*(mu).
#'
#' In the case of binary segmentation, changepoints are estimated by recursively
#' maximizing the CUSUM statistic: the first estimated changepoint maximizes the
#' CUSUM statistic over all possible locations. Subsequent changepoints are
#' estimated at the location that maximizes the CUSUM statistic of the data when
#' regions between previously estimated changepoints are considered. This
#' results in estimated changepoints \eqn{M(y) = {\hat{\tau}_1, ...,
#' \hat{\tau}_K}}, the order each changepoint entered the model \eqn{O(y)}, and
#' the sign of the change in mean due to each changepoint \eqn{\Delta(y)}. See
#' Vostrikova, L. (1981) for additional details.
#'
#' @examples
#'
#' ### Generate sample data
#' set.seed(1)
#' mu <- rep(c(1, 3, -4), each = 100)
#' dat <- mu + rnorm(length(mu))
#'
#'
#' ## L0 segmentation
#' fit <- changepoint_estimates(dat, "L0", 4, functional_pruning_out = TRUE)
#' print(fit)
#'
#' # Estimated changepoints
#' fit$change_pts
#'
#' # Plot optimal cost function at s = 300
#' # Colors denote the most recent changepoint associated with optimal cost at each mu
#' plot(fit, s = 300)
#' eval_cost_at.ChangepointInference_L0_estimated_changes(fit, s = 300, mu = -3)
#'
#' ### Binary segmentation
#' fit <- changepoint_estimates(dat, "BS", 2)
#' print(fit)
#'
#' # sorted changepoints, order changepoints were estimated,
#' # sign of change in mean due to changepoint (in same order as fit$ordered_change_pts)
#' fit$change_pts
#' fit$ordered_change_pts
#' fit$change_pt_signs
#' @seealso \code{\link{changepoint_estimates}},
#'   \code{\link{changepoint_inference}},
#'
#' @references
#'
#' Jewell, S., Fearnhead, P., and Witten, D. (2019+). Testing for a change in
#' mean after changepoint detection. Technical report.
#'
#' Maidstone, R., Hocking, T., Rigaill, G., & Fearnhead, P. (2017). On optimal
#' multiple changepoint algorithms for large data. Statistics and Computing,
#' 27(2), 519-533.
#'
#' Rigaill, G. (2015). A pruned dynamic programming algorithm to recover the
#' best segmentations with 1 to K_max change-points. Journal de la Societe
#' Francaise de Statistique, 156(4), 180-205.
#'
#' Vostrikova, L. Y. E. (1981). Detecting 'disorder' in multidimensional random
#' processes. In Doklady Akademii Nauk (Vol. 259, No. 2, pp. 270-274). Russian
#' Academy of Sciences.
#'
#' @export
changepoint_estimates <- structure(function(dat, type, tuning_parameter, functional_pruning_out = FALSE) {
  n.data <- length(dat)
  stopifnot(3 <= n.data)
  stopifnot(is.numeric(tuning_parameter))

if (type == "L0") {

stopifnot(length(tuning_parameter) == 1)
  stopifnot(0 <= tuning_parameter)
  cost_mat_r <- double(n.data)
  end_vec_r <- integer(n.data)
  mean_vec_r <- double(n.data)
  intervals_mat_r <- integer(n.data)
  min_mean <- min(dat)
  max_mean <- max(dat)

  fpop_out <- .fpop(dat, tuning_parameter, min_mean, max_mean, cost_mat_r, end_vec_r, mean_vec_r, intervals_mat_r)
  end_vec_r <- end_vec_r +  1

  piecewise_square_losses <- NULL
  if (functional_pruning_out) {
    for (d in fpop_out) {
      piecewise_square_losses <- rbind(piecewise_square_losses, data.frame(d))
    }
    colnames(piecewise_square_losses) <- c("square", "linear", "constant", "min_mean", "max_mean", "prev_last_mean", "data_i", "s")  
  }
  
  out <-
      list(
        estimated_means = rev(mean_vec_r),
        dat = dat,
        type = type, 
        change_pts = rev(unique(end_vec_r[end_vec_r > 0])),
        call = match.call(),
        tuning_parameter = tuning_parameter,
        cost = cost_mat_r,
        n_intervals = intervals_mat_r,
        end_vec = end_vec_r,
        piecewise_square_losses = piecewise_square_losses
      )
    class(out) <- "ChangepointInference_L0_estimated_changes"
    return(out)


  } else if (type == "BS") {
    stopifnot(tuning_parameter > 0)
    
    fit <- .k_step_bs(dat, tuning_parameter) 
    
    aug_changepoints <- c(0, sort(fit[ , 1]), length(dat))
    fitted_means <- 0 * numeric(length(dat))
    for (i in 1:(tuning_parameter + 1)) { 
      fitted_means[(aug_changepoints[i] + 1):aug_changepoints[i + 1]] <- mean(dat[(aug_changepoints[i] + 1):aug_changepoints[i + 1]])
    }
    
    out <-
      list(
        dat = dat,
        type = type, 
        change_pts = sort(fit[ , 1]),
        estimated_means = fitted_means, 
        ordered_change_pts = fit[ , 1],
        change_pt_signs = fit[, 2],
        call = match.call(),
        tuning_parameter = tuning_parameter
      )
    class(out) <- "ChangepointInference_BS_estimated_changes"
    
    return(out)
    
  } else {
    stop(paste0("Type = ", type, " not supported"))
  }
})

