#' Estimation and inference for changepoints estimated via L0 and binary
#' segmentation
#'
#' @param dat observed data
#' @param type specifies algorithm and conditioning event; see details.
#' @param tuning_parameter L0 segmentation: tuning parameter \eqn{\lambda};
#'   Binary segmentation: number of segments K.
#' @param window_size window size for fixed window hypothesis testing.
#' @param sig noise variance. If unknown, the sample variance is used for to
#'   compute p-values.
#' @param approximation_threshold For binary segmentation, this threshold
#'   specifies how accurately S is calculated. In particular, the set S is
#'   calculated exactly up to +/- max( approximation_threshold * || v ||_2
#'   sqrt(sig), |sum(v * dat)| ) for a contrast v depending on inference type.
#'
#' @return Returns a list with elements:
#' @return \code{change_pts} the set of changepoints
#' @return \code{pvals}  p-values associated with each changepoint
#' @return \code{pval_approximation_error} constants C such that |pval - true_p|
#'   <= C where true_p is the oracle p-value obtained with
#'   approximation_threshold = Inf
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
#' Inference:
#'
#' For each estimated changepoint \eqn{\hat{\tau}_j}, we test the null
#' hypothesis that there is no change in mean *near* \eqn{\hat{\tau}_j}. In
#' particular, for a linear contrast \eqn{v in R^T} defined as \eqn{v_t = 0}, if
#' \eqn{t < tL} or \eqn{t > tR}, \eqn{v_t = 1 / (\hat{\tau}_j - tL + 1)}, if
#' \eqn{tL \le t \le \hat{\tau}_j}, and \eqn{v_t = -1 / (tR - \hat{\tau}_j)}, if
#' \eqn{\hat{\tau}_j + 1 \le t \le tR}, we test the null hypothesis \eqn{H0:
#' \sum_t v_t * \mu_t = 0}.
#'
#' We define adaptive and fixed defintions of *near*. In the adaptive case, near
#' is defined based on the neighboring changepoints: \deqn{tL = max(1,
#' \hat{\tau}_{j-1} + 1),} and \deqn{tR = min(T, \hat{\tau}_{j + 1}).} In the
#' fixed window case, near is defined based on a fixed window around
#' \eqn{\hat{\tau}_j}: \deqn{tL = max(1, \hat{\tau}_j - window_size + 1),} and
#' \deqn{tR = min(T, \hat{\tau}_j + window_size).}
#'
#' Our framework allows us to efficiently compute p-values based on adaptive and
#' fixed test statistics and for different conditioning sets. Let \eqn{t(v) * y =
#' \sum_t v_t * y_t} be the observed test statistic and let \eqn{\phi = \sum_t
#' v_t * Y_t}. In Jewell et al. (2019+), we show that the p-value corresponding
#' to the test \eqn{H0: \sum_t v_t * \mu_t = 0} can be written as
#' \deqn{Pr(|\phi| \ge |t(v) * y| | \phi in S)} for a conditioning set S.
#'
#' This software computes p-values for the following test statistics and
#' conditioning sets S. In what follows, \eqn{y'(\phi)} is a perturbation of the
#' observed data y. See Theorem 1 of Jewell et al. (2019+) for additional
#' details, and \code{\link{construct_perturbed_data_tL_tR}}.
#'
#' Type = 'L0-fixed': \deqn{Pr(|\phi| >= |t(v) * y| | \hat{\tau}_j in
#' M(y'(\phi))), for fixed v.}
#'
#' Type = 'BS-fixed': \deqn{Pr(|\phi| >= |t(v) * y| | \hat{\tau}_j in
#' M(y'(\phi))), for fixed v.}
#'
#' Type = 'BS-adaptive-M-O-D':
#'
#' \deqn{Pr(|\phi| >= |t(v) * y| | M(y) = M(y'(\phi)), O(y) = O(y'(\phi)),
#' \Delta(y) = \Delta(y'(\phi)), for adaptive v.}
#'
#' Type = 'BS-adaptive-M-O':
#'
#' \deqn{Pr(|\phi| >= |t(v) * y| | M(y) = M(y'(\phi)), O(y) = O(y'(\phi))), for
#' adaptive v.}
#'
#' Type = 'BS-adaptive-M':
#'
#' \deqn{Pr(|\phi| >= |t(v) * y| | M(y) = M(y'(\phi))), for adaptive v.}
#'
#' Since Y_t is Gaussian, \eqn{\phi | S} is a Gaussian truncated to S.
#' Therefore, to calculate each of these probabilities, we must first determine
#' S. For L0 segmentation, we compute S exactly; for binary segmentation we
#' compute S exactly up to +/- max( approximation_threshold * || v ||_2
#' sqrt(sig), |t(v) * y| ). The p-value using this approximation is a
#' conservative estimate of the true p-value; we also output a constant C such
#' that \eqn{|pval - true_p| \le C} where true_p is the p-value obtained with
#' approximation_threshold = Inf. If the noise variance is unknown, we use the
#' sample variance based on the segmentation of the data.
#' @examples
#'
#' ### Generate sample data
#' set.seed(1)
#' mu <- rep(c(1, 3, -4), each = 100)
#' dat <- mu + rnorm(length(mu))
#'
#' ### Basic usage
#'
#' ## L0 segmentation with fixed window size
#' fit <- changepoint_inference(dat, "L0-fixed", 4, window_size = 10, sig = 1)
#' print(fit)
#'
#' # Estimated changepoints and p-values
#' data.frame(change_pts = fit$change_pts,pvals = fit$pvals)
#'
#' ## Binary segmentation
#' k <- 2
#'
#' # with fixed window size, and conditioning on just tau-hat-j in M(y'(phi))
#' fit <- changepoint_inference(dat, "BS-fixed", k, window_size = 10, sig = 1)
#' print(fit)
#' data.frame(change_pts = fit$change_pts,pvals = fit$pvals)
#'
#' # with adaptive window size, and conditioning on the estimated changepoints
#' fit <- changepoint_inference(dat, "BS-adaptive-M", k, sig = 1)
#' print(fit)
#' data.frame(change_pts = fit$change_pts,pvals = fit$pvals)
#'
#'
#' # with adaptive window size, and conditioning on the estimated changepoints, orders, and signs
#' fit <- changepoint_inference(dat, "BS-adaptive-M-O-D", k, sig = 1)
#' print(fit)
#' data.frame(change_pts = fit$change_pts,pvals = fit$pvals)
#'
#'
#' ### More advanced features
#'
#' ## L0 segmentation with fixed window size
#' ## we now set return_conditioning_sets = TRUE
#' fit <- changepoint_inference(dat, "L0-fixed", 4, window_size = 10, sig = 1, return_conditioning_sets = TRUE)
#' print(fit)
#'
#' ## Examine the set S for the first estimated changepoint
#' ## Plot the optimal cost of segmenting y'(phi)
#' ## where colors indicate whether or not phi is in S.
#' plot(fit, fit$change_pts[1])
#'
#' ## Just plot the set S for the first changepoint
#' plot(fit, fit$change_pts[1], FALSE)
#'
#' ## The sets S for all changepoints are obtained through the following list
#' ## The ith element of the list corresponds to the ith changepoint in fit$change_pts
#' ## Each row corresponds to the optimal cost (as a quadratic with square, linear, and constant coefficients)
#' ## of segmenting y'(phi) for phi in [min_mean, max_mean]. Contained = 1 denotes that the interval is contained in S
#' fit$conditioning_sets
#'
#' ## Similar functions for BS
#' fit <- changepoint_inference(dat, "BS-fixed", k, window_size = 10, sig = 1, return_conditioning_sets = TRUE)
#' print(fit)
#'
#' ## Again, the sets S for all changepoints are obtained through the following list
#' ## The ith element of the list corresponds to the ith changepoint in fit$change_pt
#' ## Each row corresponds to whether or not the interval [min_mean, max_mean] is contained in S (contained = 1)
#' ## or not (contained = 0)
#' fit$conditioning_sets
#'
#' ## Plot the set S for the first changepoint
#' plot(fit, fit$change_pts[1])
#'
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
changepoint_inference <- structure(function(dat, type, tuning_parameter, window_size = NULL, sig = NULL, approximation_threshold = 10, return_conditioning_sets = FALSE) {
  stopifnot(tuning_parameter > 0)
  stopifnot(approximation_threshold > 1)
  stopifnot(type %in% c("L0-fixed", "BS-fixed", "BS-adaptive-M-O-D", "BS-adaptive-M-O", "BS-adaptive-M"))
  ## L0 segmentation
  if (is.null(sig)) {
    estimated = TRUE
  } else {
    estimated = FALSE
  }
  if (type == "L0-fixed") {
    if (is.null(sig)) {
      fit <- changepoint_estimates(dat, "L0", tuning_parameter)
      sig <- var(dat - fit$estimated_means)
    }
    
    # warning("Approximation threshold not yet implemented for L0 segmentation. You can safely ignore this error.")
    out_fpop_inference <- .fpop_inference(dat, tuning_parameter, window_size, sig, return_conditioning_sets)

    if (return_conditioning_sets) { 
      conditioning_sets = fpop_inference_intervals_formatter(out_fpop_inference[[2]])
    } else { 
      conditioning_sets = NA
    }
    
    out <- list(
    dat = dat,
    type = type,
    call = match.call(),
    tuning_parameter = tuning_parameter,
    window_size = window_size,
    sig = sig,
    change_pts = out_fpop_inference[[1]][, 1],
    pvals = out_fpop_inference[[1]][, 2],
    approximation_error = out_fpop_inference[[1]][, 3],
    conditioning_sets = conditioning_sets, 
    estimated_variance = estimated
    )
    class(out) <- "ChangepointInference_L0_changepoints_pvals"
    return(out)
  } else {
    ## Binary segmentation
    if (is.null(sig)) {
      fit <- changepoint_estimates(dat, "BS", tuning_parameter)
      sig <- var(dat - fit$estimated_means)
    }
      if (type == "BS-fixed" && window_size <= 0) stop("BS-fixed requires window size > 0.")

      if (type == "BS-fixed") {
        type_cpp <- 0
      } else if (type == "BS-adaptive-M-O-D") {
        type_cpp <- 1
      } else if (type == "BS-adaptive-M-O") {
        type_cpp <- 2
      } else if (type == "BS-adaptive-M") {
        type_cpp <- 3
      } else {
        stop("invalid type")
      }
    
      if (is.null(window_size)) { 
        window_size <- 0
      }  
    
      out_bs_inference <- .k_step_bs_inference(dat, tuning_parameter, type_cpp, sig, window_size, approximation_threshold, return_conditioning_sets)

      if (return_conditioning_sets) { 
        conditioning_sets = bs_inference_intervals_formatter(out_bs_inference[[2]])
      } else { 
        conditioning_sets = NA
      }
      
      out <- list(
      dat = dat,
      type = type,
      call = match.call(),
      tuning_parameter = tuning_parameter,
      window_size = window_size,
      sig = sig,
      change_pts = out_bs_inference[[1]][, 1],
      pvals = out_bs_inference[[1]][, 2],
      approximation_error = out_bs_inference[[1]][, 3], 
      ordered_change_pts = out_bs_inference[[3]],
      conditioning_sets = conditioning_sets, 
      estimated_variance = estimated
      )

      class(out) <- "ChangepointInference_BS_changepoints_pvals"
      return(out)

  }
})
