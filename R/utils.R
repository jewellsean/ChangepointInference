#' Construct contrast around an estimated changepoint
#'
#' @param tL left end point
#' @param thj estimated changepoint to construct contrast vector around
#' @param tR right end point
#' @param n number of observed datapoints
#'
#' @return \code{nu} contrast vector
#'
#' @details
#'
#' The ith entry of the n-vector contrast is defined as:
#'
#' 0,                    i < tL or i > tR
#'
#' 1 / \eqn{(\hat{\tau}_j - tL + 1)},    \eqn{tL \le i \le \hat{\tau}_j}
#'
#' -1 / \eqn{tR - \hat{\tau}_j},       \eqn{\hat{\tau}_j + 1 \le i \le tR}
#'
#' where tL = max(1, \eqn{\hat{\tau}}_j - h + 1) and tR = min(T,
#' \eqn{\hat{\tau}_j} + h) in the fixed window case with window size = h, and
#' where tL = max(1, \eqn{\hat{\tau}_{j-1}} + 1) and tR = min(T,
#' \eqn{\hat{\tau}_{j + 1}}) in the adaptive case.
#'
#'
#' @examples
#'
#' ### Generate sample data
#' set.seed(1)
#' mu <- rep(c(1, 3, -4), each = 100)
#' dat <- mu + rnorm(length(mu))
#'
#' # Suppose a changepoint is estimated at thj = 100. Then a fixed window contrast around
#' # thj with window size = 10 is constructed as follows:
#' thj <- 100
#' window_size <- 10
#' v <- construct_v_tL_tR(thj - window_size + 1, thj, thj + window_size, length(dat))
#'
#' # The change in mean using the underlying mean around thj is (1 - 3 = -2)
#' sum(v * mu)
#'
#' # Using the observed data
#' sum(v * dat)
#'
#' @seealso 
#' \code{\link{changepoint_estimates}}, \code{\link{changepoint_inference}}
#'
#' @references 
#' Jewell, S., Fearnhead, P., and Witten, D. (2019+). Testing for a change in mean after changepoint detection. Technical report.
#'
#' @export
construct_v_tL_tR <- function(tL, thj, tR, n) {
  ind1 <- tL:thj
  ind2 <- (thj+1):tR
  
  v <- numeric(n)
  v[ind1] <- 1 / (thj - tL + 1)
  v[ind2] <- -1 / (tR - thj)
  return(v)
}


#' Construct perturbed data around an estimated changepoint
#'
#' @param y observed data
#' @param tL left end point
#' @param thj estimated changepoint
#' @param tR right end point
#' @param phi perturb constant
#'
#' @return \code{yphi} perturbed data
#'
#' @details
#'
#' Perturbs the observed data y around thj based on a scalar phi:
#'
#' \deqn{y'(\phi) = y - v * t(v) * y / ||v||_2^2 + v * \phi / ||v||_2^2}
#'
#' See \code{\link{construct_v_tL_tR}} for description of tL and tR.
#'
#' @examples
#'
#' ### Generate sample data
#' set.seed(1)
#' mu <- rep(c(1, 3, -4), each = 100)
#' dat <- mu + rnorm(length(mu))
#'
#' # Suppose a changepoint is estimated at thj = 100. Then a fixed window contrast around
#' # thj with window size = 10 is constructed as follows:
#' thj <- 100
#' window_size <- 10
#' v <- construct_v_tL_tR(thj - window_size + 1, thj, thj + window_size, length(dat))
#'
#' # Perturb the data around thj
#' phi <- 2
#' yphi <- dat - v %*% t(v) %*% dat / (sum(v * v)) + v * phi / sum(v * v)
#'
#' # construct_perturbed_data_tL_tR is just a utility function for the above calculation:
#' sum(yphi - construct_perturbed_data_tL_tR(dat, thj - window_size + 1, thj, thj + window_size, phi))
#'
#' @seealso \code{\link{changepoint_estimates}},
#' \code{\link{changepoint_inference}}, \code{\link{construct_v_tL_tR}}
#'
#' @references Jewell, S., Fearnhead, P., and Witten, D. (2019+). Testing for a
#' change in mean after changepoint detection. Technical report.
#'
#' @export
construct_perturbed_data_tL_tR <- function(y, tL, thj, tR, phi) {
  n <- length(y)
  yphi <- y
  nu <- construct_v_tL_tR(tL, thj, tR, n)
  ind1 <- tL:thj
  ind2 <- (thj+1):tR
  n1 <- length(ind1)
  n2 <- length(ind2)
  yphi[ind1] <- yphi[ind1] + (phi - as.numeric(t(nu) %*% y)) / (1 + n1 / n2)
  yphi[ind2] <- yphi[ind2] - (phi - as.numeric(t(nu) %*% y)) / (1 + n2 / n1)
  return(yphi)
}


### 
### Formatting utilities
###
fpop_inference_intervals_formatter <- function(phi_list) {
  phi_outs <- list()
  for (i in 1:length(phi_list)) {
    phi_i <- phi_list[[i]]
    phi_out <- data.frame(phi_i)
    colnames(phi_out) <- c("square", "linear", "constant", "min_mean", "max_mean", "prev_last_mean", "contained", "s")
    phi_out$contained <- phi_out$contained - 1
    phi_outs[[i]] <- phi_out[,c("square", "linear", "constant", "min_mean", "max_mean", "contained")]
  }
  return(phi_outs) 
}


bs_inference_intervals_formatter <- function(phi_list) {
  phi_outs <- list()
  for (i in 1:length(phi_list)) {
    phi_i <- phi_list[[i]]
    phi_out <- data.frame(phi_i)
    colnames(phi_out) <- c("min_mean", "max_mean", "contained")
    phi_outs[[i]] <- phi_out
  }
  return(phi_outs) 
}