## Class utils 

#' Summary of estimated changepoints and p-values for L0 segmentation
#' 
#' @param  x ChangepointInference_L0_changepoints_pvals object
#'
#' @export
print.ChangepointInference_L0_changepoints_pvals <- function(x, ...) {
  cat("\n Output: \n")
  cat("Type \t\t\t\t\t", x$type, "\n")
  cat("Number of estimated changepoints \t", length(x$change_pts), "\n")
  cat("\n Settings: \n")
  cat("Data length \t\t\t\t", length(x$dat), "\n")
  cat("Tuning parameter \t\t\t", x$tuning_parameter, "\n")
  cat("Window size \t\t\t\t", x$window_size, "\n") 
  if (x$estimated_variance) {
    cat("Estimated variance \t\t\t", x$sig, "\n")   
  } else { 
    cat("Variance \t\t\t\t", x$sig, "\n")   
  }
}

#' Summary of estimated changepoints and p-values for binary segmentation
#' 
#' @param  x ChangepointInference_BS_changepoints_pvals object
#'
#' @export
print.ChangepointInference_BS_changepoints_pvals <- function(x, ...) {
  cat("\n Output: \n")
  cat("Type \t\t\t\t\t", x$type, "\n")
  cat("Number of estimated changepoints \t", length(x$change_pts), "\n")
  cat("\n Settings: \n")
  cat("Data length \t\t\t\t", length(x$dat), "\n")
  cat("Window size \t\t\t\t", x$window_size, "\n") 
  if (x$estimated_variance) {
    cat("Estimated variance \t\t\t", x$sig, "\n")   
  } else { 
    cat("Variance \t\t\t\t", x$sig, "\n")   
  }
}


#' Summary of estimated changepoints for L0 segmentation
#' 
#' @param  x ChangepointInference_L0_estimated_changes object
#'
#' @export
print.ChangepointInference_L0_estimated_changes <- function(x, ...) {
  cat("\n Output: \n")
  cat("Type \t\t\t\t\t", x$type, "\n")
  cat("Number of estimated changepoints \t", length(x$change_pts), "\n")
  cat("\n Settings: \n")
  cat("Data length \t\t\t\t", length(x$dat), "\n")
  cat("Tuning parameter \t\t\t", x$tuning_parameter, "\n")
}

#' Summary of estimated changepoints for binary segmentation
#' 
#' @param  x ChangepointInference_BS_estimated_changes object
#'
#' @export
print.ChangepointInference_BS_estimated_changes <- function(x, ...) {
  cat("\n Output: \n")
  cat("Type \t\t\t\t\t", x$type, "\n")
  cat("Number of estimated changepoints \t", length(x$change_pts), "\n")
  cat("\n Settings: \n")
  cat("Data length \t\t\t\t", length(x$dat), "\n")
}


expand_fpop_intervals <- function(df, PLOT_MIN, PLOT_MAX, ni) {
  n_segments <- dim(df)[[1]]
  df_plot <- NULL
  
  for (seg_i in 1:n_segments) {
    if (df$max_mean[seg_i] >= PLOT_MIN && df$min_mean[seg_i] <= PLOT_MAX) {
      xs <- seq(from = max(PLOT_MIN, df$min_mean[seg_i]), 
                to = min(df$max_mean[seg_i], PLOT_MAX), 
                length.out = ni)
      ys <- df$square[seg_i] * xs ^2 + df$linear[seg_i] * xs + df$constant[seg_i]
      if (!is.null(df$contained)) {
        df_tmp <- data.frame(x = xs, y = ys, contained = df$contained[seg_i])  
      } else { 
        df_tmp <- data.frame(x = xs, y = ys, data_i = df$data_i[seg_i])
      }
      df_plot <- rbind(df_plot, df_tmp)
    }
  }
  return(df_plot)
}

#' Plot optimal cost and/or S as a function of \eqn{\phi}
#' 
#' @param  x ChangepointInference_L0_changepoints_pvals
#' @param thj changepoint 
#' @param plot_cost plot the optimal cost (TRUE), or plot S (FALSE)
#' @param PLOT_MIN minimum phi 
#' @param PLOT_MAX maximum phi 
#' @param ni number of values to calculate the optimal cost at in each segment 
#'
#' @importFrom magrittr %>%
#' @export
plot.ChangepointInference_L0_changepoints_pvals <- function(x, thj, plot_cost = TRUE, PLOT_MIN = -10, PLOT_MAX = 10, ni = 1000, ...) {
  if (is.null(x$conditioning_sets)) { 
    stop("Re-run changepoint_inference with parameter 'return_conditioning_sets' set to true")
    }
  ind <- which(x$change_pts == thj)
  stopifnot(ind > 0)
  
  col_red <- "#d95f02"
  col_blue <- "#1f78b4"
  
  if (plot_cost) {
    df <- x$conditioning_sets[[ind]]
    df_plot <- expand_fpop_intervals(df, PLOT_MIN, PLOT_MAX, ni)
    par(mfrow = c(1, 1))
    plot(df_plot$x, df_plot$y, pch = 20, cex = 0.1,
         col = ifelse(df_plot$contained == 1, col_blue, col_red),
         xlab = latex2exp::TeX("$\\phi$"), 
         ylab = latex2exp::TeX("Cost$(\\phi)$"))
    legend("bottomright", 
           pch = 20,
           col = c(col_blue, col_red),
           c(latex2exp::TeX("$\\phi$ in  S"), latex2exp::TeX("$\\phi$ not in S")))
    
  } else {
    x$conditioning_sets[[ind]] %>% dplyr::mutate(y = 1, contained = factor(contained, levels = c(0, 1), labels = c("Not in S", "In S"))) %>% 
      ggplot2::ggplot() + 
      ggplot2::geom_rect(ggplot2::aes(xmin = min_mean, xmax = max_mean, ymin = -10, ymax = 10, fill = contained)) +
      ggplot2::coord_cartesian(xlim = c(PLOT_MIN, PLOT_MAX)) + 
      ggplot2::xlab(latex2exp::TeX("$\\phi$")) + 
      ggplot2::ylab('') + 
      ggplot2::scale_fill_manual(values=c(col_red, col_blue)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), 
            panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"), 
            legend.position="bottom", legend.title = ggplot2::element_blank())
  }
}

#' Plot S as a function of \eqn{\phi}
#' 
#' @param  x ChangepointInference_BS_changepoints_pvals
#' @param thj changepoint 
#' @param PLOT_MIN minimum phi 
#' @param PLOT_MAX maximum phi 
#'
#' @importFrom magrittr %>%
#' @export
plot.ChangepointInference_BS_changepoints_pvals <- function(x, thj, PLOT_MIN = -10, PLOT_MAX = 10, ...) {
  if (is.null(x$conditioning_sets)) { 
    stop("Re-run changepoint_inference with parameter 'return_conditioning_sets' set to true")
  }
  ind <- which(x$change_pts == thj)
  stopifnot(ind > 0)
  
  col_red <- "#d95f02"
  col_blue <- "#1f78b4"
  
  x$conditioning_sets[[ind]] %>% dplyr::mutate(y = 1, contained = factor(contained, levels = c(0, 1), labels = c("Not in S", "In S"))) %>% 
    ggplot2::ggplot() + 
    ggplot2::geom_rect(ggplot2::aes(xmin = min_mean, xmax = max_mean, ymin = -10, ymax = 10, fill = contained)) +
    ggplot2::coord_cartesian(xlim = c(PLOT_MIN, PLOT_MAX)) + 
    ggplot2::xlab(latex2exp::TeX("$\\phi$")) + 
    ggplot2::ylab('') + 
    ggplot2::scale_fill_manual(values=c(col_red, col_blue)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), 
            panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"), 
            legend.position="bottom", legend.title = ggplot2::element_blank())
}


#' Plot optimal cost Cost_s*(mu) for any s  
#' 
#' @param x ChangepointInference_L0_estimated_changes
#' @param s s in Cost_s*(mu), s = 1,..., length(dat)
#' @param PLOT_MIN minimum phi to calculate cost 
#' @param PLOT_MAX maximum phi to calculate cost 
#' @param ni number of values to calculate the optimal cost at in each segment 
#'
#' @export
plot.ChangepointInference_L0_estimated_changes <- function(x, s = 1, PLOT_MIN = -10, PLOT_MAX = 10, ni = 1000, ...) {
  
  if(is.null(x$piecewise_square_losses)) {
    stop("Rerun after setting functional_pruning_out = TRUE in estimate_changepoints function.")
  }
  
  df <- x$piecewise_square_losses[x$piecewise_square_losses$s == s, ]
  df_plot <- expand_fpop_intervals(df, PLOT_MIN, PLOT_MAX, ni)
  colnames(df_plot) <- c("x", "y", "tau")
  df_plot$tau <- as.factor(df_plot$tau)
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x, y = y, color = tau)) + 
    ggplot2::geom_point(size = 0.5) + 
    ggplot2::xlab(latex2exp::TeX("$\\mu$")) + 
    ggplot2::ylab(latex2exp::TeX(paste0("Cost$_{", s, "}^{*}(\\mu)$"))) + 
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 0)) 
  print(p)
  return(p)
}

eval_cost_at_helper <- function(df, pt) {
  n_segments <- dim(df)[[1]]
  for (seg_i in 1:n_segments) { 
    if (df$min_mean[seg_i] <= pt && df$max_mean[seg_i] >= pt) {
      return(df$square[seg_i] * pt ^2 + df$linear[seg_i] * pt + df$constant[seg_i])
    } 
  }
  warning("Queried point not contained in any interval; returning +Inf.")
  return(Inf) ## pt not contained in any interval 
}



#' Evaluate Cost_s*(mu) for any s and mu 
#' 
#' @param x ChangepointInference_L0_estimated_changes
#' @param s s in Cost_s*(mu), s = 1,..., length(dat)
#' @param mu mean cost to evaluate cost at 
#'
#' @export
eval_cost_at.ChangepointInference_L0_estimated_changes <- function(x, s = 1, mu = 0) {
  df <- x$piecewise_square_losses[x$piecewise_square_losses$s == s, ]
  return(eval_cost_at_helper(df, mu))
}


#' Evaluate Cost(phi) for any phi 
#' 
#' @param x ChangepointInference_L0_changepoints_pvals
#' @param thj changepoint 
#' @param phi evaluate cost at pertubation phi 
#' 
#' @export
eval_cost_at.ChangepointInference_L0_changepoints_pvals <- function(x, thj, phi) {
  ind <- which(x$change_pts == thj)
  stopifnot(ind > 0)
  df <- x$conditioning_sets[[ind]]
  return(eval_cost_at_helper(df, phi))
}

