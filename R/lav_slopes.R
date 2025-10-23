#' Simple slopes and interaction plots for fitted \code{lavaan} models
#'
#' Computes conditional (simple) slopes of a focal predictor across values
#' of a moderator from a fitted \code{lavaan} model that includes their explicit 
#' product term. Plots predicted lines with Wald confidence ribbons, and print 
#' an APA-style test of the interaction for easy reporting and interpretation,
#' plus a simple-slopes table. 
#' 
#' The model should include a main effect for the predictor, a main effect for the moderator,
#' and their product term. The simple slope of the predictor at a given moderator value
#' combines the predictor main effect with the interaction term. The moderator can 
#' be continuous or categorical. Standard errors use the delta method with the 
#' model covariance matrix of the estimates. 
#'
#' @param fit A fitted \code{lavaan} object that includes the product term (required).
#' @param outcome Character. Name of the dependent variable in \code{fit} (required).
#' @param pred Character. Name of the focal predictor whose simple slopes are probed (required).
#' @param modx Character. Name of the moderator (required).
#' @param interaction Character. Name of the product term in \code{fit} (e.g., \code{"X_Z"}) (required).
#' @param data \code{data.frame}. Raw data. If \code{NULL}, the function tries to pull
#' data from \code{fit} via \code{lavInspect}.
#' @param modx.values Numeric or character vector. Values or levels of the moderator
#' at which to compute slopes; derived automatically when \code{NULL}.
#' @param modx.labels Character vector. Legend/table labels for \code{modx.values}
#' (default: the character form of \code{modx.values}).
#' @param pred.range Numeric length-2. Range \code{c(min, max)} for the x-axis;
#' uses observed range in \code{data} when available, else \code{c(-2, 2)}.
#' @param conf.level Numeric in (0,1). Confidence level for CIs and ribbons (default: 0.95).
#' @param x.label Character. X-axis label (default: \code{pred}).
#' @param y.label Character. Y-axis label (default: \code{outcome}).
#' @param legend.title Character. Legend title; if \code{NULL}, the legend shows only levels (default: NULL).
#' @param colors Character vector. Colors for lines and ribbons; named vector recommended with names matching \code{modx.labels} (default: Okabe–Ito palette).
#' @param line.size Numeric > 0. Line width (default: 0.80).
#' @param alpha Numeric in \code{(0,1)}. Ribbon opacity (default 0.20).
#' @param table Logical. Print APA-style interaction test and simple-slopes table (default: \code{TRUE}).
#' @param digits Integer \code{>= 0}. Decimal digits in printed output (default: 2).
#' @param modx_n_unique_cutoff Integer \code{>= 1}. Threshold for treating a numeric moderator
#' as continuous and using mean ± SD (default: 4).
#' @param return_data Logical. If \code{TRUE}, include the plotting data.frame in the returned list (default: FALSE).
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{plot}}{\code{ggplot} object with lines and confidence ribbons.}
#'   \item{\code{slope_table}}{Data frame with moderator levels, simple slopes, SE, z, and CI.}
#'   \item{\code{plot_data}}{Only when \code{return_data = TRUE}: data used to build the plot.}
#' }
#'
#' @section Notes:
#' Estimates are unstandardized; a standardized beta for the interaction is also reported
#' for reference. Wald tests assume large-sample normality of estimates.
#'
#' @examples
#' \dontrun{
#' library(lavaan)
#' # Fit a model with X, Z, and their product X_Z predicting Y
#' # model <- ' Y ~ X + Z + X_Z '
#' # fit <- sem(model, data = your_data)
#' # lav_slopes(fit, outcome = "Y", pred = "X", modx = "Z",
#' #            interaction = "X_Z", legend.title = "Moderator")
#' }
#'
#' @importFrom lavaan parameterEstimates parTable lavInspect
#' @importFrom stats vcov qnorm sd pnorm
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon scale_colour_manual scale_fill_manual labs
#' @export
#' 
#' @usage
#' lav_slopes(
#'   fit,
#'   outcome,
#'   pred,
#'   modx,
#'   interaction,
#'   data = NULL,
#'   modx.values = NULL,
#'   modx.labels = NULL,
#'   pred.range = NULL,
#'   conf.level = 0.95,
#'   x.label = NULL,
#'   y.label = NULL,
#'   legend.title = NULL,
#'   colors = NULL,
#'   line.size = 0.80,
#'   alpha = 0.20,
#'   table = TRUE,
#'   digits = 2,
#'   modx_n_unique_cutoff = 4L,
#'   return_data = FALSE
#' )

lav_slopes <- function(
    fit,
    outcome,
    pred,
    modx,
    interaction,
    data = NULL,
    modx.values = NULL,
    modx.labels = NULL,
    pred.range = NULL,
    conf.level = 0.95,
    x.label = NULL,
    y.label = NULL,
    legend.title = NULL,
    colors = NULL,
    line.size = 0.80,
    alpha = 0.20,
    table = TRUE,
    digits = 2,
    modx_n_unique_cutoff = 4L,
    return_data = FALSE
) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  
  if (!inherits(fit, "lavaan")) stop("'fit' must be a lavaan object.", call. = FALSE)
  if (missing(interaction) || is.null(interaction))
    stop("'interaction' (product-term name) must be supplied.", call. = FALSE)
  
# cerca di recuperare i dati!!! perchè non funziona?
  fetch_data <- function(f) {
    out <- tryCatch(lavaan::lavInspect(f, "data"), error = function(e) NULL)
    if (is.null(out))
      out <- tryCatch(lavaan::lavInspect(f, "data.original"), error = function(e) NULL)
    if (is.null(out))
      out <- tryCatch(as.data.frame(f@Data@X), error = function(e) NULL)
    out
  }
  dat <- data %||% fetch_data(fit)
  
# estrai valori moderatori e label 
  if (is.null(modx.values)) {
    if (is.null(dat) || !(modx %in% names(dat)))
      stop("'modx.values' missing and moderator data unavailable.\n",
           "Supply 'modx.values' or pass raw data via 'data'.", call. = FALSE)
    z <- dat[[modx]]
    if (is.numeric(z) && length(unique(z)) > modx_n_unique_cutoff) {
      m  <- mean(z, na.rm = TRUE); sd <- stats::sd(z, na.rm = TRUE)
      modx.values <- round(c(m - sd, m, m + sd), 2)
      modx.labels <- modx.labels %||% c("-1 SD", "Mean", "+1 SD")
      message("Derived modx.values as Mean +/- 1 SD: ",
              paste(modx.values, collapse = ", "))
    } else {
      modx.values <- sort(unique(z))
      modx.labels <- modx.labels %||% as.character(modx.values)
      message("Derived modx.values from moderator levels: ",
              paste(modx.labels, collapse = ", "))
    }
  }
  if (!is.null(modx.labels) && length(modx.labels) != length(modx.values))
    stop("'modx.labels' must match length of 'modx.values'.", call. = FALSE)
  modx.labels <- modx.labels %||% as.character(modx.values)
  
# range del predittore focale
  if (is.null(pred.range)) {
    if (!is.null(dat) && pred %in% names(dat))
      pred.range <- range(dat[[pred]], na.rm = TRUE)
    else
      pred.range <- c(-2, 2)
  }
  x_seq <- seq(pred.range[1L], pred.range[2L], length.out = 1000L)
  
# stime dei parametri e vcov
  pe <- lavaan::parameterEstimates(fit, standardized = FALSE)
  pt <- lavaan::parTable(fit)
  vc <- tryCatch(stats::vcov(fit),
                 error = function(e) lavaan::lavInspect(fit, "vcov"))
  if (!is.matrix(vc)) stop("Could not retrieve covariance matrix from 'fit'.", call. = FALSE)
  
  find_par <- function(lhs, rhs, op = "~", rhs_regex = NULL,
                       required = TRUE, default_est = 0) {
    idx_pe <- which(pe$lhs == lhs & pe$rhs == rhs & pe$op == op)
    if (length(idx_pe) == 0L && !is.null(rhs_regex))
      idx_pe <- which(pe$lhs == lhs & pe$op == op & grepl(rhs_regex, pe$rhs))
    if (length(idx_pe) == 0L) {
      if (!required) { est_val <- default_est; idx_pe <- NA_integer_ } else {
        stop(sprintf("Parameter '%s %s %s' not found.", lhs, op, rhs), call. = FALSE)
      }
    } else {
      idx_pe  <- idx_pe[1L]; est_val <- pe$est[idx_pe]
    }
    idx_pt <- which(pt$lhs == lhs & pt$rhs == rhs & pt$op == op)
    if (length(idx_pt) == 0L && !is.null(rhs_regex))
      idx_pt <- which(pt$lhs == lhs & pt$op == op & grepl(rhs_regex, pt$rhs))
    free_idx <- if (length(idx_pt) == 0L) 0L else {
      fi <- pt$free[idx_pt[1L]]; if (is.na(fi)) 0L else as.integer(fi)
    }
    list(est = est_val, free = free_idx, idx_pe = idx_pe)
  }
  
# estrai inetercetta + effetti principali + interazione
  p0 <- find_par(outcome, "",   op = "~1", required = FALSE, default_est = 0)
  p1 <- find_par(outcome, pred, op = "~")
  p2 <- find_par(outcome, modx, op = "~")
  p3 <- find_par(outcome, interaction, op = "~",
                 rhs_regex = paste0("(", pred, ").*(", modx, ")|(",
                                    modx, ").*(", pred, ")"))
  
  safe_free <- function(x) { if (is.na(x) || length(x) == 0L) 0L else as.integer(x) }
  p0$free <- safe_free(p0$free); p1$free <- safe_free(p1$free)
  p2$free <- safe_free(p2$free); p3$free <- safe_free(p3$free)
  
  vc_ <- function(a, b) {
    a <- safe_free(a); b <- safe_free(b)
    if (a <= 0L || b <= 0L) return(0)
    if (a > nrow(vc) || b > ncol(vc)) return(0)
    val <- vc[a, b]; if (is.na(val)) 0 else val
  }
  
  var_b1   <- vc_(p1$free, p1$free)
  var_b3   <- vc_(p3$free, p3$free)
  cov_b1b3 <- vc_(p1$free, p3$free)
  
  var_b0   <- vc_(p0$free, p0$free)
  var_b2   <- vc_(p2$free, p2$free)
  cov_b0b2 <- vc_(p0$free, p2$free)
  cov_b0b1 <- vc_(p0$free, p1$free)
  cov_b0b3 <- vc_(p0$free, p3$free)
  cov_b1b2 <- vc_(p1$free, p2$free)
  cov_b2b3 <- vc_(p2$free, p3$free)
  
# slope e anche predizioni
  z_crit <- stats::qnorm(1 - (1 - conf.level) / 2)
  slope_tbl <- data.frame(
    Moderator = modx.labels,
    Z_value   = modx.values,
    Slope = NA_real_, SE = NA_real_, z = NA_real_,
    CI_low = NA_real_, CI_high = NA_real_,
    stringsAsFactors = FALSE
  )
  plot_df <- NULL
  
  for (i in seq_along(modx.values)) {
    z0 <- modx.values[i]
    
    slope     <- p1$est + p3$est * z0
    var_slope <- var_b1 + 2 * z0 * cov_b1b3 + (z0 ^ 2) * var_b3
    se_slope  <- sqrt(max(var_slope, 0))
    
    slope_tbl$Slope[i]   <- slope
    slope_tbl$SE[i]      <- se_slope
    slope_tbl$z[i]       <- if (se_slope > 0) slope / se_slope else NA_real_
    slope_tbl$CI_low[i]  <- slope - z_crit * se_slope
    slope_tbl$CI_high[i] <- slope + z_crit * se_slope
    
    intercept <- p0$est + p2$est * z0
    var_int   <- max(var_b0 + 2 * z0 * cov_b0b2 + (z0 ^ 2) * var_b2, 0)
    
    y_hat   <- intercept + slope * x_seq
    se_yhat <- sqrt(pmax((x_seq ^ 2) * var_slope + var_int +
                           2 * x_seq * (cov_b0b1 + z0 * cov_b0b3 +
                                          z0 * cov_b1b2 + (z0 ^ 2) * cov_b2b3), 0))
    
    plot_df <- rbind(plot_df,
                     data.frame(
                       X = x_seq, Y_hat = y_hat,
                       CI_low = y_hat - z_crit * se_yhat,
                       CI_high= y_hat + z_crit * se_yhat,
                       LegendGroup = factor(modx.labels[i], levels = modx.labels),
                       stringsAsFactors = FALSE
                     ))
  }
  
# colori okabe ito
  if (is.null(colors)) {
    okabe_ito <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    colors <- okabe_ito[seq_along(modx.labels)]
  }
  if (is.null(names(colors)) || !all(modx.labels %in% names(colors)))
    names(colors) <- modx.labels
  
# plotting con ggplot 
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = X, y = Y_hat, colour = LegendGroup, fill = LegendGroup)) +
    ggplot2::geom_line(linewidth = line.size) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = CI_low, ymax = CI_high),
      alpha = alpha, colour = NA) +
    ggplot2::scale_colour_manual(values = colors, name = legend.title) +
    ggplot2::scale_fill_manual(values = colors, name = legend.title) +
    ggplot2::labs(x = x.label %||% pred, y = y.label %||% outcome)
  
# output in console
  if (table) {
    int_se <- if (!is.na(p3$idx_pe)) pe$se[p3$idx_pe] else NA_real_
    int_z  <- if (!is.na(int_se) && int_se > 0) p3$est / int_se else NA_real_
    int_p  <- if (!is.na(int_z)) 2 * stats::pnorm(abs(int_z), lower.tail = FALSE) else NA_real_
    int_ci <- c(p3$est - z_crit * int_se, p3$est + z_crit * int_se)
    
    pe_std <- lavaan::parameterEstimates(fit, standardized = TRUE)
    idx_std <- which(pe_std$lhs == outcome & pe_std$rhs == interaction & pe_std$op == "~")
    beta_std <- if (length(idx_std)) pe_std$std.all[idx_std[1L]] else NA_real_
    
    fmt  <- function(x, k = digits) if (is.na(x)) "NA" else formatC(x, format = "f", digits = k)
    fmtp <- function(p) ifelse(is.na(p), "NA", ifelse(p < .001, "< .001",
                                                      formatC(p, digits = 3, format = "f")))
    
    cat("\nInteraction effect (", pred, " * ", modx, " → ", outcome, "): ",
        "b = ", fmt(p3$est), ", SE = ", fmt(int_se),
        ", β = ", fmt(beta_std), ", z = ", if (is.na(int_z)) "NA" else fmt(int_z, k = max(1, digits+1)),
        ", p = ", fmtp(int_p),
        ", ", sprintf("%.0f%% CI", conf.level*100), " [", fmt(int_ci[1L]), ", ", fmt(int_ci[2L]), "]",
        "\n", sep = "")
  }
  
  if (table) {
    tbl <- slope_tbl
    fmt  <- function(x, k = digits) formatC(x, format = "f", digits = k)
    cat("\nSimple slopes of", pred, "predicting", outcome,
        "at levels of", modx, sprintf("(%.0f%% CI)", conf.level*100), "\n", sep = " ")
    cat(rep("-", 72), "\n", sep = "")
    cat(sprintf("%-18s %-10s %10s %8s %8s %13s\n",
                "Moderator", "Z_value", "Slope", "SE", "z",
                sprintf("%.0f%% CI", conf.level*100)))
    cat(rep("-", 72), "\n", sep = "")
    for (i in seq_len(nrow(tbl))) {
      cat(sprintf("%-18s %-10s %10s %8s %8s [%s, %s]\n",
                  as.character(tbl$Moderator[i]),
                  as.character(tbl$Z_value[i]),
                  fmt(tbl$Slope[i]),
                  fmt(tbl$SE[i]),
                  if (is.na(tbl$z[i])) "NA" else fmt(tbl$z[i], k = max(1, digits+1)),
                  fmt(tbl$CI_low[i]),
                  fmt(tbl$CI_high[i])
      ))
    }
    cat(rep("-", 72), "\n\n", sep = "")
  }
  
  print(p)
  out <- list(plot = p, slope_table = slope_tbl)
  if (return_data) out$plot_data <- plot_df
  invisible(out)
}
