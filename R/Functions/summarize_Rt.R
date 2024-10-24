summarize_Rt <- function(rt_samples, bci_level = 0.95,
                         true_rt = NULL, plot = FALSE,
                         plot_title = NULL,
                         ylim = c(0, 3)) {
  # True Rt same length as estimated Rt
  if (!is.null(true_rt)) stopifnot(length(true_rt) == ncol(rt_samples))
  duration <- seq_len(ncol(rt_samples))

  # Get quantiles
  quantiles <- sort(c(0.5, (1 - bci_level) / 2, (1 + bci_level) / 2))
  rt_quantiles <- matrix(ncol = ncol(rt_samples), nrow = length(quantiles))
  rownames(rt_quantiles) <- paste0("quantile_", quantiles)
  for (i in 1:3) {
    rt_quantiles[i, ] <- apply(rt_samples, MARGIN = 2,
                               FUN = function(samples) {
                                 quantile(samples, probs = quantiles[i])
                               })
  }

  # Get summary table
  MCIW <- round(mean(abs(rt_quantiles[1, ] - rt_quantiles[3, ])), 3)
  MASV <- round(mean(abs(diff(rt_quantiles[2, ]))), 3)
  summary_tbl <- data.frame(MCIW = MCIW, MASV = MASV)
  if (!is.null(true_rt)) {
    summary_tbl$MAD <- round(mean(abs(rt_quantiles[2, ] - true_rt)),3)
    summary_tbl$coverage <- round(mean(rt_quantiles[1, ] < true_rt &
                                         true_rt < rt_quantiles[3, ]), 3)
    summary_tbl$true_MASV <- round(mean(abs(diff(true_rt))), 3)
  }
  rownames(summary_tbl) <- plot_title

  if (plot) {
    par(mfrow = c(1, 1),
        mar = c(4, 6, 3, 1))
    plot.new()
    plot.window(xlim = c(min(duration), max(duration)), ylim = ylim)
    axis(1, at = duration, cex.axis = 1.5); axis(2, cex.axis = 1.5)
    polygon(x = c(duration, rev(duration)), y = c(rt_quantiles[1, ],
                                                  rev(rt_quantiles[3, ])),
            col = "slategray2", border = "white")
    lines(duration, rt_quantiles[1, ], lty = "dashed", lwd = 1.5)
    lines(duration, rt_quantiles[2, ], lwd = 1.5)
    lines(duration, rt_quantiles[3, ], lty = "dashed", lwd = 1.5)
    if (!is.null(true_rt)) lines(duration, true_rt, col = "red", lwd = 1.5)
    title(plot_title,
          xlab = "Week",
          ylab = "Effective Reproduction\n Number",
          cex.main = 1.5,
          cex.lab = 1.5)
  }
  return(summary_tbl)
}