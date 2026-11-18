# Figure 5: breakpoint in mean and maximum canopy cover along the CHI bins (strucchange)
#
# Intercept-only model (canopy ~ 1), one breakpoint. 
#
# Input:  data/bins_w2.csv
# Output: results/breakpoint_results.csv
#         figures/average_*.tiff, figures/max_*.tiff
#         figures/fig5_panel_average.png, figures/fig5_panel_max.png (Figure 5a, b)

library(strucchange)
library(ggplot2)
library(dplyr)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

bins <- read.csv("data/bins_w2.csv")

# strucchange bin numbers = row numbers
bins$bin_number <- 1:nrow(bins)
print(bins[, c("bin_number", "CH", "CHI_min", "CHI_max", "n_sites", "canopy", "canopy_max")])


run_breakpoint_analysis <- function(bins, y_col, y_label, file_prefix, make_plots = TRUE) {

  y <- bins[[y_col]]
  x <- bins$CH
  n <- length(y)

  canopy_ts <- ts(y)

  # supF test
  fstat     <- Fstats(canopy_ts ~ 1)
  sc_result <- sctest(fstat)

  # breakpoint location and 95% CI (one breakpoint)
  bp    <- breakpoints(canopy_ts ~ 1, breaks = 1)
  ci_bp <- confint(bp, breaks = 1)$confint
  i     <- breakpoints(bp, breaks = 1)$breakpoints

  # BIC: no break vs one break
  bic <- summary(bp)$RSS["BIC", ]

  # strucchange output
  cat("\n==========", y_col, "==========\n")
  print(sc_result)
  print(summary(bp))
  print(confint(bp, breaks = 1))
  cat("Breakpoint at bin number", i, "-> CHI bin", x[i], "\n")
  cat("95% CI from bin number", ci_bp[1], "to", ci_bp[3],
      "-> CHI bins", x[ci_bp[1]], "to", x[ci_bp[3]], "\n")

  # segment means
  mean_before <- mean(y[1:i])
  mean_after  <- mean(y[(i + 1):n])

  # CHI values
  breakpoint <- x[i]
  ci_lower   <- x[ci_bp[1]]
  ci_upper   <- x[ci_bp[3]]

  results <- data.frame(
    response          = y_col,
    n_bins            = n,
    supF              = round(unname(sc_result$statistic), 2),
    p_value           = signif(sc_result$p.value, 3),
    BIC_no_break      = round(unname(bic[1]), 1),
    BIC_one_break     = round(unname(bic[2]), 1),
    bp_bin_number     = i,
    ci_bin_lower      = ci_bp[1],
    ci_bin_upper      = ci_bp[3],
    breakpoint_CH     = breakpoint,
    last_site_before  = round(bins$CHI_max[i], 2),
    first_site_after  = round(bins$CHI_min[i + 1], 2),
    ci_lower          = ci_lower,
    ci_upper          = ci_upper,
    mean_before       = round(mean_before, 1),
    mean_after        = round(mean_after, 1)
  )

  if (make_plots) {

    # F statistic plot
    k <- as.numeric(time(fstat$Fstats))
    fstat_df <- data.frame(CHI = x[k], fstat = as.numeric(fstat$Fstats))

    p_fstat <- ggplot(fstat_df, aes(x = CHI, y = fstat)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = breakpoint, linetype = "dotted", linewidth = 1) +
      labs(x = expression("Cumulative Herbivory Index (g day"^-1*" m"^-2*")"),
           y = "F statistic") +
      theme_classic(base_size = 16) +
      theme(axis.text    = element_text(color = "black"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8))

    ggsave(paste0("figures/", file_prefix, "_fstat_plot.tiff"),
           p_fstat, width = 6, height = 5, dpi = 300)

    # canopy cover plot
    plot_df <- data.frame(CHI = x, canopy = y)

    p_canopy <- ggplot(plot_df, aes(x = CHI, y = canopy)) +
      annotate("rect", xmin = ci_lower, xmax = ci_upper, ymin = -Inf, ymax = Inf,
               fill = "lightgrey", alpha = 0.4) +
      geom_point(size = 2.5) +
      geom_vline(xintercept = breakpoint, linetype = "dotted", linewidth = 1) +
      annotate("segment", x = min(x), xend = breakpoint, y = mean_before, yend = mean_before,
               color = "red", linewidth = 1.2) +
      annotate("segment", x = breakpoint, xend = max(x), y = mean_after, yend = mean_after,
               color = "red", linewidth = 1.2) +
      labs(x = expression("Cumulative Herbivory Index (g day"^-1*" m"^-2*")"),
           y = y_label) +
      coord_cartesian(ylim = c(0, 60)) +
      theme_classic(base_size = 16) +
      theme(axis.text    = element_text(color = "black"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8))

    ggsave(paste0("figures/", file_prefix, "_breakpoint_plot.tiff"),
           p_canopy, width = 6, height = 5, dpi = 300)

    # Figure 5 panel: same axes as panel c (03_gam.R), 140 x 60.8 mm
    # panels a + b + c stacked = 140 x 190 mm
    p_panel <- ggplot(plot_df, aes(x = CHI, y = canopy)) +
      annotate("rect", xmin = ci_lower, xmax = ci_upper, ymin = -Inf, ymax = Inf,
               fill = "lightgrey", alpha = 0.4) +
      geom_point(size = 1.5) +
      geom_vline(xintercept = breakpoint, linetype = "dotted", linewidth = 0.5) +
      annotate("segment", x = min(x), xend = breakpoint, y = mean_before, yend = mean_before,
               color = "red", linewidth = 0.6) +
      annotate("segment", x = breakpoint, xend = max(x), y = mean_after, yend = mean_after,
               color = "red", linewidth = 0.6) +
      scale_x_continuous(breaks = seq(0, 40, 10)) +
      scale_y_continuous(breaks = seq(0, 60, 20)) +
      coord_cartesian(xlim = c(-1, 40), ylim = c(-3, 63), expand = FALSE) +
      labs(x = NULL, y = y_label) +
      theme_classic(base_size = 8) +
      theme(axis.text    = element_text(color = "black"),
            axis.text.x  = element_blank(),
            axis.title   = element_text(face = "bold", size = 9),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
            plot.margin  = margin(2, 3, 2, 2, "mm"))

    ggsave(paste0("figures/fig5_panel_", file_prefix, ".png"),
           p_panel, width = 140, height = 60.8, units = "mm", dpi = 600)
  }

  results
}


# --- Average and maximum canopy cover ------------------------------

res_avg <- run_breakpoint_analysis(bins, y_col = "canopy",
                                   y_label = "Average canopy cover (%)",
                                   file_prefix = "average")

res_max <- run_breakpoint_analysis(bins, y_col = "canopy_max",
                                   y_label = "Maximum canopy cover (%)",
                                   file_prefix = "max")

breakpoint_results <- rbind(res_avg, res_max)
print(breakpoint_results)

write.csv(breakpoint_results, "results/breakpoint_results.csv", row.names = FALSE)

