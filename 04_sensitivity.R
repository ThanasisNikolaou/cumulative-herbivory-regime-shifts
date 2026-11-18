# Sensitivity of the CHI and of the threshold (Supplementary)
#
# Part 1. Monte Carlo: consumption-rate parameters drawn within the range reported in each source
# Part 2. Alternative indices: (a) equal grazing rate, (b) untransformed F, (c) leave one species out
# Part 3. Binning: (a) bin width and starting point, (b) no merging (2 sites rule), (c) equal-count bins, (d) no binning
#
# In Parts 1 and 2 the bin width is scaled to the range of each index (2 x max / baseline max),
# and the breakpoint location is given as the share of sites below the break.
#
# Input:  data/dataset_eastMed_reefs_abundance_density.csv (through 00_CHI.R)
# Output: results/mc_draws.csv, mc_summary.csv, mc_site_CHI.csv,
#         alternative_indices.csv, sensitivity_binning.csv
#         figures/supplementary figures/FigS3_MC_site_CHI.png

library(strucchange)
library(mgcv)
library(dplyr)
library(ggplot2)

source("00_CHI.R")    # parameters, consumption_rates(), calc_chi()
source("01_bins.R")   # make_bins()

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

sites <- read.csv("data/sites_CHI.csv")
n_sites_total <- nrow(sites)


# --- Functions -----------------------------------------------------

# breakpoint on a set of bins (h = 15% of bins; 2 bins if < 14 bins)
breakpoint_summary <- function(bins, y_col) {

  y <- ts(bins[[y_col]])
  n <- length(y)
  h <- if (n >= 14) floor(0.15 * n) else 2

  bp  <- breakpoints(y ~ 1, breaks = 1, h = h)
  i   <- breakpoints(bp, breaks = 1)$breakpoints
  ci  <- suppressWarnings(confint(bp, breaks = 1)$confint)
  bic <- summary(bp)$RSS["BIC", ]

  # CH of bin k (NA if outside the bins)
  ch_at <- function(k) if (is.na(k) || k < 1 || k > n) NA else bins$CH[k]

  data.frame(
    response          = y_col,
    n_bins            = n,
    h                 = h,
    bp_bin_number     = i,
    breakpoint_CH     = ch_at(i),
    ci_lower_CH       = ch_at(ci[1]),
    ci_upper_CH       = ch_at(ci[3]),
    sites_below       = sum(bins$n_sites[1:i]),
    share_below       = round(sum(bins$n_sites[1:i]) / n_sites_total, 3),
    BIC_prefers_break = unname(bic[2] < bic[1]),
    p_value           = signif(sctest(Fstats(y ~ 1, from = h))$p.value, 3)
  )
}

# bin width scaled to the range of the index
scaled_bins <- function(s) make_bins(s, width = 2 * max(s$CHI) / max(sites$CHI))

# site-level GAM, as in 03_gam.R
fit_gam <- function(s) gam(canopy_avg ~ s(CHI, k = 4), family = tw(link = "log"),
                           data = s, method = "REML")

# breakpoints and GAM for one index
run_index <- function(index_values) {
  s <- sites
  s$CHI <- index_values
  b <- scaled_bins(s)
  g <- fit_gam(s)
  bp_avg <- breakpoint_summary(b, "canopy")
  bp_max <- breakpoint_summary(b, "canopy_max")

  data.frame(
    rank_corr       = round(cor(index_values, sites$CHI, method = "spearman"), 3),
    n_bins          = bp_avg$n_bins,
    avg_breakpoint  = round(bp_avg$breakpoint_CH, 2),
    avg_share_below = bp_avg$share_below,
    avg_p           = bp_avg$p_value,
    max_breakpoint  = round(bp_max$breakpoint_CH, 2),
    max_share_below = bp_max$share_below,
    max_p           = bp_max$p_value,
    gam_p           = signif(summary(g)$s.table[1, "p-value"], 3),
    gam_dev_expl    = round(summary(g)$dev.expl, 3),
    gam_AIC         = round(AIC(g), 1)
  )
}


# --- Part 1. Monte Carlo --------------------------------------------
# uniform draws within the range of each parameter

# B (bites/min): lower and upper quartile (Minguito-Frutos et al., 2025, Fig. 2a)
B_low  <- c(sl = 19.14,  sr = 22.07,  ss = 13.867, sc = 13.086)
B_high <- c(sl = 36.718, sr = 47.852, ss = 30.469, sc = 20.703)

# F (%): 95% Wilson CI from the number of stomachs; S. luridus: range of the three studies
wilson <- function(F, n) {
  k <- round(F / 100 * n)
  100 * prop.test(k, n, correct = FALSE)$conf.int
}
F_low  <- c(sl = 10,   sr = wilson(7, 253)[1], ss = wilson(83.33, 18)[1], sc = wilson(19.44, 36)[1])
F_high <- c(sl = 71.6, sr = wilson(7, 253)[2], ss = wilson(83.33, 18)[2], sc = wilson(19.44, 36)[2])

# I (g/bite): range across the three canopy algae (Verges et al., 2014a, Fig. 5)
# S. rivulatus = S. luridus; S. cretense: lowest S. luridus to highest S. salpa
I_sl <- c(0.0048, 0.0098)
I_ss <- c(0.0145, 0.0215)
I_sc <- c(0.0048, 0.0215)

# D (feeding hours/day)
D_range <- c(6, 10)

# P. lividus (g ind-1 day-1): range of the three trials (Pusceddu et al., 2021)
pl_range <- c(0.34, 0.58)

# A. lixula / P. lividus: means +/- SE (Agnetta et al., 2015)
al_num <- c(2.16 - 0.46, 2.16 + 0.46)
al_den <- c(26.47 - 2.07, 26.47 + 2.07)

n_draws <- 1000
set.seed(123)

mc_results <- list()
mc_CHI     <- matrix(NA, nrow = n_sites_total, ncol = n_draws)

for (d in 1:n_draws) {

  par   <- fish_par
  par$B <- runif(4, B_low, B_high)
  par$F <- runif(4, F_low, F_high)
  I_luridus <- runif(1, I_sl[1], I_sl[2])
  par$I <- c(I_luridus, I_luridus, runif(1, I_ss[1], I_ss[2]), runif(1, I_sc[1], I_sc[2]))
  par$D <- runif(1, D_range[1], D_range[2])

  pl <- runif(1, pl_range[1], pl_range[2])
  al <- runif(1, al_num[1], al_num[2]) / runif(1, al_den[1], al_den[2])

  C <- consumption_rates(par, pl, al)
  mc_CHI[, d] <- calc_chi(sites, C)

  mc_results[[d]] <- cbind(draw = d, run_index(mc_CHI[, d]))
}

mc_draws <- do.call(rbind, mc_results)
write.csv(mc_draws, "results/mc_draws.csv", row.names = FALSE)

# summary of the draws
mc_summary <- data.frame(
  n_draws                 = n_draws,
  rank_corr_median        = median(mc_draws$rank_corr),
  rank_corr_2.5           = quantile(mc_draws$rank_corr, 0.025),
  avg_p_below_0.05        = mean(mc_draws$avg_p < 0.05),
  max_p_below_0.05        = mean(mc_draws$max_p < 0.05),
  gam_p_below_0.05        = mean(mc_draws$gam_p < 0.05),
  avg_share_below_median  = median(mc_draws$avg_share_below),
  avg_share_below_2.5     = quantile(mc_draws$avg_share_below, 0.025),
  avg_share_below_97.5    = quantile(mc_draws$avg_share_below, 0.975),
  max_share_below_median  = median(mc_draws$max_share_below),
  max_share_below_2.5     = quantile(mc_draws$max_share_below, 0.025),
  max_share_below_97.5    = quantile(mc_draws$max_share_below, 0.975),
  row.names = NULL
)
print(t(mc_summary))
write.csv(mc_summary, "results/mc_summary.csv", row.names = FALSE)

# 95% range of the CHI per site
mc_site_CHI <- data.frame(
  id        = sites$id,
  site      = sites$site,
  CHI       = round(sites$CHI, 2),
  CHI_2.5   = round(apply(mc_CHI, 1, quantile, 0.025), 2),
  CHI_97.5  = round(apply(mc_CHI, 1, quantile, 0.975), 2)
)
write.csv(mc_site_CHI, "results/mc_site_CHI.csv", row.names = FALSE)

# figure: site CHI with its 95% range
fig_data <- mc_site_CHI[order(mc_site_CHI$CHI), ]
fig_data$rank <- 1:nrow(fig_data)

# sites below the baseline break (average cover)
base_break <- breakpoint_summary(make_bins(sites), "canopy")$sites_below

fig_mc <- ggplot(fig_data, aes(x = rank, y = CHI)) +
  geom_linerange(aes(ymin = CHI_2.5, ymax = CHI_97.5), colour = "steelblue", alpha = 0.5, linewidth = 1) +
  geom_point(colour = "steelblue", size = 1.5) +
  geom_vline(xintercept = base_break + 0.5, linetype = "dashed", colour = "grey40") +
  annotate("text", x = base_break, y = max(fig_data$CHI_97.5),
           label = paste0("baseline break\n(", base_break, " sites below)"),
           hjust = 1.05, vjust = 1, size = 3, colour = "grey40") +
  labs(x = "Sites ordered by baseline CHI", y = expression(CHI~(g~day^-1~m^-2))) +
  theme_classic()
fig_mc

# Supplementary Figure 3
dir.create("figures/supplementary figures", showWarnings = FALSE)
ggsave("figures/supplementary figures/FigS3_MC_site_CHI.png", fig_mc, width = 16, height = 8, units = "cm", dpi = 300)


# --- Part 2. Alternative indices ------------------------------------

C_base  <- consumption_rates(fish_par, pl_rate, al_ratio)
species <- names(C_base)

alternatives <- list()

alternatives[["CHI (baseline)"]] <- sites$CHI

# equal weighting: density / max density, summed over species
alternatives[["Equal weighting"]] <- rowSums(sapply(species, function(s) sites[[s]] / max(sites[[s]])))

# untransformed F (F/100 instead of F*)
C_raw <- C_base
C_raw[fish_par$species] <- fish_par$B * 60 * fish_par$D * (fish_par$F / 100) * fish_par$I
alternatives[["Untransformed F"]] <- calc_chi(sites, C_raw)

# leave one species out
for (s in species) {
  C_loo <- C_base
  C_loo[s] <- 0
  alternatives[[paste("Without", s)]] <- calc_chi(sites, C_loo)
}

alternative_indices <- do.call(rbind, lapply(names(alternatives), function(a)
  cbind(index = a, run_index(alternatives[[a]]))))

# mean share of each species in the CHI
share <- colMeans(sweep(as.matrix(sites[, species]), 2, C_base, "*")) / mean(sites$CHI)
alternative_indices$share_of_CHI <- NA
for (s in species) alternative_indices$share_of_CHI[alternative_indices$index == paste("Without", s)] <- round(share[s], 3)

print(alternative_indices)
write.csv(alternative_indices, "results/alternative_indices.csv", row.names = FALSE)


# --- Part 3. Binning (baseline CHI) ---------------------------------

binning_results <- list()

# fixed width: widths and origins
for (w in c(1.5, 2, 2.5, 3)) {
  for (o in seq(0, w - 0.5, by = 0.5)) {
    b <- make_bins(sites, width = w, origin = o)
    for (v in c("canopy", "canopy_max")) {
      binning_results[[length(binning_results) + 1]] <-
        cbind(test = "fixed width", width = w, origin = o, breakpoint_summary(b, v))
    }
  }
}

# no merging (width 2, origin 0)
b <- make_bins(sites, width = 2, origin = 0, min_sites = 1)
for (v in c("canopy", "canopy_max")) {
  binning_results[[length(binning_results) + 1]] <-
    cbind(test = "no merging", width = 2, origin = 0, breakpoint_summary(b, v))
}

# equal-count bins (k sites per bin; a smaller last bin is merged with the previous one)
for (k in 3:6) {
  s <- sites %>% arrange(CHI)
  s$bin <- ceiling(seq_len(nrow(s)) / k)
  if (sum(s$bin == max(s$bin)) < k) s$bin[s$bin == max(s$bin)] <- max(s$bin) - 1

  b <- s %>%
    group_by(bin) %>%
    summarise(CH = min(CHI), CHI_min = min(CHI), CHI_max = max(CHI), n_sites = n(),
              canopy = mean(canopy_avg), canopy_max = max(canopy_avg))

  for (v in c("canopy", "canopy_max")) {
    binning_results[[length(binning_results) + 1]] <-
      cbind(test = paste("equal count, k =", k), width = NA, origin = NA, breakpoint_summary(b, v))
  }
}

# no binning (sites ordered by CHI)
s <- sites %>% arrange(CHI)
s_bins <- data.frame(CH = s$CHI, CHI_max = s$CHI, n_sites = 1, canopy = s$canopy_avg)
binning_results[[length(binning_results) + 1]] <-
  cbind(test = "no binning (sites)", width = NA, origin = NA, breakpoint_summary(s_bins, "canopy"))

sensitivity_binning <- do.call(rbind, binning_results)
print(sensitivity_binning)

write.csv(sensitivity_binning, "results/sensitivity_binning.csv", row.names = FALSE)
