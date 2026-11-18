# 05_environment.R
# Environmental differences between regions (SST, chl-a, wave height) and their relation to the CHI and canopy cover
#
# 1. regional differences (Shapiro-Wilk, Wilcoxon)
# 2. Spearman correlations with the CHI and among the variables
# 3. site-level GAMs of canopy cover with chl-a and Hs as covariates
#
# Input:  data/sites_CHI.csv, data/site_environment.csv
# Output: results/env_region_comparison.csv, env_correlations.csv,
#         env_gam_comparison.csv, env_gam_summary.txt

library(mgcv)

dir.create("results", showWarnings = FALSE)

sites <- read.csv("data/sites_CHI.csv")
env   <- read.csv("data/site_environment.csv")
sites <- merge(sites, env[, c("id", "chl_mean", "Hs_mean", "chl_dist_km", "wave_dist_km")], by = "id")
sites$SST    <- sites$avg_temp
sites$chl    <- sites$chl_mean
sites$Hs     <- sites$Hs_mean
sites$region <- factor(sites$region)

# sites whose nearest model cell is > 10 km away are left out (steps 1-3)
max_dist_km <- 10
far <- sites$chl_dist_km > max_dist_km | sites$wave_dist_km > max_dist_km
cat("Left out of the environmental analyses (model cell >", max_dist_km, "km away):",
    paste(sites$site[far], collapse = ", "), "\n")
env_sites <- sites[!far, ]

env_vars <- c("SST", "chl", "Hs")


# --- 1. Regional differences ---------------------------------------

region_comparison <- do.call(rbind, lapply(env_vars, function(v) {
  nas <- env_sites[[v]][env_sites$region == "NAS"]
  sal <- env_sites[[v]][env_sites$region == "SAL"]
  data.frame(variable     = v,
             NAS_median   = round(median(nas), 3),
             NAS_IQR      = paste(round(quantile(nas, c(0.25, 0.75)), 3), collapse = "-"),
             NAS_mean     = round(mean(nas), 3),
             NAS_SD       = round(sd(nas), 3),
             SAL_median   = round(median(sal), 3),
             SAL_IQR      = paste(round(quantile(sal, c(0.25, 0.75)), 3), collapse = "-"),
             SAL_mean     = round(mean(sal), 3),
             SAL_SD       = round(sd(sal), 3),
             shapiro_p_NAS = signif(shapiro.test(nas)$p.value, 3),
             shapiro_p_SAL = signif(shapiro.test(sal)$p.value, 3),
             W            = unname(wilcox.test(nas, sal, exact = FALSE)$statistic),
             p_value      = signif(wilcox.test(nas, sal, exact = FALSE)$p.value, 3))
}))
print(region_comparison)
write.csv(region_comparison, "results/env_region_comparison.csv", row.names = FALSE)


# --- 2. Correlations -----------------------------------------------

pairs <- list(c("CHI", "SST"), c("CHI", "chl"), c("CHI", "Hs"),
              c("SST", "chl"), c("SST", "Hs"), c("chl", "Hs"))
correlations <- do.call(rbind, lapply(pairs, function(p) {
  test <- cor.test(env_sites[[p[1]]], env_sites[[p[2]]], method = "spearman", exact = FALSE)
  data.frame(variable_1 = p[1], variable_2 = p[2], rho = round(unname(test$estimate), 2),
             p_value = signif(test$p.value, 3))
}))
print(correlations)
write.csv(correlations, "results/env_correlations.csv", row.names = FALSE)


# --- 3. GAMs with environmental covariates --------------------------

# chl on log10 scale (right-skewed)
env_sites$log_chl <- log10(env_sites$chl)

fit <- function(formula) gam(formula, family = tw(link = "log"), data = env_sites, method = "REML")

gam_chi    <- fit(canopy_avg ~ s(CHI, k = 4))
gam_env    <- fit(canopy_avg ~ s(log_chl, k = 3) + s(Hs, k = 3))
gam_full   <- fit(canopy_avg ~ s(CHI, k = 4) + s(log_chl, k = 3) + s(Hs, k = 3))

# region not included: collinear with SST, chl and Hs

# without the site with the highest chl
top_chl    <- env_sites$site[which.max(env_sites$chl)]
gam_no_top <- gam(canopy_avg ~ s(CHI, k = 4) + s(log_chl, k = 3) + s(Hs, k = 3),
                  family = tw(link = "log"), data = env_sites[env_sites$site != top_chl, ], method = "REML")

term_p <- function(model, term) {
  st <- summary(model)$s.table
  if (term %in% rownames(st)) signif(st[term, "p-value"], 3) else NA
}

models <- list(gam_chi, gam_env, gam_full, gam_no_top)
gam_comparison <- data.frame(
  model         = c("CHI", "chl + Hs", "CHI + chl + Hs",
                    paste0("CHI + chl + Hs, without ", top_chl)),
  n_sites       = sapply(models, function(m) nrow(m$model)),
  AIC           = round(sapply(models, AIC), 1),
  dev_explained = round(sapply(models, function(m) summary(m)$dev.expl), 3),
  p_CHI         = sapply(models, term_p, term = "s(CHI)"),
  p_chl         = sapply(models, term_p, term = "s(log_chl)"),
  p_Hs          = sapply(models, term_p, term = "s(Hs)")
)
print(gam_comparison)
write.csv(gam_comparison, "results/env_gam_comparison.csv", row.names = FALSE)

# concurvity
sink("results/env_gam_summary.txt")
print(summary(gam_full))
cat("\nConcurvity (full model)\n")
print(round(concurvity(gam_full, full = TRUE), 2))
sink()

# CHI effect with and without covariates
new_data <- data.frame(CHI = seq(0, max(env_sites$CHI), length.out = 100),
                       log_chl = median(env_sites$log_chl), Hs = median(env_sites$Hs))
curves <- data.frame(CHI         = new_data$CHI,
                     CHI_only    = predict(gam_chi,  new_data, type = "response"),
                     with_covars = predict(gam_full, new_data, type = "response"))
plot(curves$CHI, curves$CHI_only, type = "l", lwd = 2, ylim = c(0, max(curves[, -1]) * 1.1),
     xlab = "CHI (g day-1 m-2)", ylab = "Predicted canopy cover (%)")
lines(curves$CHI, curves$with_covars, lwd = 2, lty = 2)
legend("topright", c("CHI only", "CHI + chl + Hs (at their medians)"), lty = 1:2, lwd = 2, bty = "n")
