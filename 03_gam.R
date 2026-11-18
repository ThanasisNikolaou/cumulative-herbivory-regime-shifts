# Figure 5: site-level GAM of canopy cover against the CHI (n = 60)
#
# canopy_avg ~ s(CHI, k = 4), Tweedie (log link)
#
# Input:  data/sites_CHI.csv
# Output: results/gam_summary.txt, results/gam_region_comparison.csv,
#         results/gam_significant_change.csv, results/gam_curve.csv,
#         figures/GAM_*.png, figures/fig5_panel_gam.png (Figure 5c)
#         figures/supplementary figures/FigS4_DHARMa_diagnostics.png

library(mgcv)
library(DHARMa)
library(gratia)
library(ggplot2)
library(dplyr)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

sites <- read.csv("data/sites_CHI.csv")
sites$region <- factor(sites$region)


# --- Main model ----------------------------------------------------

gam_model <- gam(canopy_avg ~ s(CHI, k = 4),
                 family = tw(link = "log"),
                 data = sites, method = "REML")

summary(gam_model)

capture.output(summary(gam_model), file = "results/gam_summary.txt")

# residual diagnostics
set.seed(123)
sim_res <- simulateResiduals(fittedModel = gam_model, n = 1000)

png("figures/GAM_DHARMa_diagnostics.png", width = 2400, height = 1200, res = 300)
plot(sim_res)
dev.off()

# Supplementary Figure 4
dir.create("figures/supplementary figures", showWarnings = FALSE)
png("figures/supplementary figures/FigS4_DHARMa_diagnostics.png", width = 2400, height = 1200, res = 300)
plot(sim_res)
dev.off()

testUniformity(sim_res)
testDispersion(sim_res)
testZeroInflation(sim_res)


# --- Region --------------------------------------------------------

gam_region <- gam(canopy_avg ~ s(CHI, k = 4) + region,
                  family = tw(link = "log"),
                  data = sites, method = "REML")

gam_interaction <- gam(canopy_avg ~ s(CHI, k = 4, by = region) + region,
                       family = tw(link = "log"),
                       data = sites, method = "REML")

summary(gam_region)
summary(gam_interaction)

region_comparison <- data.frame(
  model = c("CHI", "CHI + region", "CHI x region"),
  AIC   = round(AIC(gam_model, gam_region, gam_interaction)$AIC, 1),
  deviance_explained = round(c(summary(gam_model)$dev.expl,
                               summary(gam_region)$dev.expl,
                               summary(gam_interaction)$dev.expl), 3)
)
region_comparison

write.csv(region_comparison, "results/gam_region_comparison.csv", row.names = FALSE)


# --- First derivative (SiZer) ---------------------------------------
# simultaneous intervals; derivative on the link scale

set.seed(123)
deriv <- derivatives(gam_model,
                     select   = "s(CHI)",
                     type     = "central",
                     n        = 200,
                     interval = "simultaneous",
                     n_sim    = 10000)

draw(deriv, add_change = TRUE, change_type = "sizer")

# significant slope: interval does not include 0
deriv$significant <- (deriv$.lower_ci > 0) | (deriv$.upper_ci < 0)

# stretches of significant change
runs   <- rle(deriv$significant)
ends   <- cumsum(runs$lengths)
starts <- ends - runs$lengths + 1
significant_change <- data.frame(CHI_from = deriv$CHI[starts],
                                 CHI_to   = deriv$CHI[ends],
                                 significant = runs$values) %>%
  filter(significant)
significant_change

write.csv(significant_change, "results/gam_significant_change.csv", row.names = FALSE)


# --- Figure --------------------------------------------------------

newdata <- data.frame(CHI = seq(min(sites$CHI), max(sites$CHI), length.out = 200))
pred <- predict(gam_model, newdata, type = "link", se.fit = TRUE)

newdata$fit   <- exp(pred$fit)
newdata$upper <- exp(pred$fit + 1.96 * pred$se.fit)
newdata$lower <- exp(pred$fit - 1.96 * pred$se.fit)
write.csv(newdata, "results/gam_curve.csv", row.names = FALSE)   # fitted curve (Figure 5c)
p_gam <- ggplot() +
  geom_rect(data = significant_change,
            aes(xmin = CHI_from, xmax = CHI_to, ymin = -Inf, ymax = Inf),
            fill = "#de2d26", alpha = 0.2) +
  geom_ribbon(data = newdata, aes(x = CHI, ymin = lower, ymax = upper),
              fill = "grey70", alpha = 0.4) +
  geom_line(data = newdata, aes(x = CHI, y = fit), linewidth = 1.2) +
  geom_point(data = sites, aes(x = CHI, y = canopy_avg, colour = region),
             size = 2, alpha = 0.8) +
  scale_colour_manual(values = c(NAS = "#7aa3cd", SAL = "#d17480")) +
  coord_cartesian(ylim = c(0, 60)) +
  labs(x = expression("Cumulative Herbivory Index (g day"^-1*" m"^-2*")"),
       y = "Canopy cover (%)", colour = "") +
  theme_classic(base_size = 14) +
  theme(axis.text    = element_text(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

p_gam

ggsave("figures/GAM_canopy_vs_CHI.png", p_gam, width = 6, height = 5, dpi = 300)

# Figure 5 panel c: same axes as panels a, b (02_breakpoint.R), 140 x 68.4 mm
# panels a + b + c stacked = 140 x 190 mm
p_panel <- ggplot() +
  geom_rect(data = significant_change,
            aes(xmin = CHI_from, xmax = CHI_to, ymin = -Inf, ymax = Inf),
            fill = "#de2d26", alpha = 0.2) +
  geom_ribbon(data = newdata, aes(x = CHI, ymin = lower, ymax = upper),
              fill = "grey70", alpha = 0.4) +
  geom_line(data = newdata, aes(x = CHI, y = fit), linewidth = 0.7) +
  geom_point(data = sites, aes(x = CHI, y = canopy_avg, colour = region),
             size = 1.2, alpha = 0.85) +
  scale_colour_manual(values = c(NAS = "#7aa3cd", SAL = "#d17480")) +
  scale_x_continuous(breaks = seq(0, 40, 10)) +
  scale_y_continuous(breaks = seq(0, 60, 20)) +
  coord_cartesian(xlim = c(-1, 40), ylim = c(-3, 63), expand = FALSE) +
  labs(x = expression(bold("CHI (g consumed/day/m"^2*")")),
       y = "Canopy cover (%)") +
  theme_classic(base_size = 8) +
  theme(axis.text       = element_text(color = "black"),
        axis.title      = element_text(face = "bold", size = 9),
        legend.position = "none",
        panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.4),
        plot.margin     = margin(2, 3, 2, 2, "mm"))

ggsave("figures/fig5_panel_gam.png", p_panel, width = 140, height = 68.4, units = "mm", dpi = 600)
