# Estimate Cumulative Herbivory Index (CHI) per site (g canopy day-1 m-2)
#
# CHI = sum of C_s * A_s (C_s: consumption rate, g ind-1 day-1; A_s: density, ind m-2)
# Fish:        C_s = B_s * 60 * D * F*_s * I_s (B in bites/min, D in hours), with F*_s = log10(F_s + 1) / log10(101)
# Sea urchins: C_s from feeding experiments
#
# Input:  data/dataset_eastMed_reefs_abundance_density.csv
# Output: data/sites_CHI.csv

library(dplyr)


# --- Parameters (Table 1) --------------------------------------------

# S. luridus F: mean of three studies weighted by number of stomachs
# (Azzurro et al., 2007; Stergiou, 1988; Shakman et al., 2009)
F_luridus <- weighted.mean(c(36.84, 71.6, 10), w = c(40, 209, 394))

fish_par <- data.frame(
  species = c("sl", "sr", "ss", "sc"),
  name    = c("Siganus luridus", "Siganus rivulatus", "Sarpa salpa", "Sparisoma cretense"),
  B = c(27.69, 36.08, 23.22, 17.48),        # bites/min (Minguito-Frutos et al., 2025)
  F = c(F_luridus, 7, 83.33, 19.44),        # % stomachs with canopy algae
  I = c(0.00687, 0.00687, 0.0182, 0.0182),  # g/bite (Verges et al., 2014a)
  D = 8                                     # feeding hours/day
)

pl_rate  <- 0.44            # P. lividus, g ind-1 day-1 (Pusceddu et al., 2021)
al_ratio <- 2.16 / 26.47    # A. lixula / P. lividus, gut contents (Agnetta et al., 2015)


# --- Functions -------------------------------------------------------

# consumption rate per species (g ind-1 day-1)
consumption_rates <- function(fish_par, pl_rate, al_ratio) {
  F_star <- log10(fish_par$F + 1) / log10(101)
  C_fish <- fish_par$B * 60 * fish_par$D * F_star * fish_par$I  # x 60: B is in bites/min and D in hours

  C <- c(C_fish, pl_rate, pl_rate * al_ratio)
  names(C) <- c(fish_par$species, "pl", "al")
  C
}

# CHI per site
calc_chi <- function(sites, C) {
  densities <- as.matrix(sites[, names(C)])
  as.vector(densities %*% C)
}


# --- CHI per site ----------------------------------------------------

sites <- read.csv("data/dataset_eastMed_reefs_abundance_density.csv")

C <- consumption_rates(fish_par, pl_rate, al_ratio)
round(C, 3)

sites$CHI <- calc_chi(sites, C)

# drop the original CHI column
sites <- sites %>% select(-CHi)
write.csv(sites, "data/sites_CHI.csv", row.names = FALSE)
