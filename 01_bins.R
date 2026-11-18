# Bins of fixed CHI width with mean and maximum canopy cover per bin (input of 02_breakpoint.R)
#
# Input:  data/sites_CHI.csv
# Output: data/bins_w2.csv

library(dplyr)

sites <- read.csv("data/sites_CHI.csv")


make_bins <- function(sites, width = 2, origin = 0, min_sites = 2) {

  sites$bin <- floor((sites$CHI - origin) / width)

  # merge sparse bins
  repeat {
    counts <- table(sites$bin)
    sparse <- as.numeric(names(counts)[counts < min_sites])
    if (length(sparse) == 0 || length(counts) <= 2) break

    top   <- max(sparse)
    lower <- as.numeric(names(counts))
    lower <- lower[lower < top]
    if (length(lower) == 0) break
    sites$bin[sites$bin == top] <- max(lower)
  }

  sites %>%
    group_by(bin) %>%
    summarise(
      CH         = max(0, origin + first(bin) * width),  # lower edge of the bin (not below 0)
      CHI_min    = min(CHI),
      CHI_max    = max(CHI),
      n_sites    = n(),
      n_NAS      = sum(region == "NAS"),
      n_SAL      = sum(region == "SAL"),
      canopy     = mean(canopy_avg),
      canopy_max = max(canopy_avg)
    ) %>%
    arrange(CH) %>%
    select(-bin)
}


bins <- make_bins(sites, width = 2)
print(bins, n = Inf)

write.csv(bins, "data/bins_w2.csv", row.names = FALSE)
