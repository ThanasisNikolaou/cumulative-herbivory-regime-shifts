# Cumulative herbivory pressure and regime shifts in eastern Mediterranean macroalgal forests

Athanasios Nikolaou¹, Deevesh Ashley Hemraj², Jacob Carstensen², George Tsirtsis¹, Stelios Katsanevakis¹

¹ Department of Marine Sciences, University of the Aegean, University Hill 81100, Mytilene, Lesvos, Greece

² Department of Ecoscience, Aarhus University, Frederiksborgvej 399, Roskilde, DK-4000 Denmark

## About

This repository contains the data and R code used in the manuscript. Running the scripts reproduces the analyses and figures.

## Scripts

| Script | Description |
|---|---|
| `00_CHI.R` | Consumption rates (Table 1) and CHI of each site |
| `01_bins.R` | Groups sites into CHI bins of 2 g consumed/day/m² |
| `02_breakpoint.R` | Breakpoint analysis of canopy cover along the CHI (Figure 5) |
| `03_gam.R` | Site-level GAM of canopy cover against the CHI (Figure 5, Supplementary Figure 4) |
| `04_sensitivity.R` | Sensitivity analysis of the CHI and the breakpoint (Supplementary Tables 1-3, Supplementary Figure 3) |
| `05_environment.R` | SST, chlorophyll-a and wave height in the two regions and in the GAM (Supplementary Tables 4-5) |
| `06_barplots.R` | Canopy cover, herbivores and CHI per site (Figure 2) |
| `07a_nmds.R` | NMDS of the herbivore community (Figure 3) |
| `07b_nmds_temperature.R` | NMDS1 against temperature and changepoint analysis (Figure 3) |
| `08a_fetch_occurrences.R` | Occurrence records of the six herbivores from OBIS and GBIF |
| `08b_godas_sst.R` | GODAS SST at each occurrence record |
| `08c_CTI.R` | Thermal preference of each species (Ti) and Community Temperature Index (CTI) (Supplementary Figures 1-2) |
| `08d_raincloud.R` | CTI and CHI in the two regions (Figure 4) |

The inputs and outputs of each script are listed at the top of the script.

## Data

`data/dataset_eastMed_reefs_abundance_density.csv` has one row per site (60 sites, 30 in NAS and 30 in SAL).

| Column | Unit | Description |
|---|---|---|
| `id`, `site` | | Site number and name |
| `region` | | NAS (North Aegean Sea) or SAL (South Aegean and Levantine Sea) |
| `lon`, `lat` | decimal degrees | Coordinates |
| `canopy_avg` | % | Cover of canopy-forming macroalgae, average of 0 and 5 m |
| `pl`, `al` | ind/m² | Density of *Paracentrotus lividus* and *Arbacia lixula* |
| `ss`, `sl`, `sr`, `sc` | ind/m² | Density of *Sarpa salpa*, *Siganus luridus*, *Siganus rivulatus* and *Sparisoma cretense* |
| `avg_temp`, `min_temp`, `max_temp` | °C | Average, minimum and maximum annual SST, 2016-2021 (Copernicus Marine Service, SST_MED_SST_L4_REP_OBSERVATIONS_010_021) |
| `CTI_med` | °C | Community Temperature Index (from `08c_CTI.R`) |

`data/dataset_eastMed_reefs_BARPLOTS.csv` is the same, with fish biomass (g/m²) instead of density. It is used for Figure 2.

`data/site_environment.csv` has the average surface chlorophyll-a (`chl_mean`, mg/m³, MEDSEA_MULTIYEAR_BGC_006_008) and significant wave height (`Hs_mean`, m, MEDSEA_MULTIYEAR_WAV_006_012) of 2016-2021, from the model cell nearest to each site.

## R packages

```r
install.packages(c(
  "ggplot2", "dplyr", "patchwork", "strucchange", "mgcv", "DHARMa", "gratia",
  "vegan", "changepoint", "tidyverse", "ggdist", "cowplot", "ggtext",
  "sf", "rnaturalearth", "rnaturalearthdata", "tidync", "robis", "rgbif"
))
```

## Contact

For questions about the data or code, please contact the corresponding author.
