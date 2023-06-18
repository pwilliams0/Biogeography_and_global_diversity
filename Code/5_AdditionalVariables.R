#============================================================#
# Compare R2 values from models using additional environment #
#                 variables besides climate                  #
#                  (Extended Data Table 2)                   #
#============================================================#

library(tidyverse)
sf_use_s2(FALSE) # Avoid some issues of invalid geometries

# ----- LOAD CELLS WITH ALL COVARIATES  ----------

cells <- st_read("Data/global_cells.shp") %>%
  # Landmass area
  left_join(read.csv("Data/Raw/landmass_area_cells.csv"), "cell_id") %>%
  # Elevation (mean and range)
  left_join(read.csv("Data/elev_cells.csv"), "cell_id") %>%
  dplyr::select(-X) %>%
  # Present climate
  left_join(read.csv("Data/clim_Present_cells.csv"), "cell_id") %>%
  dplyr::select(-X) %>%
  # Climate distance compared to present
  left_join(read.csv("Data/clim_diff_cells.csv"), "cell_id") %>%
  dplyr::select(-X) %>%
  # Past ice cover
  left_join(read.csv("Data/ice_cells.csv"), "cell_id") %>%
  dplyr::select(-X) %>%
  # Years since significant land conversion
  left_join(read.csv("Data/KK10_cells.csv"), "cell_id") %>%
  dplyr::select(-X) %>%
  # Human Impact Index
  left_join(read.csv("Data/HII_cells.csv"), "cell_id") %>%
  dplyr::select(-X)

# ----- BIRDS  ----------

# ---------- Species richness ----------

# Load data
data <- read.csv("Data/bird_SR_cells.csv") %>%
  left_join(cells, by="cell_id")

# Run model for climate-only, then each variable plus climate
# Calculate additional variance explained
SR_clim <- lm(SR ~ clim_pca_1 + I(clim_pca_1^2) +
                clim_pca_2 + I(clim_pca_2^2) +
                clim_pca_3 + I(clim_pca_3^2) +
                clim_pca_4 + I(clim_pca_4^2),
              data=data,
              na.action = na.fail)

# Landmass area
SR_Landmass_area <- lm(SR ~ log(Landmass_area) +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(SR_Landmass_area)$r.squared - summary(SR_clim)$r.squared

# Mean elevation
SR_elev_mean <- lm(SR ~ elev_mean + I(elev_mean^2) +
                     clim_pca_1 + I(clim_pca_1^2) +
                     clim_pca_2 + I(clim_pca_2^2) +
                     clim_pca_3 + I(clim_pca_3^2) +
                     clim_pca_4 + I(clim_pca_4^2),
                   data=data,
                   na.action = na.fail)
summary(SR_elev_mean)$r.squared - summary(SR_clim)$r.squared

# Elevation range
SR_elev_range <- lm(SR ~ elev_range +
                      clim_pca_1 + I(clim_pca_1^2) +
                      clim_pca_2 + I(clim_pca_2^2) +
                      clim_pca_3 + I(clim_pca_3^2) +
                      clim_pca_4 + I(clim_pca_4^2),
                    data=data,
                    na.action = na.fail)
summary(SR_elev_range)$r.squared - summary(SR_clim)$r.squared

# Climate distance, present vs. LGM
SR_clim_diff_LGM <- lm(SR ~ clim_diff_LGM +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(SR_clim_diff_LGM)$r.squared - summary(SR_clim)$r.squared

# Climate distance, present vs. MidHolo
SR_clim_diff_MidHol <- lm(SR ~ clim_diff_MidHol +
                            clim_pca_1 + I(clim_pca_1^2) +
                            clim_pca_2 + I(clim_pca_2^2) +
                            clim_pca_3 + I(clim_pca_3^2) +
                            clim_pca_4 + I(clim_pca_4^2),
                          data=data,
                          na.action = na.fail)
summary(SR_clim_diff_MidHol)$r.squared - summary(SR_clim)$r.squared

# Past ice cover
SR_ice_cover <- lm(SR ~ ice_cover +
                     clim_pca_1 + I(clim_pca_1^2) +
                     clim_pca_2 + I(clim_pca_2^2) +
                     clim_pca_3 + I(clim_pca_3^2) +
                     clim_pca_4 + I(clim_pca_4^2),
                   data=data,
                   na.action = na.fail)
summary(SR_ice_cover)$r.squared - summary(SR_clim)$r.squared

# Years since significant land conversion
SR_LandConvYears <- lm(SR ~ LandConvYears +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(SR_LandConvYears)$r.squared - summary(SR_clim)$r.squared

# Human Impact Index
SR_HII <- lm(SR ~ HII +
               clim_pca_1 + I(clim_pca_1^2) +
               clim_pca_2 + I(clim_pca_2^2) +
               clim_pca_3 + I(clim_pca_3^2) +
               clim_pca_4 + I(clim_pca_4^2),
             data=data,
             na.action = na.fail)
summary(SR_HII)$r.squared - summary(SR_clim)$r.squared

# ---------- Functional richness ----------

# Load data
data <- read.csv("Data/bird_FR_cells.csv") %>%
  left_join(cells, by="cell_id")

# Run model for climate-only, then each variable plus climate
# Calculate additional variance explained
FR_clim <- lm(FR ~ clim_pca_1 + I(clim_pca_1^2) +
                clim_pca_2 + I(clim_pca_2^2) +
                clim_pca_3 + I(clim_pca_3^2) +
                clim_pca_4 + I(clim_pca_4^2),
              data=data,
              na.action = na.fail)

# Landmass area
FR_Landmass_area <- lm(FR ~ log(Landmass_area) +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(FR_Landmass_area)$r.squared - summary(FR_clim)$r.squared

# Mean elevation
FR_elev_mean <- lm(FR ~ elev_mean + I(elev_mean^2) +
                     clim_pca_1 + I(clim_pca_1^2) +
                     clim_pca_2 + I(clim_pca_2^2) +
                     clim_pca_3 + I(clim_pca_3^2) +
                     clim_pca_4 + I(clim_pca_4^2),
                   data=data,
                   na.action = na.fail)
summary(FR_elev_mean)$r.squared - summary(FR_clim)$r.squared

# Elevation range
FR_elev_range <- lm(FR ~ elev_range +
                      clim_pca_1 + I(clim_pca_1^2) +
                      clim_pca_2 + I(clim_pca_2^2) +
                      clim_pca_3 + I(clim_pca_3^2) +
                      clim_pca_4 + I(clim_pca_4^2),
                    data=data,
                    na.action = na.fail)
summary(FR_elev_range)$r.squared - summary(FR_clim)$r.squared

# Climate distance, present vs. LGM
FR_clim_diff_LGM <- lm(FR ~ clim_diff_LGM +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(FR_clim_diff_LGM)$r.squared - summary(FR_clim)$r.squared

# Climate distance, present vs. MidHolo
FR_clim_diff_MidHol <- lm(FR ~ clim_diff_MidHol +
                            clim_pca_1 + I(clim_pca_1^2) +
                            clim_pca_2 + I(clim_pca_2^2) +
                            clim_pca_3 + I(clim_pca_3^2) +
                            clim_pca_4 + I(clim_pca_4^2),
                          data=data,
                          na.action = na.fail)
summary(FR_clim_diff_MidHol)$r.squared - summary(FR_clim)$r.squared

# Past ice cover
FR_ice_cover <- lm(FR ~ ice_cover +
                     clim_pca_1 + I(clim_pca_1^2) +
                     clim_pca_2 + I(clim_pca_2^2) +
                     clim_pca_3 + I(clim_pca_3^2) +
                     clim_pca_4 + I(clim_pca_4^2),
                   data=data,
                   na.action = na.fail)
summary(FR_ice_cover)$r.squared - summary(FR_clim)$r.squared

# Years since significant land conversion
FR_LandConvYears <- lm(FR ~ LandConvYears +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(FR_LandConvYears)$r.squared - summary(FR_clim)$r.squared

# Human Impact Index
FR_HII <- lm(FR ~ HII +
               clim_pca_1 + I(clim_pca_1^2) +
               clim_pca_2 + I(clim_pca_2^2) +
               clim_pca_3 + I(clim_pca_3^2) +
               clim_pca_4 + I(clim_pca_4^2),
             data=data,
             na.action = na.fail)
summary(FR_HII)$r.squared - summary(FR_clim)$r.squared

# ----- MAMMALS  ----------

# ---------- Species richness ----------

# Load data
data <- read.csv("Data/mamm_SR_cells.csv") %>%
  left_join(cells, by="cell_id")

# Run model for climate-only, then each variable plus climate
# Calculate additional variance explained
SR_clim <- lm(SR ~ clim_pca_1 + I(clim_pca_1^2) +
                clim_pca_2 + I(clim_pca_2^2) +
                clim_pca_3 + I(clim_pca_3^2) +
                clim_pca_4 + I(clim_pca_4^2),
              data=data,
              na.action = na.fail)

# Landmass area
SR_Landmass_area <- lm(SR ~ log(Landmass_area) +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(SR_Landmass_area)$r.squared - summary(SR_clim)$r.squared

# Mean elevation
SR_elev_mean <- lm(SR ~ elev_mean + I(elev_mean^2) +
                     clim_pca_1 + I(clim_pca_1^2) +
                     clim_pca_2 + I(clim_pca_2^2) +
                     clim_pca_3 + I(clim_pca_3^2) +
                     clim_pca_4 + I(clim_pca_4^2),
                   data=data,
                   na.action = na.fail)
summary(SR_elev_mean)$r.squared - summary(SR_clim)$r.squared

# Elevation range
SR_elev_range <- lm(SR ~ elev_range +
                      clim_pca_1 + I(clim_pca_1^2) +
                      clim_pca_2 + I(clim_pca_2^2) +
                      clim_pca_3 + I(clim_pca_3^2) +
                      clim_pca_4 + I(clim_pca_4^2),
                    data=data,
                    na.action = na.fail)
summary(SR_elev_range)$r.squared - summary(SR_clim)$r.squared

# Climate distance, present vs. LGM
SR_clim_diff_LGM <- lm(SR ~ clim_diff_LGM +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(SR_clim_diff_LGM)$r.squared - summary(SR_clim)$r.squared

# Climate distance, present vs. MidHolo
SR_clim_diff_MidHol <- lm(SR ~ clim_diff_MidHol +
                            clim_pca_1 + I(clim_pca_1^2) +
                            clim_pca_2 + I(clim_pca_2^2) +
                            clim_pca_3 + I(clim_pca_3^2) +
                            clim_pca_4 + I(clim_pca_4^2),
                          data=data,
                          na.action = na.fail)
summary(SR_clim_diff_MidHol)$r.squared - summary(SR_clim)$r.squared

# Past ice cover
SR_ice_cover <- lm(SR ~ ice_cover +
                     clim_pca_1 + I(clim_pca_1^2) +
                     clim_pca_2 + I(clim_pca_2^2) +
                     clim_pca_3 + I(clim_pca_3^2) +
                     clim_pca_4 + I(clim_pca_4^2),
                   data=data,
                   na.action = na.fail)
summary(SR_ice_cover)$r.squared - summary(SR_clim)$r.squared

# Years since significant land conversion
SR_LandConvYears <- lm(SR ~ LandConvYears +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(SR_LandConvYears)$r.squared - summary(SR_clim)$r.squared

# Human Impact Index
SR_HII <- lm(SR ~ HII +
               clim_pca_1 + I(clim_pca_1^2) +
               clim_pca_2 + I(clim_pca_2^2) +
               clim_pca_3 + I(clim_pca_3^2) +
               clim_pca_4 + I(clim_pca_4^2),
             data=data,
             na.action = na.fail)
summary(SR_HII)$r.squared - summary(SR_clim)$r.squared

# ---------- Functional richness ----------

# Load data
data <- read.csv("Data/mamm_FR_cells.csv") %>%
  left_join(cells, by="cell_id")

# Run model for climate-only, then each variable plus climate
# Calculate additional variance explained
FR_clim <- lm(FR ~ clim_pca_1 + I(clim_pca_1^2) +
                clim_pca_2 + I(clim_pca_2^2) +
                clim_pca_3 + I(clim_pca_3^2) +
                clim_pca_4 + I(clim_pca_4^2),
              data=data,
              na.action = na.fail)

# Landmass area
FR_Landmass_area <- lm(FR ~ log(Landmass_area) +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(FR_Landmass_area)$r.squared - summary(FR_clim)$r.squared

# Mean elevation
FR_elev_mean <- lm(FR ~ elev_mean + I(elev_mean^2) +
                     clim_pca_1 + I(clim_pca_1^2) +
                     clim_pca_2 + I(clim_pca_2^2) +
                     clim_pca_3 + I(clim_pca_3^2) +
                     clim_pca_4 + I(clim_pca_4^2),
                   data=data,
                   na.action = na.fail)
summary(FR_elev_mean)$r.squared - summary(FR_clim)$r.squared

# Elevation range
FR_elev_range <- lm(FR ~ elev_range +
                      clim_pca_1 + I(clim_pca_1^2) +
                      clim_pca_2 + I(clim_pca_2^2) +
                      clim_pca_3 + I(clim_pca_3^2) +
                      clim_pca_4 + I(clim_pca_4^2),
                    data=data,
                    na.action = na.fail)
summary(FR_elev_range)$r.squared - summary(FR_clim)$r.squared

# Climate distance, present vs. LGM
FR_clim_diff_LGM <- lm(FR ~ clim_diff_LGM +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(FR_clim_diff_LGM)$r.squared - summary(FR_clim)$r.squared

# Climate distance, present vs. MidHolo
FR_clim_diff_MidHol <- lm(FR ~ clim_diff_MidHol +
                            clim_pca_1 + I(clim_pca_1^2) +
                            clim_pca_2 + I(clim_pca_2^2) +
                            clim_pca_3 + I(clim_pca_3^2) +
                            clim_pca_4 + I(clim_pca_4^2),
                          data=data,
                          na.action = na.fail)
summary(FR_clim_diff_MidHol)$r.squared - summary(FR_clim)$r.squared

# Past ice cover
FR_ice_cover <- lm(FR ~ ice_cover +
                     clim_pca_1 + I(clim_pca_1^2) +
                     clim_pca_2 + I(clim_pca_2^2) +
                     clim_pca_3 + I(clim_pca_3^2) +
                     clim_pca_4 + I(clim_pca_4^2),
                   data=data,
                   na.action = na.fail)
summary(FR_ice_cover)$r.squared - summary(FR_clim)$r.squared

# Years since significant land conversion
FR_LandConvYears <- lm(FR ~ LandConvYears +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(FR_LandConvYears)$r.squared - summary(FR_clim)$r.squared

# Human Impact Index
FR_HII <- lm(FR ~ HII +
               clim_pca_1 + I(clim_pca_1^2) +
               clim_pca_2 + I(clim_pca_2^2) +
               clim_pca_3 + I(clim_pca_3^2) +
               clim_pca_4 + I(clim_pca_4^2),
             data=data,
             na.action = na.fail)
summary(FR_HII)$r.squared - summary(FR_clim)$r.squared

# ----- BATS  ----------

# ---------- Species richness ----------

# Load data
data <- read.csv("Data/bat_SR_cells.csv") %>%
  left_join(cells, by="cell_id")

# Run model for climate-only, then each variable plus climate
# Calculate additional variance explained
SR_clim <- lm(SR ~ clim_pca_1 + I(clim_pca_1^2) +
                clim_pca_2 + I(clim_pca_2^2) +
                clim_pca_3 + I(clim_pca_3^2) +
                clim_pca_4 + I(clim_pca_4^2),
              data=data,
              na.action = na.fail)

# Landmass area
SR_Landmass_area <- lm(SR ~ log(Landmass_area) +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(SR_Landmass_area)$r.squared - summary(SR_clim)$r.squared

# Mean elevation
SR_elev_mean <- lm(SR ~ elev_mean + I(elev_mean^2) +
                     clim_pca_1 + I(clim_pca_1^2) +
                     clim_pca_2 + I(clim_pca_2^2) +
                     clim_pca_3 + I(clim_pca_3^2) +
                     clim_pca_4 + I(clim_pca_4^2),
                   data=data,
                   na.action = na.fail)
summary(SR_elev_mean)$r.squared - summary(SR_clim)$r.squared

# Elevation range
SR_elev_range <- lm(SR ~ elev_range +
                      clim_pca_1 + I(clim_pca_1^2) +
                      clim_pca_2 + I(clim_pca_2^2) +
                      clim_pca_3 + I(clim_pca_3^2) +
                      clim_pca_4 + I(clim_pca_4^2),
                    data=data,
                    na.action = na.fail)
summary(SR_elev_range)$r.squared - summary(SR_clim)$r.squared

# Climate distance, present vs. LGM
SR_clim_diff_LGM <- lm(SR ~ clim_diff_LGM +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(SR_clim_diff_LGM)$r.squared - summary(SR_clim)$r.squared

# Climate distance, present vs. MidHolo
SR_clim_diff_MidHol <- lm(SR ~ clim_diff_MidHol +
                            clim_pca_1 + I(clim_pca_1^2) +
                            clim_pca_2 + I(clim_pca_2^2) +
                            clim_pca_3 + I(clim_pca_3^2) +
                            clim_pca_4 + I(clim_pca_4^2),
                          data=data,
                          na.action = na.fail)
summary(SR_clim_diff_MidHol)$r.squared - summary(SR_clim)$r.squared

# Past ice cover
SR_ice_cover <- lm(SR ~ ice_cover +
                     clim_pca_1 + I(clim_pca_1^2) +
                     clim_pca_2 + I(clim_pca_2^2) +
                     clim_pca_3 + I(clim_pca_3^2) +
                     clim_pca_4 + I(clim_pca_4^2),
                   data=data,
                   na.action = na.fail)
summary(SR_ice_cover)$r.squared - summary(SR_clim)$r.squared

# Years since significant land conversion
SR_LandConvYears <- lm(SR ~ LandConvYears +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(SR_LandConvYears)$r.squared - summary(SR_clim)$r.squared

# Human Impact Index
SR_HII <- lm(SR ~ HII +
               clim_pca_1 + I(clim_pca_1^2) +
               clim_pca_2 + I(clim_pca_2^2) +
               clim_pca_3 + I(clim_pca_3^2) +
               clim_pca_4 + I(clim_pca_4^2),
             data=data,
             na.action = na.fail)
summary(SR_HII)$r.squared - summary(SR_clim)$r.squared

# ---------- Functional richness ----------

# Load data
data <- read.csv("Data/bat_FR_cells.csv") %>%
  na.omit() %>%
  left_join(cells, by="cell_id")

# Run model for climate-only, then each variable plus climate
# Calculate additional variance explained
FR_clim <- lm(FR ~ clim_pca_1 + I(clim_pca_1^2) +
                clim_pca_2 + I(clim_pca_2^2) +
                clim_pca_3 + I(clim_pca_3^2) +
                clim_pca_4 + I(clim_pca_4^2),
              data=data,
              na.action = na.fail)

# Landmass area
FR_Landmass_area <- lm(FR ~ log(Landmass_area) +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(FR_Landmass_area)$r.squared - summary(FR_clim)$r.squared

# Mean elevation
FR_elev_mean <- lm(FR ~ elev_mean + I(elev_mean^2) +
                     clim_pca_1 + I(clim_pca_1^2) +
                     clim_pca_2 + I(clim_pca_2^2) +
                     clim_pca_3 + I(clim_pca_3^2) +
                     clim_pca_4 + I(clim_pca_4^2),
                   data=data,
                   na.action = na.fail)
summary(FR_elev_mean)$r.squared - summary(FR_clim)$r.squared

# Elevation range
FR_elev_range <- lm(FR ~ elev_range +
                      clim_pca_1 + I(clim_pca_1^2) +
                      clim_pca_2 + I(clim_pca_2^2) +
                      clim_pca_3 + I(clim_pca_3^2) +
                      clim_pca_4 + I(clim_pca_4^2),
                    data=data,
                    na.action = na.fail)
summary(FR_elev_range)$r.squared - summary(FR_clim)$r.squared

# Climate distance, present vs. LGM
FR_clim_diff_LGM <- lm(FR ~ clim_diff_LGM +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(FR_clim_diff_LGM)$r.squared - summary(FR_clim)$r.squared

# Climate distance, present vs. MidHolo
FR_clim_diff_MidHol <- lm(FR ~ clim_diff_MidHol +
                            clim_pca_1 + I(clim_pca_1^2) +
                            clim_pca_2 + I(clim_pca_2^2) +
                            clim_pca_3 + I(clim_pca_3^2) +
                            clim_pca_4 + I(clim_pca_4^2),
                          data=data,
                          na.action = na.fail)
summary(FR_clim_diff_MidHol)$r.squared - summary(FR_clim)$r.squared

# Past ice cover
FR_ice_cover <- lm(FR ~ ice_cover +
                     clim_pca_1 + I(clim_pca_1^2) +
                     clim_pca_2 + I(clim_pca_2^2) +
                     clim_pca_3 + I(clim_pca_3^2) +
                     clim_pca_4 + I(clim_pca_4^2),
                   data=data,
                   na.action = na.fail)
summary(FR_ice_cover)$r.squared - summary(FR_clim)$r.squared

# Years since significant land conversion
FR_LandConvYears <- lm(FR ~ LandConvYears +
                         clim_pca_1 + I(clim_pca_1^2) +
                         clim_pca_2 + I(clim_pca_2^2) +
                         clim_pca_3 + I(clim_pca_3^2) +
                         clim_pca_4 + I(clim_pca_4^2),
                       data=data,
                       na.action = na.fail)
summary(FR_LandConvYears)$r.squared - summary(FR_clim)$r.squared

# Human Impact Index
FR_HII <- lm(FR ~ HII +
               clim_pca_1 + I(clim_pca_1^2) +
               clim_pca_2 + I(clim_pca_2^2) +
               clim_pca_3 + I(clim_pca_3^2) +
               clim_pca_4 + I(clim_pca_4^2),
             data=data,
             na.action = na.fail)
summary(FR_HII)$r.squared - summary(FR_clim)$r.squared
