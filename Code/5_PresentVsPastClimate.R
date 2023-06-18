#===================================================================#
# Compare R2 values from models using climate data from present day #
#  (Present), last glacial maximum (LGM), & mid-Holocene (MidHol)   #
#                      (Extended Data Table 3)                      #
#===================================================================#

library(tidyverse)
sf_use_s2(FALSE) # Avoid some issues of invalid geometries

# ----- BIRDS  ----------

# Load cells with covariates
cells <- st_read("Data/global_cells.shp") %>%
  # Landmass area
  left_join(read.csv("Data/Raw/landmass_area_cells.csv"), "cell_id") %>%
  # Elevation (mean and range)
  left_join(read.csv("Data/elev_cells.csv"), "cell_id") %>%
  dplyr::select(-X) %>%
  # Phylobetadiversity NMDS axes
  left_join(read.csv("Data/bird_pb_NMDS.csv"), "cell_id") %>%
  dplyr::select(-X)

# Present climate
cells_Present <- cells %>%
  left_join(read.csv("Data/clim_Present_cells.csv"), "cell_id") %>%
  dplyr::select(-X)

# LGM climate
cells_LGM <- cells %>%
  left_join(read.csv("Data/clim_LGM_cells.csv"), "cell_id") %>%
  dplyr::select(-X)

# MidHol climate
cells_MidHol <- cells %>%
  left_join(read.csv("Data/clim_MidHol_cells.csv"), "cell_id") %>%
  dplyr::select(-X)

# ---------- Species richness ----------

# Load SR data for each climate
data_Present <- read.csv("Data/bird_SR_cells.csv") %>%
  left_join(cells_Present, by="cell_id")
data_LGM <- read.csv("Data/bird_SR_cells.csv") %>%
  left_join(cells_LGM, by="cell_id")
data_MidHol <- read.csv("Data/bird_SR_cells.csv") %>%
  left_join(cells_MidHol, by="cell_id")

# Run model for each climate, environment variables only
SR_env_Present <- lm(SR ~ log(Landmass_area) +
                       elev_mean + I(elev_mean^2) +
                       elev_range +
                       clim_pca_1 + I(clim_pca_1^2) +
                       clim_pca_2 + I(clim_pca_2^2) +
                       clim_pca_3 + I(clim_pca_3^2) +
                       clim_pca_4 + I(clim_pca_4^2),
                     data=data_Present,
                     na.action = na.fail)
SR_env_LGM <- lm(SR ~ log(Landmass_area) +
                   elev_mean + I(elev_mean^2) +
                   elev_range +
                   clim_pca_1 + I(clim_pca_1^2) +
                   clim_pca_2 + I(clim_pca_2^2) +
                   clim_pca_3 + I(clim_pca_3^2) +
                   clim_pca_4 + I(clim_pca_4^2),
                 data=data_LGM,
                 na.action = na.fail)
SR_env_MidHol <- lm(SR ~ log(Landmass_area) +
                      elev_mean + I(elev_mean^2) +
                      elev_range +
                      clim_pca_1 + I(clim_pca_1^2) +
                      clim_pca_2 + I(clim_pca_2^2) +
                      clim_pca_3 + I(clim_pca_3^2) +
                      clim_pca_4 + I(clim_pca_4^2),
                    data=data_MidHol,
                    na.action = na.fail)
summary(SR_env_Present)$r.squared
summary(SR_env_LGM)$r.squared
summary(SR_env_MidHol)$r.squared

# Run model for each climate, environment and phylobetadiversity
SR_all_Present <- lm(SR ~ log(Landmass_area) +
                       elev_mean + I(elev_mean^2) +
                       elev_range +
                       clim_pca_1 + I(clim_pca_1^2) +
                       clim_pca_2 + I(clim_pca_2^2) +
                       clim_pca_3 + I(clim_pca_3^2) +
                       clim_pca_4 + I(clim_pca_4^2) +
                       MDS1 + I(MDS1^2) + MDS2 +
                       I(MDS2^2) + MDS3 + I(MDS3^2) +
                       MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
                     data=data_Present,
                     na.action = na.fail)
SR_all_LGM <- lm(SR ~ log(Landmass_area) +
                   elev_mean + I(elev_mean^2) +
                   elev_range +
                   clim_pca_1 + I(clim_pca_1^2) +
                   clim_pca_2 + I(clim_pca_2^2) +
                   clim_pca_3 + I(clim_pca_3^2) +
                   clim_pca_4 + I(clim_pca_4^2) +
                   MDS1 + I(MDS1^2) + MDS2 +
                   I(MDS2^2) + MDS3 + I(MDS3^2) +
                   MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
                 data=data_LGM,
                 na.action = na.fail)
SR_all_MidHol <- lm(SR ~ log(Landmass_area) +
                      elev_mean + I(elev_mean^2) +
                      elev_range +
                      clim_pca_1 + I(clim_pca_1^2) +
                      clim_pca_2 + I(clim_pca_2^2) +
                      clim_pca_3 + I(clim_pca_3^2) +
                      clim_pca_4 + I(clim_pca_4^2) +
                      MDS1 + I(MDS1^2) + MDS2 +
                      I(MDS2^2) + MDS3 + I(MDS3^2) +
                      MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
                    data=data_MidHol,
                    na.action = na.fail)
summary(SR_all_Present)$r.squared
summary(SR_all_LGM)$r.squared
summary(SR_all_MidHol)$r.squared

# ---------- Functional richness ----------

# Load FR data for each climate
data_Present <- read.csv("Data/bird_FR_cells.csv") %>%
  left_join(cells_Present, by="cell_id")
data_LGM <- read.csv("Data/bird_FR_cells.csv") %>%
  left_join(cells_LGM, by="cell_id")
data_MidHol <- read.csv("Data/bird_FR_cells.csv") %>%
  left_join(cells_MidHol, by="cell_id")

# Run model for each climate, environment variables only
FR_env_Present <- lm(FR ~ log(Landmass_area) +
                       elev_mean + I(elev_mean^2) +
                       elev_range +
                       clim_pca_1 + I(clim_pca_1^2) +
                       clim_pca_2 + I(clim_pca_2^2) +
                       clim_pca_3 + I(clim_pca_3^2) +
                       clim_pca_4 + I(clim_pca_4^2),
                     data=data_Present,
                     na.action = na.fail)
FR_env_LGM <- lm(FR ~ log(Landmass_area) +
                   elev_mean + I(elev_mean^2) +
                   elev_range +
                   clim_pca_1 + I(clim_pca_1^2) +
                   clim_pca_2 + I(clim_pca_2^2) +
                   clim_pca_3 + I(clim_pca_3^2) +
                   clim_pca_4 + I(clim_pca_4^2),
                 data=data_LGM,
                 na.action = na.fail)
FR_env_MidHol <- lm(FR ~ log(Landmass_area) +
                      elev_mean + I(elev_mean^2) +
                      elev_range +
                      clim_pca_1 + I(clim_pca_1^2) +
                      clim_pca_2 + I(clim_pca_2^2) +
                      clim_pca_3 + I(clim_pca_3^2) +
                      clim_pca_4 + I(clim_pca_4^2),
                    data=data_MidHol,
                    na.action = na.fail)
summary(FR_env_Present)$r.squared
summary(FR_env_LGM)$r.squared
summary(FR_env_MidHol)$r.squared

# Run model for each climate, environment and phylobetadiversity
FR_all_Present <- lm(FR ~ log(Landmass_area) +
                       elev_mean + I(elev_mean^2) +
                       elev_range +
                       clim_pca_1 + I(clim_pca_1^2) +
                       clim_pca_2 + I(clim_pca_2^2) +
                       clim_pca_3 + I(clim_pca_3^2) +
                       clim_pca_4 + I(clim_pca_4^2) +
                       MDS1 + I(MDS1^2) + MDS2 +
                       I(MDS2^2) + MDS3 + I(MDS3^2) +
                       MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
                     data=data_Present,
                     na.action = na.fail)
FR_all_LGM <- lm(FR ~ log(Landmass_area) +
                   elev_mean + I(elev_mean^2) +
                   elev_range +
                   clim_pca_1 + I(clim_pca_1^2) +
                   clim_pca_2 + I(clim_pca_2^2) +
                   clim_pca_3 + I(clim_pca_3^2) +
                   clim_pca_4 + I(clim_pca_4^2) +
                   MDS1 + I(MDS1^2) + MDS2 +
                   I(MDS2^2) + MDS3 + I(MDS3^2) +
                   MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
                 data=data_LGM,
                 na.action = na.fail)
FR_all_MidHol <- lm(FR ~ log(Landmass_area) +
                      elev_mean + I(elev_mean^2) +
                      elev_range +
                      clim_pca_1 + I(clim_pca_1^2) +
                      clim_pca_2 + I(clim_pca_2^2) +
                      clim_pca_3 + I(clim_pca_3^2) +
                      clim_pca_4 + I(clim_pca_4^2) +
                      MDS1 + I(MDS1^2) + MDS2 +
                      I(MDS2^2) + MDS3 + I(MDS3^2) +
                      MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
                    data=data_MidHol,
                    na.action = na.fail)
summary(FR_all_Present)$r.squared
summary(FR_all_LGM)$r.squared
summary(FR_all_MidHol)$r.squared

# ----- MAMMALS  ----------

# Load cells with covariates
cells <- st_read("Data/global_cells.shp") %>%
  # Landmass area
  left_join(read.csv("Data/Raw/landmass_area_cells.csv"), "cell_id") %>%
  # Elevation (mean and range)
  left_join(read.csv("Data/elev_cells.csv"), "cell_id") %>%
  dplyr::select(-X) %>%
  # Phylobetadiversity NMDS axes
  left_join(read.csv("Data/mamm_pb_NMDS.csv"), "cell_id") %>%
  dplyr::select(-X)

# Present climate
cells_Present <- cells %>%
  left_join(read.csv("Data/clim_Present_cells.csv"), "cell_id") %>%
  dplyr::select(-X)

# LGM climate
cells_LGM <- cells %>%
  left_join(read.csv("Data/clim_LGM_cells.csv"), "cell_id") %>%
  dplyr::select(-X)

# MidHol climate
cells_MidHol <- cells %>%
  left_join(read.csv("Data/clim_MidHol_cells.csv"), "cell_id") %>%
  dplyr::select(-X)

# ---------- Species richness ----------

# Load SR data for each climate
data_Present <- read.csv("Data/mamm_SR_cells.csv") %>%
  left_join(cells_Present, by="cell_id")
data_LGM <- read.csv("Data/mamm_SR_cells.csv") %>%
  left_join(cells_LGM, by="cell_id")
data_MidHol <- read.csv("Data/mamm_SR_cells.csv") %>%
  left_join(cells_MidHol, by="cell_id")

# Run model for each climate, environment variables only
SR_env_Present <- lm(SR ~ log(Landmass_area) +
                       elev_mean + I(elev_mean^2) +
                       elev_range +
                       clim_pca_1 + I(clim_pca_1^2) +
                       clim_pca_2 + I(clim_pca_2^2) +
                       clim_pca_3 + I(clim_pca_3^2) +
                       clim_pca_4 + I(clim_pca_4^2),
                     data=data_Present,
                     na.action = na.fail)
SR_env_LGM <- lm(SR ~ log(Landmass_area) +
                   elev_mean + I(elev_mean^2) +
                   elev_range +
                   clim_pca_1 + I(clim_pca_1^2) +
                   clim_pca_2 + I(clim_pca_2^2) +
                   clim_pca_3 + I(clim_pca_3^2) +
                   clim_pca_4 + I(clim_pca_4^2),
                 data=data_LGM,
                 na.action = na.fail)
SR_env_MidHol <- lm(SR ~ log(Landmass_area) +
                      elev_mean + I(elev_mean^2) +
                      elev_range +
                      clim_pca_1 + I(clim_pca_1^2) +
                      clim_pca_2 + I(clim_pca_2^2) +
                      clim_pca_3 + I(clim_pca_3^2) +
                      clim_pca_4 + I(clim_pca_4^2),
                    data=data_MidHol,
                    na.action = na.fail)
summary(SR_env_Present)$r.squared
summary(SR_env_LGM)$r.squared
summary(SR_env_MidHol)$r.squared

# Run model for each climate, environment and phylobetadiversity
SR_all_Present <- lm(SR ~ log(Landmass_area) +
               elev_mean + I(elev_mean^2) +
               elev_range +
               clim_pca_1 + I(clim_pca_1^2) +
               clim_pca_2 + I(clim_pca_2^2) +
               clim_pca_3 + I(clim_pca_3^2) +
               clim_pca_4 + I(clim_pca_4^2) +
               MDS1 + I(MDS1^2) + MDS2 +
               I(MDS2^2) + MDS3 + I(MDS3^2) +
               MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
             data=data_Present,
             na.action = na.fail)
SR_all_LGM <- lm(SR ~ log(Landmass_area) +
                       elev_mean + I(elev_mean^2) +
                       elev_range +
                       clim_pca_1 + I(clim_pca_1^2) +
                       clim_pca_2 + I(clim_pca_2^2) +
                       clim_pca_3 + I(clim_pca_3^2) +
                       clim_pca_4 + I(clim_pca_4^2) +
                       MDS1 + I(MDS1^2) + MDS2 +
                       I(MDS2^2) + MDS3 + I(MDS3^2) +
                       MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
                     data=data_LGM,
                     na.action = na.fail)
SR_all_MidHol <- lm(SR ~ log(Landmass_area) +
                   elev_mean + I(elev_mean^2) +
                   elev_range +
                   clim_pca_1 + I(clim_pca_1^2) +
                   clim_pca_2 + I(clim_pca_2^2) +
                   clim_pca_3 + I(clim_pca_3^2) +
                   clim_pca_4 + I(clim_pca_4^2) +
                   MDS1 + I(MDS1^2) + MDS2 +
                   I(MDS2^2) + MDS3 + I(MDS3^2) +
                   MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
                 data=data_MidHol,
                 na.action = na.fail)
summary(SR_all_Present)$r.squared
summary(SR_all_LGM)$r.squared
summary(SR_all_MidHol)$r.squared

# ---------- Functional richness ----------

# Load FR data for each climate
data_Present <- read.csv("Data/mamm_FR_cells.csv") %>%
  left_join(cells_Present, by="cell_id")
data_LGM <- read.csv("Data/mamm_FR_cells.csv") %>%
  left_join(cells_LGM, by="cell_id")
data_MidHol <- read.csv("Data/mamm_FR_cells.csv") %>%
  left_join(cells_MidHol, by="cell_id")

# Run model for each climate, environment variables only
FR_env_Present <- lm(FR ~ log(Landmass_area) +
                       elev_mean + I(elev_mean^2) +
                       elev_range +
                       clim_pca_1 + I(clim_pca_1^2) +
                       clim_pca_2 + I(clim_pca_2^2) +
                       clim_pca_3 + I(clim_pca_3^2) +
                       clim_pca_4 + I(clim_pca_4^2),
                     data=data_Present,
                     na.action = na.fail)
FR_env_LGM <- lm(FR ~ log(Landmass_area) +
                   elev_mean + I(elev_mean^2) +
                   elev_range +
                   clim_pca_1 + I(clim_pca_1^2) +
                   clim_pca_2 + I(clim_pca_2^2) +
                   clim_pca_3 + I(clim_pca_3^2) +
                   clim_pca_4 + I(clim_pca_4^2),
                 data=data_LGM,
                 na.action = na.fail)
FR_env_MidHol <- lm(FR ~ log(Landmass_area) +
                      elev_mean + I(elev_mean^2) +
                      elev_range +
                      clim_pca_1 + I(clim_pca_1^2) +
                      clim_pca_2 + I(clim_pca_2^2) +
                      clim_pca_3 + I(clim_pca_3^2) +
                      clim_pca_4 + I(clim_pca_4^2),
                    data=data_MidHol,
                    na.action = na.fail)
summary(FR_env_Present)$r.squared
summary(FR_env_LGM)$r.squared
summary(FR_env_MidHol)$r.squared

# Run model for each climate, environment and phylobetadiversity
FR_all_Present <- lm(FR ~ log(Landmass_area) +
                       elev_mean + I(elev_mean^2) +
                       elev_range +
                       clim_pca_1 + I(clim_pca_1^2) +
                       clim_pca_2 + I(clim_pca_2^2) +
                       clim_pca_3 + I(clim_pca_3^2) +
                       clim_pca_4 + I(clim_pca_4^2) +
                       MDS1 + I(MDS1^2) + MDS2 +
                       I(MDS2^2) + MDS3 + I(MDS3^2) +
                       MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
                     data=data_Present,
                     na.action = na.fail)
FR_all_LGM <- lm(FR ~ log(Landmass_area) +
                   elev_mean + I(elev_mean^2) +
                   elev_range +
                   clim_pca_1 + I(clim_pca_1^2) +
                   clim_pca_2 + I(clim_pca_2^2) +
                   clim_pca_3 + I(clim_pca_3^2) +
                   clim_pca_4 + I(clim_pca_4^2) +
                   MDS1 + I(MDS1^2) + MDS2 +
                   I(MDS2^2) + MDS3 + I(MDS3^2) +
                   MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
                 data=data_LGM,
                 na.action = na.fail)
FR_all_MidHol <- lm(FR ~ log(Landmass_area) +
                      elev_mean + I(elev_mean^2) +
                      elev_range +
                      clim_pca_1 + I(clim_pca_1^2) +
                      clim_pca_2 + I(clim_pca_2^2) +
                      clim_pca_3 + I(clim_pca_3^2) +
                      clim_pca_4 + I(clim_pca_4^2) +
                      MDS1 + I(MDS1^2) + MDS2 +
                      I(MDS2^2) + MDS3 + I(MDS3^2) +
                      MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
                    data=data_MidHol,
                    na.action = na.fail)
summary(FR_all_Present)$r.squared
summary(FR_all_LGM)$r.squared
summary(FR_all_MidHol)$r.squared
