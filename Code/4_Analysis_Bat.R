#===================================================#
# Analyze species richness, phylogenetic diversity, #
#   functional richness, and mean functional beta   #
#            diversity turnover for bats            #
#===================================================#

library(tidyverse)
library(sf)
library(colorspace)
library(rnaturalearth)
library(rnaturalearthdata)
library(spdep)
sf_use_s2(FALSE) # Avoid some issues of invalid geometries

# Load outlines of landmasses for mapping
world <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_crop(c(xmin=-180, xmax=180, ymin=-60, ymax=85))

# ----- LOAD CELLS WITH COVARIATES ----------

cells <- st_read("Data/global_cells.shp") %>%
  # Landmass area
  left_join(read.csv("Data/Raw/landmass_area_cells.csv"), "cell_id") %>%
  # Elevation (mean and range)
  left_join(read.csv("Data/elev_cells.csv"), "cell_id") %>%
  dplyr::select(-X) %>%
  # Phylobetadiversity NMDS axes
  left_join(read.csv("Data/bat_pb_NMDS.csv"), "cell_id") %>%
  dplyr::select(-X) %>%
  # Present climate
  left_join(read.csv("Data/clim_Present_cells.csv"), "cell_id") %>%
  dplyr::select(-X)
  
# ----- SPECIES RICHNESS ----------

# ---------- Run models, save R2 ----------

# Load data
data <- read.csv("Data/bat_SR_cells.csv") %>%
  left_join(cells, by="cell_id")

# Run model, all covariates
SR_all <- lm(SR ~ log(Landmass_area) +
               elev_mean + I(elev_mean^2) +
               elev_range +
               clim_pca_1 + I(clim_pca_1^2) +
               clim_pca_2 + I(clim_pca_2^2) +
               clim_pca_3 + I(clim_pca_3^2) +
               clim_pca_4 + I(clim_pca_4^2) +
               MDS1 + I(MDS1^2) + MDS2 +
               I(MDS2^2) + MDS3 + I(MDS3^2) +
               MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
             data=data,
             na.action = na.fail)
# Run model, environment
SR_env <- lm(SR ~ log(Landmass_area) +
               elev_mean + I(elev_mean^2) +
               elev_range +
               clim_pca_1 + I(clim_pca_1^2) +
               clim_pca_2 + I(clim_pca_2^2) +
               clim_pca_3 + I(clim_pca_3^2) +
               clim_pca_4 + I(clim_pca_4^2),
             data=data,
             na.action = na.fail)
# Run model, phylobetadiversity
SR_pb <- lm(SR ~ MDS1 + I(MDS1^2) + MDS2 +
              I(MDS2^2) + MDS3 + I(MDS3^2) +
              MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
            data=data,
            na.action = na.fail)
# Save R2 of each model
R2_all <- summary(SR_all)$adj.r.squared
R2_env <- summary(SR_env)$adj.r.squared
R2_pb <- summary(SR_pb)$adj.r.squared
# Phylobetadiversity only
R2_pb_only <- R2_all - R2_env
# Environment only
R2_env_only <- R2_all - R2_pb
# Shared
R2_shared <- R2_all - R2_pb_only - R2_env_only

# Save results
bat_SR_results <- data.frame(matrix(ncol = 4, nrow = 3))
colnames(bat_SR_results) <- c("Taxon","Metric","varpart","R2")
bat_SR_results$Taxon <- "bat"
bat_SR_results$Metric <- "SR"
bat_SR_results$varpart <- factor(c("Shared","Environment only",
                                    "Phylobetadiversity only"),
                                  levels=c("Shared","Environment only",
                                           "Phylobetadiversity only"))
bat_SR_results$R2 <- c(R2_shared, R2_env_only, R2_pb_only)
write.csv(bat_SR_results, "Results/bat_SR_results.csv")

# ---------- Map differences in residuals (Extended Data Fig. 1g) ----------

# Difference with/without phylobetadiversity
res_df <- data.frame(
  res = abs(resid(SR_env)) - abs(resid(SR_all)),
  cell_id = as.numeric(data$cell_id))
res_map <- cells %>% right_join(res_df, by="cell_id")
# Save figure
svg("Results/diffres_bat_SR.svg",
    width = 2.5, height = 1.5)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data = res_map, aes(color = res, fill = res),
          lwd=.1, linewidth=.1) +
  scale_colour_continuous_diverging(palette="Blue-Red") +
  scale_fill_continuous_diverging(palette="Blue-Red") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.17, "cm"),
        legend.key.width = unit(1.1,"cm"))
dev.off()

# ----- PHYLOGENETIC DIVERSITY ----------

# ---------- Run models, save R2 ----------

# Load data
data <- read.csv("Data/bat_PD_cells.csv") %>%
  left_join(cells, by="cell_id")

# Run model, all covariates
PD_all <- lm(PD ~ log(Landmass_area) +
               elev_mean + I(elev_mean^2) +
               elev_range +
               clim_pca_1 + I(clim_pca_1^2) +
               clim_pca_2 + I(clim_pca_2^2) +
               clim_pca_3 + I(clim_pca_3^2) +
               clim_pca_4 + I(clim_pca_4^2) +
               MDS1 + I(MDS1^2) + MDS2 +
               I(MDS2^2) + MDS3 + I(MDS3^2) +
               MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
             data=data,
             na.action = na.fail)
# Run model, environment
PD_env <- lm(PD ~ log(Landmass_area) +
               elev_mean + I(elev_mean^2) +
               elev_range +
               clim_pca_1 + I(clim_pca_1^2) +
               clim_pca_2 + I(clim_pca_2^2) +
               clim_pca_3 + I(clim_pca_3^2) +
               clim_pca_4 + I(clim_pca_4^2),
             data=data,
             na.action = na.fail)
# Run model, phylobetadiversity
PD_pb <- lm(PD ~ MDS1 + I(MDS1^2) + MDS2 +
              I(MDS2^2) + MDS3 + I(MDS3^2) +
              MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
            data=data,
            na.action = na.fail)
# Save R2 of each model
R2_all <- summary(PD_all)$adj.r.squared
R2_env <- summary(PD_env)$adj.r.squared
R2_pb <- summary(PD_pb)$adj.r.squared
# Phylobetadiversity only
R2_pb_only <- R2_all - R2_env
# Environment only
R2_env_only <- R2_all - R2_pb
# Shared
R2_shared <- R2_all - R2_pb_only - R2_env_only

# Save results
bat_PD_results <- data.frame(matrix(ncol = 4, nrow = 3))
colnames(bat_PD_results) <- c("Taxon","Metric","varpart","R2")
bat_PD_results$Taxon <- "bat"
bat_PD_results$Metric <- "PD"
bat_PD_results$varpart <- factor(c("Shared","Environment only",
                                    "Phylobetadiversity only"),
                                  levels=c("Shared","Environment only",
                                           "Phylobetadiversity only"))
bat_PD_results$R2 <- c(R2_shared, R2_env_only, R2_pb_only)
write.csv(bat_PD_results, "Results/bat_PD_results.csv")

# ---------- Map differences in residuals (Extended Data Fig. 1h) ----------

# Difference with/without phylobetadiversity
res_df <- data.frame(
  res = abs(resid(PD_env)) - abs(resid(PD_all)),
  cell_id = as.numeric(data$cell_id))
res_map <- cells %>% right_join(res_df, by="cell_id")
# Save figure
svg("Results/diffres_bat_PD.svg",
    width = 2.5, height = 1.5)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data = res_map, aes(color = res, fill = res),
          lwd=.1, linewidth=.1) +
  scale_colour_continuous_diverging(palette="Blue-Red") +
  scale_fill_continuous_diverging(palette="Blue-Red") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.17, "cm"),
        legend.key.width = unit(1.1,"cm"))
dev.off()

# ----- FUNCTIONAL RICHNESS ----------

# ---------- Run models, save R2 ----------

# Load data
data <- read.csv("Data/bat_FR_cells.csv") %>%
  na.omit() %>%
  left_join(cells, by="cell_id")

# Run model, all covariates
FR_all <- lm(FR ~ log(Landmass_area) +
               elev_mean + I(elev_mean^2) +
               elev_range +
               clim_pca_1 + I(clim_pca_1^2) +
               clim_pca_2 + I(clim_pca_2^2) +
               clim_pca_3 + I(clim_pca_3^2) +
               clim_pca_4 + I(clim_pca_4^2) +
               MDS1 + I(MDS1^2) + MDS2 +
               I(MDS2^2) + MDS3 + I(MDS3^2) +
               MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
             data=data,
             na.action = na.fail)
# Run model, environment
FR_env <- lm(FR ~ log(Landmass_area) +
               elev_mean + I(elev_mean^2) +
               elev_range +
               clim_pca_1 + I(clim_pca_1^2) +
               clim_pca_2 + I(clim_pca_2^2) +
               clim_pca_3 + I(clim_pca_3^2) +
               clim_pca_4 + I(clim_pca_4^2),
             data=data,
             na.action = na.fail)
# Run model, phylobetadiversity
FR_pb <- lm(FR ~ MDS1 + I(MDS1^2) + MDS2 +
              I(MDS2^2) + MDS3 + I(MDS3^2) +
              MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
            data=data,
            na.action = na.fail)
# Save R2 of each model
R2_all <- summary(FR_all)$adj.r.squared
R2_env <- summary(FR_env)$adj.r.squared
R2_pb <- summary(FR_pb)$adj.r.squared
# Phylobetadiversity only
R2_pb_only <- R2_all - R2_env
# Environment only
R2_env_only <- R2_all - R2_pb
# Shared
R2_shared <- R2_all - R2_pb_only - R2_env_only

# Save results
bat_FR_results <- data.frame(matrix(ncol = 4, nrow = 3))
colnames(bat_FR_results) <- c("Taxon","Metric","varpart","R2")
bat_FR_results$Taxon <- "bat"
bat_FR_results$Metric <- "FR"
bat_FR_results$varpart <- factor(c("Shared","Environment only",
                                    "Phylobetadiversity only"),
                                  levels=c("Shared","Environment only",
                                           "Phylobetadiversity only"))
bat_FR_results$R2 <- c(R2_shared, R2_env_only, R2_pb_only)
write.csv(bat_FR_results, "Results/bat_FR_results.csv")

# ---------- Map differences in residuals (Extended Data Fig. 1i) ----------

# Difference with/without phylobetadiversity
res_df <- data.frame(
  res = abs(resid(FR_env)) - abs(resid(FR_all)),
  cell_id = as.numeric(data$cell_id))
res_map <- cells %>% right_join(res_df, by="cell_id")
# Save figure
svg("Results/diffres_bat_FR.svg",
    width = 2.5, height = 1.5)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data = res_map, aes(color = res, fill = res),
          lwd=.1, linewidth=.1) +
  scale_colour_continuous_diverging(palette="Blue-Red") +
  scale_fill_continuous_diverging(palette="Blue-Red") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.17, "cm"),
        legend.key.width = unit(1.1,"cm"))
dev.off()

# ----- MEAN FUNC BETA TURNOVER ----------

# ---------- Run models, save R2 ----------

# Load data
data <- read.csv("Data/bat_mean_fb.csv") %>%
  left_join(cells, by="cell_id")

# Run model, all covariates
FB_all <- lm(mean_fb ~ log(Landmass_area) +
               elev_mean + I(elev_mean^2) +
               elev_range +
               clim_pca_1 + I(clim_pca_1^2) +
               clim_pca_2 + I(clim_pca_2^2) +
               clim_pca_3 + I(clim_pca_3^2) +
               clim_pca_4 + I(clim_pca_4^2) +
               MDS1 + I(MDS1^2) + MDS2 +
               I(MDS2^2) + MDS3 + I(MDS3^2) +
               MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
             data=data,
             na.action = na.fail)
# Run model, environment
FB_env <- lm(mean_fb ~ log(Landmass_area) +
               elev_mean + I(elev_mean^2) +
               elev_range +
               clim_pca_1 + I(clim_pca_1^2) +
               clim_pca_2 + I(clim_pca_2^2) +
               clim_pca_3 + I(clim_pca_3^2) +
               clim_pca_4 + I(clim_pca_4^2),
             data=data,
             na.action = na.fail)
# Run model, phylobetadiversity
FB_pb <- lm(mean_fb ~ MDS1 + I(MDS1^2) + MDS2 +
              I(MDS2^2) + MDS3 + I(MDS3^2) +
              MDS1:MDS2 + MDS1:MDS3 + MDS2:MDS3,
            data=data,
            na.action = na.fail)
# Save R2 of each model
R2_all <- summary(FB_all)$adj.r.squared
R2_env <- summary(FB_env)$adj.r.squared
R2_pb <- summary(FB_pb)$adj.r.squared
# Phylobetadiversity only
R2_pb_only <- R2_all - R2_env
# Environment only
R2_env_only <- R2_all - R2_pb
# Shared
R2_shared <- R2_all - R2_pb_only - R2_env_only

# Save results
bat_FB_results <- data.frame(matrix(ncol = 4, nrow = 3))
colnames(bat_FB_results) <- c("Taxon","Metric","varpart","R2")
bat_FB_results$Taxon <- "bat"
bat_FB_results$Metric <- "FB"
bat_FB_results$varpart <- factor(c("Shared","Environment only",
                                    "Phylobetadiversity only"),
                                  levels=c("Shared","Environment only",
                                           "Phylobetadiversity only"))
bat_FB_results$R2 <- c(R2_shared, R2_env_only, R2_pb_only)
write.csv(bat_FB_results, "Results/bat_FB_results.csv")

# ---------- Map differences in residuals (Extended Data Fig. 1j) ----------

# Difference with/without phylobetadiversity
res_df <- data.frame(
  res = abs(resid(FB_env)) - abs(resid(FB_all)),
  cell_id = as.numeric(data$cell_id))
res_map <- cells %>% right_join(res_df, by="cell_id")
# Save figure
svg("Results/diffres_bat_FB.svg",
    width = 2.5, height = 1.5)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data = res_map, aes(color = res, fill = res),
          lwd=.1, linewidth=.1) +
  scale_colour_continuous_diverging(palette="Blue-Red") +
  scale_fill_continuous_diverging(palette="Blue-Red") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.17, "cm"),
        legend.key.width = unit(1.1,"cm"))
dev.off()
