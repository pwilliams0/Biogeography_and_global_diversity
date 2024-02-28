#===================================#
# Maps of global diversity patterns #
#===================================#

library(tidyverse)
library(sf)
library(colorspace)
library(rnaturalearth)
library(rnaturalearthdata)
sf_use_s2(FALSE)
# Set color palette
colorspace::sequential_hcl(n = 12, h = 0, c = c(8, 80, NA),
                           l = c(5, 95), power = 1,
                           register = "Biogeography")

cells <- st_read("Data/global_cells.shp")
world <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_crop(c(xmin=-180, xmax=180, ymin=-60, ymax=85))

# ----- BIRDS ----------

# ---------- Species richness (Supplementary Fig. 3a) ----------
# Load data
data <- cells %>%
  right_join(read.csv("Data/bird_SR_cells.csv"), "cell_id")

# Save figure
svg("Results/SR_map_bird.svg",
    width = 2.5, height = 1.5)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data=data, aes(color = SR, fill = SR),
          lwd=.1, linewidth=.1) +
  scale_color_continuous_sequential("Biogeography") +
  scale_fill_continuous_sequential("Biogeography") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.17, "cm"),
        legend.key.width = unit(1.1,"cm"))
dev.off()

# ---------- Phylogenetic diversity (Supplementary Fig. 3b) ----------
# Load data
data <- cells %>%
  right_join(read.csv("Data/bird_PD_cells.csv"), "cell_id")

# Save figure
svg("Results/PD_map_bird.svg",
    width = 2.5, height = 1.5)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data=data, aes(color = PD, fill = PD),
          lwd=.1, linewidth=.1) +
  scale_color_continuous_sequential("Biogeography") +
  scale_fill_continuous_sequential("Biogeography") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.17, "cm"),
        legend.key.width = unit(1.1,"cm"))
dev.off()

# ---------- Functional richness (Supplementary Fig. 3c) ----------
# Load data
data <- cells %>%
  right_join(read.csv("Data/bird_FR_cells.csv"), "cell_id")

# Save figure
svg("Results/FR_map_bird.svg",
    width = 2.5, height = 1.5)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data=data, aes(color = FR, fill = FR),
          lwd=.1, linewidth=.1) +
  scale_color_continuous_sequential("Biogeography") +
  scale_fill_continuous_sequential("Biogeography") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.17, "cm"),
        legend.key.width = unit(1.1,"cm"))
dev.off()

# ---------- Mean func beta turnover (Fig. 4a) ----------
# Load data
data <- cells %>%
  right_join(read.csv("Data/bird_mean_fb.csv"), "cell_id")

# Save figure
svg("Results/FB_map_bird.svg",
    width = 3, height = 1.8)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data=data, aes(color = mean_fb, fill = mean_fb),
          lwd=.1, linewidth=.1) +
  scale_color_continuous_sequential("Biogeography", limits=c(0,.51)) +
  scale_fill_continuous_sequential("Biogeography", limits=c(0,.51)) +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.2, "cm"),
        legend.key.width = unit(1.3,"cm"))
dev.off()

# ----- MAMMALS ----------

# ---------- Species richness (Supplementary Fig. 3d) ----------
# Load data
data <- cells %>%
  right_join(read.csv("Data/mamm_SR_cells.csv"), "cell_id")

# Save figure
svg("Results/SR_map_mamm.svg",
    width = 2.5, height = 1.5)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data=data, aes(color = SR, fill = SR),
          lwd=.1, linewidth=.1) +
  scale_color_continuous_sequential("Biogeography") +
  scale_fill_continuous_sequential("Biogeography") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.17, "cm"),
        legend.key.width = unit(1.1,"cm"))
dev.off()

# ---------- Phylogenetic diversity (Supplementary Fig. 3e) ----------
# Load data
data <- cells %>%
  right_join(read.csv("Data/mamm_PD_cells.csv"), "cell_id")

# Save figure
svg("Results/PD_map_mamm.svg",
    width = 2.5, height = 1.5)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data=data, aes(color = PD, fill = PD),
          lwd=.1, linewidth=.1) +
  scale_color_continuous_sequential("Biogeography") +
  scale_fill_continuous_sequential("Biogeography") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.17, "cm"),
        legend.key.width = unit(1.1,"cm"))
dev.off()

# ---------- Functional richness (Supplementary Fig. 3f) ----------
# Load data
data <- cells %>%
  right_join(read.csv("Data/mamm_FR_cells.csv"), "cell_id")

# Save figure
svg("Results/FR_map_mamm.svg",
    width = 2.5, height = 1.5)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data=data, aes(color = FR, fill = FR),
          lwd=.1, linewidth=.1) +
  scale_color_continuous_sequential("Biogeography") +
  scale_fill_continuous_sequential("Biogeography") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.17, "cm"),
        legend.key.width = unit(1.1,"cm"))
dev.off()

# ---------- Mean func beta turnover (Fig. 4b) ----------
# Load data
data <- cells %>%
  right_join(read.csv("Data/mamm_mean_fb.csv"), "cell_id")

# Save figure
svg("Results/FB_map_mamm.svg",
    width = 3, height = 1.8)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data=data, aes(color = mean_fb, fill = mean_fb),
          lwd=.1, linewidth=.1) +
  scale_color_continuous_sequential("Biogeography", limits=c(0,.51)) +
  scale_fill_continuous_sequential("Biogeography", limits=c(0,.51)) +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.2, "cm"),
        legend.key.width = unit(1.3,"cm"))
dev.off()

# ----- BATS ----------

# ---------- Species richness (Supplementary Fig. 3g) ----------
# Load data
data <- cells %>%
  right_join(read.csv("Data/bat_SR_cells.csv"), "cell_id")

# Save figure
svg("Results/SR_map_bat.svg",
    width = 2.5, height = 1.5)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data=data, aes(color = SR, fill = SR),
          lwd=.1, linewidth=.1) +
  scale_color_continuous_sequential("Biogeography") +
  scale_fill_continuous_sequential("Biogeography") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.17, "cm"),
        legend.key.width = unit(1.1,"cm"))
dev.off()

# ---------- Phylogenetic diversity (Supplementary Fig. 3h) ----------
# Load data
data <- cells %>%
  right_join(read.csv("Data/bat_PD_cells.csv"), "cell_id")

# Save figure
svg("Results/PD_map_bat.svg",
    width = 2.5, height = 1.5)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data=data, aes(color = PD, fill = PD),
          lwd=.1, linewidth=.1) +
  scale_color_continuous_sequential("Biogeography") +
  scale_fill_continuous_sequential("Biogeography") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.17, "cm"),
        legend.key.width = unit(1.1,"cm"))
dev.off()

# ---------- Functional richness (Fig. 6b) ----------
# Load data
data <- cells %>%
  right_join(read.csv("Data/bat_FR_cells.csv"), "cell_id") %>%
  na.omit()

# Save figure
svg("Results/FR_map_bat.svg",
    width = 3, height = 1.8)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data=data, aes(color = FR, fill = FR),
          lwd=.1, linewidth=.1) +
  scale_color_continuous_sequential("Biogeography") +
  scale_fill_continuous_sequential("Biogeography") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.2, "cm"),
        legend.key.width = unit(1.3,"cm"))
dev.off()

# ---------- Mean func beta turnover (Supplementary Fig. 3i) ----------
# Load data
data <- cells %>%
  right_join(read.csv("Data/bat_mean_fb.csv"), "cell_id") %>%
  na.omit()

# Save figure
svg("Results/FB_map_bat.svg",
    width = 2.5, height = 1.5)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data=data, aes(color = mean_fb, fill = mean_fb),
          lwd=.1, linewidth=.1) +
  scale_color_continuous_sequential("Biogeography") +
  scale_fill_continuous_sequential("Biogeography") +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.17, "cm"),
        legend.key.width = unit(1.1,"cm"))
dev.off()
