#============================#
# Phylobetadiversity figures #
#============================#

library(tidyverse)
library(sf)
library(colorspace)
library(rnaturalearth)
library(rnaturalearthdata)
library(vegan)
sf_use_s2(FALSE)

# Set palette for realms
pal <- c("#6929c4","#1192e8","#ee538b","#9f1853",
         "#fa4d56","#570408","#198038","#002d9c",
         "#005d5d","#b28600","#009d9a")

# Load cells for maps and for Realm
cells <- st_read("Data/global_cells.shp")

# ----- REALM MAP (Fig. 1 legend) ----------

# Load coastline for map
world <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_crop(c(xmin=-180, xmax=180, ymin=-60, ymax=85))

# Save map
svg("Results/realm_map.svg",
      width = 2.3, height = 1)
ggplot() +
  geom_sf(data = world, colour="grey85", lwd=.1, linewidth=.1) +
  geom_sf(data = cells, aes(color = Realm, fill = Realm),
          lwd=.1, linewidth=.1) +
  scale_color_manual(values=pal) +
  scale_fill_manual(values=pal) +
  theme_void() +
  theme(legend.position = "none")
dev.off()

# ----- NMDS PLOTS, 2 AXES ----------

# ---------- Birds (Fig. 1a) ----------

# Load phylobetadiversity turnover distance matrix
pb_dist <- readRDS("Data/bird_pb_dist.RDS")

# Run NMDS, 2 axes (k=2)
set.seed(123)
start <- Sys.time()
pb_mds <- vegan::metaMDS(pb_dist, k=2) # 15 minutes
Sys.time() - start
pb_mds$stress # stress = 0.1873869 
pb_scores <- as.data.frame(pb_mds$points)  # Extract the site scores and convert to a data.frame
pb_scores$cell_id <- as.integer(rownames(pb_scores)) # create a column of site names, from the rownames of data.scores
# Join with cells to get Realm
# Re-scale axes
pb_plot <- cells %>%
  right_join(pb_scores, by="cell_id") %>%
  mutate(cell_id = as.character(cell_id),
         MDS1 = MDS1-min(MDS1),
         MDS1 = MDS1/max(MDS1),
         MDS2 = MDS2-min(MDS2),
         MDS2 = MDS2/max(MDS2))

# Save NMDS plot
svg("Results/pb_bird_nmds.svg",
    width = 2, height = 2)
ggplot(pb_plot) +
  aes(x = MDS1, y = MDS2, color = Realm) +
  geom_point(shape = 19, size=.3) +
  scale_color_manual(values=pal) +
  coord_fixed(1) +
  theme_classic() +
  labs(x = "NMDS 1", y ="NMDS 2") +
  theme(legend.position = "none",
        axis.text = element_text(color="black",
                                 size=6),
        axis.title = element_text(size=7),
        line = element_line(size = .25))
dev.off()

# ---------- Mammals (Fig. 1b) ----------

# Load phylobetadiversity turnover distance matrix
pb_dist <- readRDS("Data/mamm_pb_dist.RDS")

# Run NMDS, 2 axes (k=2)
set.seed(123)
start <- Sys.time()
pb_mds <- vegan::metaMDS(pb_dist, k=2) # 14 minutes
Sys.time() - start
pb_mds$stress # stress = 0.1952713
pb_scores <- as.data.frame(pb_mds$points)  # Extract the site scores and convert to a data.frame
pb_scores$cell_id <- as.integer(rownames(pb_scores)) # create a column of site names, from the rownames of data.scores
# Join with cells to get Realm
# Re-scale axes
# Flip axes to be more interpretable
pb_plot <- cells %>%
  right_join(pb_scores, by="cell_id") %>%
  mutate(cell_id = as.character(cell_id),
         MDS1 = -MDS1,
         MDS1 = MDS1-min(MDS1),
         MDS1 = MDS1/max(MDS1),
         MDS2 = MDS2-min(MDS2),
         MDS2 = MDS2/max(MDS2))

# Save NMDS plot
svg("Results/pb_mamm_nmds.svg",
    width = 2, height = 2)
ggplot(pb_plot) +
  aes(x = MDS1, y = MDS2, color = Realm) +
  geom_point(shape = 19, size=.3) +
  scale_color_manual(values=pal) +
  coord_fixed(1) +
  theme_classic() +
  labs(x = "NMDS 1", y ="NMDS 2") +
  theme(legend.position = "none",
        axis.text = element_text(color="black",
                                 size=6),
        axis.title = element_text(size=7),
        line = element_line(size = .25))
dev.off()

# ----- NMDS PLOTS, 3 AXES ----------

# ---------- Birds (Extended Data Fig. 3a-c) ----------

# Load phylobetadiversity NMDS coordinates
pb_NMDS <- read.csv("Data/bird_pb_NMDS.csv")

# Join with cells to get Realm
# Flip axes so x and y axes can be same for all plots
pb_plot <- cells %>%
  right_join(pb_NMDS, by="cell_id") %>%
  mutate(cell_id = as.character(cell_id),
         MDS2 = -MDS2)

# Save NMDS plot, axes 1 and 2
svg("Results/pb_bird_nmds_12.svg",
    width = 2, height = 2)
ggplot(pb_plot) +
  aes(x = MDS1, y = MDS2, color = Realm) +
  geom_point(shape = 19, size=.3) +
  scale_color_manual(values=pal) +
  coord_cartesian(xlim=c(-.51,.38), ylim=c(-.51,.38)) +
  theme_classic() +
  labs(x = "NMDS 1", y ="NMDS 2") +
  theme(legend.position = "none",
        axis.text = element_text(color="black",
                                 size=6),
        axis.title = element_text(size=7),
        line = element_line(size = .25))
dev.off()

# Save NMDS plot, axes 1 and 3
svg("Results/pb_bird_nmds_13.svg",
    width = 2, height = 2)
ggplot(pb_plot) +
  aes(x = MDS1, y = MDS3, color = Realm) +
  geom_point(shape = 19, size=.3) +
  scale_color_manual(values=pal) +
  coord_cartesian(xlim=c(-.51,.38), ylim=c(-.51,.38)) +
  theme_classic() +
  labs(x = "NMDS 1", y ="NMDS 3") +
  theme(legend.position = "none",
        axis.text = element_text(color="black",
                                 size=6),
        axis.title = element_text(size=7),
        line = element_line(size = .25))
dev.off()

# Save NMDS plot, axes 2 and 3
svg("Results/pb_bird_nmds_23.svg",
    width = 2, height = 2)
ggplot(pb_plot) +
  aes(x = MDS2, y = MDS3, color = Realm) +
  geom_point(shape = 19, size=.3) +
  scale_color_manual(values=pal) +
  coord_cartesian(xlim=c(-.51,.38), ylim=c(-.51,.38)) +
  theme_classic() +
  labs(x = "NMDS 2", y ="NMDS 3") +
  theme(legend.position = "none",
        axis.text = element_text(color="black",
                                 size=6),
        axis.title = element_text(size=7),
        line = element_line(size = .25))
dev.off()

# ---------- Mammals (Extended Data Fig. 3d-f)  ----------

# Load phylobetadiversity NMDS coordinates
pb_NMDS <- read.csv("Data/mamm_pb_NMDS.csv")

# Join with cells to get Realm
# Flip axes so x and y axes can be same for all plots
pb_plot <- cells %>%
  right_join(pb_NMDS, by="cell_id") %>%
  mutate(cell_id = as.character(cell_id),
         MDS2 = -MDS2,
         MDS3 = -MDS3)

# Save NMDS plot, axes 1 and 2
svg("Results/pb_mamm_nmds_12.svg",
    width = 2, height = 2)
ggplot(pb_plot) +
  aes(x = MDS1, y = MDS2, color = Realm) +
  geom_point(shape = 19, size=.3) +
  scale_color_manual(values=pal) +
  coord_cartesian(xlim=c(-.61,.33), ylim=c(-.61,.33)) +
  theme_classic() +
  labs(x = "NMDS 1", y ="NMDS 2") +
  theme(legend.position = "none",
        axis.text = element_text(color="black",
                                 size=6),
        axis.title = element_text(size=7),
        line = element_line(size = .25))
dev.off()

# Save NMDS plot, axes 1 and 3
svg("Results/pb_mamm_nmds_13.svg",
    width = 2, height = 2)
ggplot(pb_plot) +
  aes(x = MDS1, y = MDS3, color = Realm) +
  geom_point(shape = 19, size=.3) +
  scale_color_manual(values=pal) +
  coord_cartesian(xlim=c(-.61,.33), ylim=c(-.61,.33)) +
  theme_classic() +
  labs(x = "NMDS 1", y ="NMDS 3") +
  theme(legend.position = "none",
        axis.text = element_text(color="black",
                                 size=6),
        axis.title = element_text(size=7),
        line = element_line(size = .25))
dev.off()

# Save NMDS plot, axes 2 and 3
svg("Results/pb_mamm_nmds_23.svg",
    width = 2, height = 2)
ggplot(pb_plot) +
  aes(x = MDS2, y = MDS3, color = Realm) +
  geom_point(shape = 19, size=.3) +
  scale_color_manual(values=pal) +
  coord_cartesian(xlim=c(-.61,.33), ylim=c(-.61,.33)) +
  theme_classic() +
  labs(x = "NMDS 2", y ="NMDS 3") +
  theme(legend.position = "none",
        axis.text = element_text(color="black",
                                 size=6),
        axis.title = element_text(size=7),
        line = element_line(size = .25))
dev.off()

# ---------- Bats (Extended Data Fig. 3g-i)  ----------

# Load phylobetadiversity NMDS coordinates
pb_NMDS <- read.csv("Data/bat_pb_NMDS.csv")

# Join with cells to get Realm
# Flip axes to be more interpretable and consistent with other taxa
pb_plot <- cells %>%
  right_join(pb_NMDS, by="cell_id") %>%
  mutate(cell_id = as.character(cell_id),
         MDS1 = -MDS1,
         MDS2 = -MDS2)

# Save NMDS plot, axes 1 and 2
svg("Results/pb_bat_nmds_12.svg",
    width = 2, height = 2)
ggplot(pb_plot) +
  aes(x = MDS1, y = MDS2, color = Realm) +
  geom_point(shape = 19, size=.3) +
  scale_color_manual(values=pal) +
  coord_cartesian(xlim=c(-.61,.39), ylim=c(-.61,.39)) +
  theme_classic() +
  labs(x = "NMDS 1", y ="NMDS 2") +
  theme(legend.position = "none",
        axis.text = element_text(color="black",
                                 size=6),
        axis.title = element_text(size=7),
        line = element_line(size = .25))
dev.off()

# Save NMDS plot, axes 1 and 3
svg("Results/pb_bat_nmds_13.svg",
    width = 2, height = 2)
ggplot(pb_plot) +
  aes(x = MDS1, y = MDS3, color = Realm) +
  geom_point(shape = 19, size=.3) +
  scale_color_manual(values=pal) +
  coord_cartesian(xlim=c(-.61,.39), ylim=c(-.61,.39)) +
  theme_classic() +
  labs(x = "NMDS 1", y ="NMDS 3") +
  theme(legend.position = "none",
        axis.text = element_text(color="black",
                                 size=6),
        axis.title = element_text(size=7),
        line = element_line(size = .25))
dev.off()

# Save NMDS plot, axes 2 and 3
svg("Results/pb_bat_nmds_23.svg",
    width = 2, height = 2)
ggplot(pb_plot) +
  aes(x = MDS2, y = MDS3, color = Realm) +
  geom_point(shape = 19, size=.3) +
  scale_color_manual(values=pal) +
  coord_cartesian(xlim=c(-.61,.39), ylim=c(-.61,.39)) +
  theme_classic() +
  labs(x = "NMDS 2", y ="NMDS 3") +
  theme(legend.position = "none",
        axis.text = element_text(color="black",
                                 size=6),
        axis.title = element_text(size=7),
        line = element_line(size = .25))
dev.off()

# ----- SCREE PLOTS ----------

# ---------- Birds (Extended Data Fig. 5a) ----------

# Load phylobetadiversity turnover distance matrix
pb_dist <- readRDS("Data/bird_pb_dist.RDS")

# Calculate stress values for NDMS with 1-6 axes
Stress <- c()
for(i in 1:6){
  set.seed(123)
  func_mds <- vegan::metaMDS(pb_dist, k=i)
  Stress[i] <- func_mds$stress
}

# Save scree plot
svg("Results/pb_bird_scree.svg",
    width = 2.5, height = 1.5, pointsize = 7)
par(mar=c(3,3,1,.5))
plot(1:6, Stress, axes=FALSE, xlab="", ylab="", ylim=c(0,.4))
abline(h=.2, lty="longdash", col="red")
lines(1:6, Stress, lwd=.5)
axis(1, lwd=.5)
axis(2, lwd=.5)
title(ylab="Stress", xlab="Number of Dimensions", line=2)
dev.off()

# ---------- Mammals (Extended Data Fig. 5b)  ----------

# Load phylobetadiversity turnover distance matrix
pb_dist <- readRDS("Data/mamm_pb_dist.RDS")

# Calculate stress values for NDMS with 1-6 axes
Stress <- c()
for(i in 1:6){
  set.seed(123)
  func_mds <- vegan::metaMDS(pb_dist, k=i)
  Stress[i] <- func_mds$stress
}

# Save scree plot
svg("Results/pb_mamm_scree.svg",
    width = 2.5, height = 1.5, pointsize = 7)
par(mar=c(3,3,1,.5))
plot(1:6, Stress, axes=FALSE, xlab="", ylab="", ylim=c(0,.4))
abline(h=.2, lty="longdash", col="red")
lines(1:6, Stress, lwd=.5)
axis(1, lwd=.5)
axis(2, lwd=.5)
title(ylab="Stress", xlab="Number of Dimensions", line=2)
dev.off()

# ---------- Bats (Extended Data Fig. 5c)  ----------

# Load phylobetadiversity turnover distance matrix
pb_dist <- readRDS("Data/bat_pb_dist.RDS")

# Calculate stress values for NDMS with 1-6 axes
Stress <- c()
for(i in 1:6){
  set.seed(123)
  func_mds <- vegan::metaMDS(pb_dist, k=i)
  Stress[i] <- func_mds$stress
}

# Save scree plot
svg("Results/pb_bat_scree.svg",
    width = 2.5, height = 1.5, pointsize = 7)
par(mar=c(3,3,1,.5))
plot(1:6, Stress, axes=FALSE, xlab="", ylab="", ylim=c(0,.4))
abline(h=.2, lty="longdash", col="red")
lines(1:6, Stress, lwd=.5)
axis(1, lwd=.5)
axis(2, lwd=.5)
title(ylab="Stress", xlab="Number of Dimensions", line=2)
dev.off()
