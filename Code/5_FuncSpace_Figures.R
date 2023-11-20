#========================================#
# Figures of overlap in functional space #
#========================================#

library(tidyverse)
library(mFD)

# ----- MAMMALS, AUSTRALIA VS SOUTH AFRICA (Fig. 5) ----------

# Load PCoA
mamm_func_PCoA <- readRDS("Data/mamm_func_PCoA.RDS")
# Load functional beta diversity
mamm_fb <- readRDS("Data/mamm_fb.RDS")
# Load trait data
mamm_traits <- read.csv("Data/mamm_cell_df.csv") %>%
  group_by(names_IUCN) %>%
  slice(1) %>%
  ungroup() %>%
  column_to_rownames("names_IUCN")

# View correlation of traits and axes to interpret functional axes
mamm_axes <- mFD::traits.faxes.cor(
  sp_tr          = mamm_traits[,6:10], 
  sp_faxes_coord = mamm_func_PCoA[ , c("PC1", "PC4")], 
  plot           = TRUE)
View(mamm_axes$tr_faxes_stat)
mamm_axes$tr_faxes_plot
mamm_axes <- mFD::traits.faxes.cor(
  sp_tr          = mamm_traits[,11:15],
  sp_faxes_coord = mamm_func_PCoA[ , c("PC1", "PC4")], 
  plot           = TRUE)
View(mamm_axes$tr_faxes_stat)
mamm_axes$tr_faxes_plot
mamm_axes <- mFD::traits.faxes.cor(
  sp_tr          = mamm_traits[,6:10], 
  sp_faxes_coord = mamm_func_PCoA[ , c("PC2", "PC3")], 
  plot           = TRUE)
View(mamm_axes$tr_faxes_stat)
mamm_axes$tr_faxes_plot
mamm_axes <- mFD::traits.faxes.cor(
  sp_tr          = mamm_traits[,11:15],
  sp_faxes_coord = mamm_func_PCoA[ , c("PC2", "PC3")], 
  plot           = TRUE)
View(mamm_axes$tr_faxes_stat)
mamm_axes$tr_faxes_plot

# Plot overlap in functional space
# Cell 68 is in New South Wales, Australia
# Cell 83 is in Easter Cape, South Africa
# Plot particular species to aid interpretation, but don't save on final plot
beta_plot_fruits <- mFD::beta.multidim.plot(
  output_beta_fd_multidim = mamm_fb,
  plot_asb_nm             = c("X68", "X83"),
  beta_family             = c("Sorenson"),
  #plot_sp_nm              = c("Loxodonta_africana",
  #                            "Panthera_pardus",
  #                            "Trichosurus_vulpecula",
  #                            "Phascolarctos_cinereus",
  #                            "Graphiurus_murinus"),
  faxes                   = paste0("PC", 1:4),
  range_faxes             = c(NA, NA),
  color_bg                = "#E2E2E2",
  shape_sp                = c("pool" = 3.0, asb1 = 16, asb2 = 16),
  size_sp                 = c("pool" = 0, asb1 =  1, asb2 =  1),
  color_sp                = c("pool" = NA, asb1 = "#023FA5", asb2 = "#8E063B"),
  fill_sp                 = c("pool" = NA, asb1 = "#023FA5", asb2 = "#8E063B"),
  fill_vert               = c("pool" = NA, asb1 = "#023FA5", asb2 = "#8E063B"),
  color_ch                = c("pool" = NA, asb1 = "#023FA5", asb2 = "#8E063B"),
  fill_ch                 = c("pool" = "white", asb1 = "#023FA5", asb2 = "#AB5468"),
  alpha_ch                = c("pool" = 1, asb1 = 0.2, asb2 = 0.3),
  nm_size                 = 3,
  nm_color                = "black",
  nm_fontface             = "plain",
  check_input             = TRUE)

# Save figure, axes 1 and 4
svg("Results/mamm_FB_overlap1.svg",
    width = 4, height = 4, pointsize = 12)
beta_plot_fruits$PC1_PC4
dev.off()

# Save figure, axes 2 and 3
svg("Results/mamm_FB_overlap2.svg",
    width = 4, height = 4, pointsize = 12)
beta_plot_fruits$PC2_PC3
dev.off()

# ----- BATS, COLOMBIA VS MALAYSIA (Extended Data Fig. 4) ----------

# Load PCoA
bat_func_PCoA <- readRDS("Data/bat_func_PCoA.RDS")
# Load functional beta diversity
bat_fb <- readRDS("Data/bat_fb.RDS")
# Load trait data
bat_traits <- read.csv("Data/bat_cell_df.csv") %>%
  group_by(names_IUCN) %>%
  slice(1) %>%
  ungroup() %>%
  column_to_rownames("names_IUCN")

# View correlation of traits and axes to interpret functional axes
bat_axes <- mFD::traits.faxes.cor(
  sp_tr          = bat_traits[,c(6:8,10:11)], 
  sp_faxes_coord = bat_func_PCoA[ , c("PC1", "PC2", "PC3")], 
  plot           = TRUE)
View(bat_axes$tr_faxes_stat)
bat_axes$tr_faxes_plot
bat_axes <- mFD::traits.faxes.cor(
  sp_tr          = bat_traits[,12:15],
  sp_faxes_coord = bat_func_PCoA[ , c("PC1", "PC2", "PC3")], 
  plot           = TRUE)
View(bat_axes$tr_faxes_stat)
bat_axes$tr_faxes_plot

# Plot overlap in functional space
# Cell 709 is in the Colombian Amazon
# Cell 859 is in Sabah, Malaysian Borneo
# Plot particular species to aid interpretation, but don't save on final plot
beta_plot_fruits <- mFD::beta.multidim.plot(
  output_beta_fd_multidim = bat_fb,
  plot_asb_nm             = c("X709", "X859"),
  beta_family             = c("Sorenson"),
  #plot_sp_nm              = c("Pteropus_vampyrus",
  #                            "Chrotopterus_auritus",
  #                            "Coelops_robinsoni",
  #                            "Desmodus_rotundus",
  #                            "Vampyrum_spectrum",
  #                            "Phyllostomus_latifolius",
  #                            "Noctilio_leporinus",
  #                            "Lichonycteris_obscura"),
  faxes                   = paste0("PC", 1:3),
  range_faxes             = c(NA, NA),
  color_bg                = "#E2E2E2",
  shape_sp                = c("pool" = 3.0, asb1 = 16, asb2 = 16),
  size_sp                 = c("pool" = 0, asb1 =  1, asb2 =  1),
  color_sp                = c("pool" = NA, asb1 = "#023FA5", asb2 = "#8E063B"),
  fill_sp                 = c("pool" = NA, asb1 = "#023FA5", asb2 = "#8E063B"),
  fill_vert               = c("pool" = NA, asb1 = "#023FA5", asb2 = "#8E063B"),
  color_ch                = c("pool" = NA, asb1 = "#023FA5", asb2 = "#8E063B"),
  fill_ch                 = c("pool" = "white", asb1 = "#023FA5", asb2 = "#AB5468"),
  alpha_ch                = c("pool" = 1, asb1 = 0.2, asb2 = 0.3),
  nm_size                 = 3,
  nm_color                = "black",
  nm_fontface             = "plain",
  check_input             = TRUE)

# Save figure, axes 1 and 2
svg("Results/bat_FB_overlap1.svg",
    width = 4, height = 4, pointsize = 12)
beta_plot_fruits$PC1_PC2
dev.off()

# Save figure, axes 1 and 3
svg("Results/bat_FB_overlap2.svg",
    width = 4, height = 4, pointsize = 12)
beta_plot_fruits$PC1_PC3
dev.off()
