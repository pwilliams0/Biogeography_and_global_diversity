#=========================================================#
# Calculate species richness (SR), phylogenetic diversity #
#      (PD) and functional richness (FR) for mammals      #
#=========================================================#

library(tidyverse)

# ----- SPECIES RICHNESS ----------

mamm_SR_cells <- read.csv("Data/mamm_cell_df.csv") %>%
  group_by(cell_id) %>%
  mutate(SR = n()) %>%
  slice(1) %>%
  dplyr::select(cell_id, SR)

# Save file
write.csv(mamm_SR_cells, "Data/mamm_SR_cells.csv")

# ----- PHYLOGENETIC DIVERSITY ----------
library(ape)
library(picante)

# Create community matrix
mamm_comm_mat <- read.csv(file="Data/mamm_cell_df.csv",
                          stringsAsFactors=FALSE) %>%
  mutate(present = 1) %>%
  dplyr::select(cell_id, names_IUCN, present) %>%
  pivot_wider(names_from = names_IUCN,
              values_from = present) %>%
  replace(.,is.na(.),0) %>%
  mutate(cell_id = paste("X",cell_id,sep="")) %>%
  column_to_rownames("cell_id") %>%
  as.matrix()
mamm_comm_mat <- mamm_comm_mat[,order(colnames(mamm_comm_mat))]

# Load phylogeny, prune for species in community matrix
# Use ultrametric consensus phylogeny derived from 1000 phylogenies
mamm_phylo <- picante::prune.sample(
  mamm_comm_mat,
  ape::read.tree("Data/Raw/Upham_mean_phylogeny_ultrametric.tree"))

# Calculate PD
start <- Sys.time()
PD_mamm <- picante::pd(mamm_comm_mat,
                       mamm_phylo,
                       include.root = TRUE)
Sys.time() - start

# Save file
mamm_PD_cells <- PD_mamm %>%
  dplyr::select(PD) %>%
  rownames_to_column("cell_id") %>%
  mutate(cell_id = substring(cell_id, 2))
write.csv(mamm_PD_cells, "Data/mamm_PD_cells.csv")

# ----- FUNCTIONAL RICHNESS ----------
library(mFD)

# ---------- Make PCoA of traits ----------

# Load trait categories
traits_cat <- read.csv("Data/Raw/mamm_traits_cat.csv")

# Load trait data, one row per species
spp_traits <- read.csv("Data/mamm_cell_df.csv") %>%
  group_by(names_IUCN) %>% slice(1) %>% ungroup() %>%
  dplyr::select(names_IUCN, Diet.Inv, Diet.Vert,
                Diet.Fish, Diet.Scav, Diet.Fruit,
                Diet.Nect, Diet.Seed, Diet.Herb,
                ForStrat, LogBodyMass) %>%
  mutate(ForStrat = as.factor(ForStrat)) %>%
  column_to_rownames("names_IUCN")

# Create distance matrix of species using functional traits
start <- Sys.time()
sp_dist_mamm <- mFD::funct.dist(
  sp_tr         = spp_traits,
  tr_cat        = traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center", # Standardize variables
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)
Sys.time() - start # [6.0 minutes]

# Create PCoA using distance matrix
start <- Sys.time()
fspaces_quality_mamm <- mFD::quality.fspaces(
  sp_dist             = sp_dist_mamm,
  maxdim_pcoa         = 5,
  deviation_weighting = c("absolute", "squarred"),
  fdendro             = "average")
Sys.time() - start # [5.8 minutes]

# Save PCoA axes coordinates
sp_faxes_coord_mamm <- fspaces_quality_mamm$"details_fspaces"$"sp_pc_coord"

# Save files
saveRDS(fspaces_quality_mamm, "Data/mamm_func_qual.RDS")
saveRDS(sp_faxes_coord_mamm, "Data/mamm_func_PCoA.RDS")

# Check quality of PCoA
fspaces_quality_mamm$"quality_fspaces"
# Lowest mad score is 5 axes, but that gets to be very computationally demanding
# Use 4 axes going forward

# ---------- Calculate functional richness ----------

mamm_func_PCoA <- readRDS("Data/mamm_func_PCoA.RDS")

# Create community matrix
mamm_comm_mat <- read.csv(file="Data/mamm_cell_df.csv",
                          stringsAsFactors=FALSE) %>%
  mutate(present = 1) %>%
  # Remove cells with <= 4 species, need more species than axes
  group_by(cell_id) %>%
  mutate(SR = n()) %>%
  ungroup() %>%
  dplyr::filter(SR > 4) %>%
  dplyr::select(cell_id, names_IUCN, present) %>%
  pivot_wider(names_from = names_IUCN,
              values_from = present) %>%
  replace(.,is.na(.),0) %>%
  mutate(cell_id = paste("X",cell_id,sep="")) %>%
  column_to_rownames("cell_id") %>%
  as.matrix()
mamm_comm_mat <- mamm_comm_mat[,order(colnames(mamm_comm_mat))]

# Calculate functional richness
start <- Sys.time()
mamm_FD_results <- mFD::alpha.fd.multidim(
  sp_faxes_coord = mamm_func_PCoA[,c("PC1","PC2",
                                     "PC3","PC4")],
  asb_sp_w = mamm_comm_mat,
  ind_vect = c("fric"),
  scaling = TRUE, check_input = TRUE, verbose = TRUE,
  details_returned = TRUE)$functional_diversity_indices
Sys.time() - start # [1.2 hours]

# Save file
mamm_FR_cells <- mamm_FD_results %>%
  mutate(FR = fric) %>%
  dplyr::select(FR) %>%
  rownames_to_column("cell_id") %>%
  mutate(cell_id = substring(cell_id, 2))
write.csv(mamm_FR_cells, "Data/mamm_FR_cells.csv")
