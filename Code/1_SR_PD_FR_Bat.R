#=========================================================#
# Calculate species richness (SR), phylogenetic diversity #
#       (PD) and functional richness (FR) for bats        #
#=========================================================#

library(tidyverse)

# ----- SPECIES RICHNESS ----------

bat_SR_cells <- read.csv("Data/bat_cell_df.csv") %>%
  group_by(cell_id) %>%
  mutate(SR = n()) %>%
  slice(1) %>%
  dplyr::select(cell_id, SR)

# Save file
write.csv(bat_SR_cells, "Data/bat_SR_cells.csv")

# ----- PHYLOGENETIC DIVERSITY ----------
library(ape)
library(picante)

# Create community matrix
bat_comm_mat <- read.csv(file="Data/bat_cell_df.csv",
                          stringsAsFactors=FALSE) %>%
  mutate(present = 1) %>%
  dplyr::select(cell_id, names_IUCN, present) %>%
  pivot_wider(names_from = names_IUCN,
              values_from = present) %>%
  replace(.,is.na(.),0) %>%
  mutate(cell_id = paste("X",cell_id,sep="")) %>%
  column_to_rownames("cell_id") %>%
  as.matrix()
bat_comm_mat <- bat_comm_mat[,order(colnames(bat_comm_mat))]

# Load phylogeny, prune for species in community matrix
# Use ultrametric consensus phylogeny derived from 1000 phylogenies
bat_phylo <- picante::prune.sample(
  bat_comm_mat,
  ape::read.tree("Data/Raw/Upham_mean_phylogeny_ultrametric.tree"))

# Calculate PD
start <- Sys.time()
PD_bat <- picante::pd(bat_comm_mat,
                       bat_phylo,
                       include.root = TRUE)
Sys.time() - start

# Save file
bat_PD_cells <- PD_bat %>%
  dplyr::select(PD) %>%
  rownames_to_column("cell_id") %>%
  mutate(cell_id = substring(cell_id, 2))
write.csv(bat_PD_cells, "Data/bat_PD_cells.csv")

# ----- FUNCTIONAL RICHNESS ----------
library(mFD)

# ---------- Make PCoA of traits ----------

# Load trait categories
traits_cat <- read.csv("Data/Raw/mamm_traits_cat.csv")

# Load trait data, one row per species
spp_traits <- read.csv("Data/bat_cell_df.csv") %>%
  group_by(names_IUCN) %>% slice(1) %>% ungroup() %>%
  dplyr::select(names_IUCN, Diet.Inv, Diet.Vert,
                Diet.Fish, Diet.Scav, Diet.Fruit,
                Diet.Nect, Diet.Seed, Diet.Herb,
                ForStrat, LogBodyMass) %>%
  mutate(ForStrat = as.factor(ForStrat)) %>%
  column_to_rownames("names_IUCN")

# Create distance matrix of species using functional traits
start <- Sys.time()
sp_dist_bat <- mFD::funct.dist(
  sp_tr         = spp_traits,
  tr_cat        = traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center", # Standardize variables
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)
Sys.time() - start # [14 seconds]

# Create PCoA using distance matrix
start <- Sys.time()
fspaces_quality_bat <- mFD::quality.fspaces(
  sp_dist             = sp_dist_bat,
  maxdim_pcoa         = 5,
  deviation_weighting = c("absolute", "squarred"),
  fdendro             = "average")
Sys.time() - start # [5 seconds]

# Save PCoA axes coordinates
sp_faxes_coord_bat <- fspaces_quality_bat$"details_fspaces"$"sp_pc_coord"

# Save files
saveRDS(fspaces_quality_bat, "Data/bat_func_qual.RDS")
saveRDS(sp_faxes_coord_bat, "Data/bat_func_PCoA.RDS")

# Check quality of PCoA
fspaces_quality_bat$"quality_fspaces"
# Lowest mad score is 2 axes, but that's far fewer axes used for birds/mammals 
# Use 3 axes going forward

# ---------- Calculate functional richness ----------

bat_func_PCoA <- readRDS("Data/bat_func_PCoA.RDS")

# Create community matrix
bat_comm_mat <- read.csv(file="Data/bat_cell_df.csv",
                          stringsAsFactors=FALSE) %>%
  mutate(present = 1) %>%
  # Remove cells with <= 3 species, need more species than axes
  group_by(cell_id) %>%
  mutate(SR = n()) %>%
  ungroup() %>%
  dplyr::filter(SR > 3) %>%
  dplyr::select(cell_id, names_IUCN, present) %>%
  pivot_wider(names_from = names_IUCN,
              values_from = present) %>%
  replace(.,is.na(.),0) %>%
  mutate(cell_id = paste("X",cell_id,sep="")) %>%
  column_to_rownames("cell_id") %>%
  as.matrix()
bat_comm_mat <- bat_comm_mat[,order(colnames(bat_comm_mat))]

# Calculate functional richness
start <- Sys.time()
bat_FD_results <- mFD::alpha.fd.multidim(
  sp_faxes_coord = bat_func_PCoA[,c("PC1","PC2","PC3")],
  asb_sp_w = bat_comm_mat,
  ind_vect = c("fric"),
  scaling = TRUE, check_input = TRUE, verbose = TRUE,
  details_returned = TRUE)$functional_diversity_indices
Sys.time() - start # [12.4 minutes]

# Save file
bat_FR_cells <- bat_FD_results %>%
  mutate(FR = fric) %>%
  dplyr::select(FR) %>%
  rownames_to_column("cell_id") %>%
  mutate(cell_id = substring(cell_id, 2))
write.csv(bat_FR_cells, "Data/bat_FR_cells.csv")
