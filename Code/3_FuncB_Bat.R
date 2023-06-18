#=======================================================#
# Calculate functional beta diversity turnover for bats #
#=======================================================#

library(tidyverse)
library(mFD)

# Load functional traits PCoA (created in 1_SR_PD_FR_Bat.R)
bat_func_PCoA <- readRDS("Data/bat_func_PCoA.RDS")

# Create community matrix
bat_comm_mat <- read.csv(file="Data/bat_cell_df.csv",
                          stringsAsFactors=FALSE) %>%
  mutate(present = 1) %>%
  # Remove cells with <= 3 species
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

print("Done loading data")

# Calculate functional beta diversity
start<- Sys.time()
bat_fb <- mFD::beta.fd.multidim(
  sp_faxes_coord   = bat_func_PCoA[ , c("PC1", "PC2", "PC3")],
  asb_sp_occ       = bat_comm_mat,
  check_input      = TRUE,
  beta_family      = c("Sorensen"),
  betapart_para = TRUE,
  details_returned = TRUE)
Sys.time() - start # [5.3 hours, 182.65 GB]

print("Done calculating functional beta diversity")

saveRDS(bat_fb, "Data/bat_fb.RDS")
#write.table(as.matrix(bat_fb$pairasb_fbd_indices$sor_diss),
#            file="Data/bat_fb_sor.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)
write.table(as.matrix(bat_fb$pairasb_fbd_indices$sor_turn),
            file="Data/bat_fb_turn.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)
#write.table(as.matrix(bat_fb$pairasb_fbd_indices$sor_nest),
#            file="Data/bat_fb_nest.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)

# ----- MEAN FUNCTIONAL BETA DIVERSITY TURNOVER ----------

library(usedist)

# Make list of cells with more species than functional axes
cell_list <- read.csv("Data/bat_SR_cells.csv") %>%
  filter(SR > 3) %>%
  mutate(cell_id = paste("X", cell_id, sep="")) 
# Load functional beta diversity turnover distance matrix
# Select only cells with more species than functional axes
fb_dist <- dist_subset(
  as.dist(read.table("Data/bat_fb_turn.txt")),
  cell_list$cell_id) %>%
  as.matrix(dimnames=labels(dist))
# Set diagonal to NA
diag(fb_dist) <- NA
# Calculate mean values for each cell
fb_df <- fb_dist %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id) %>%
  mutate(mean_fb = apply(fb_dist, MARGIN = 1, FUN = mean, na.rm = TRUE))
# Fix "cell_id"
fb_df$cell_id <- sub("X", "", fb_df$cell_id)

# Save file
write.csv(fb_df,"Data/bat_mean_fb.csv")
