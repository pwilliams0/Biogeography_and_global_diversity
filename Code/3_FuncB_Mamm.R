#==========================================================#
# Calculate functional beta diversity turnover for mammals #
#==========================================================#

library(tidyverse)
library(mFD)

# Load functional traits PCoA (created in 1_SR_PD_FR_Mamm.R)
mamm_func_PCoA <- readRDS("Data/mamm_func_PCoA.RDS")

# Create community matrix
mamm_comm_mat <- read.csv(file="Data/mamm_cell_df.csv",
                          stringsAsFactors=FALSE) %>%
  mutate(present = 1) %>%
  # Remove cells with <= 4 species
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

print("Done loading data")

# Calculate functional beta diversity
start <- Sys.time()
mamm_fb <- mFD::beta.fd.multidim(
  sp_faxes_coord   = mamm_func_PCoA[ , c("PC1", "PC2",
                                         "PC3", "PC4")],
  asb_sp_occ       = mamm_comm_mat,
  check_input      = TRUE,
  beta_family      = c("Sorensen"),
  betapart_para = TRUE,
  details_returned = TRUE)
end <- Sys.time()
time_diff <- end - start

print("Done calculating functional beta diversity")
print(time_diff) # [2.05 days, 69.4 GB]

# Save distance matrices
saveRDS(mamm_fb, "Data/mamm_fb.RDS")
#write.table(as.matrix(mamm_fb$pairasb_fbd_indices$sor_diss),
#            file="Data/mamm_fb_sor.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)
write.table(as.matrix(mamm_fb$pairasb_fbd_indices$sor_turn),
            file="Data/mamm_fb_turn.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)
#write.table(as.matrix(mamm_fb$pairasb_fbd_indices$sor_nest),
#            file="Data/mamm_fb_nest.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)

# ----- MEAN FUNCTIONAL BETA DIVERSITY TURNOVER ----------

library(usedist)

# Make list of cells with more species than functional axes
cell_list <- read.csv("Data/mamm_SR_cells.csv") %>%
  filter(SR > 4) %>%
  mutate(cell_id = paste("X", cell_id, sep="")) 
# Load functional beta diversity turnover distance matrix
# Select only cells with more species than functional axes
fb_dist <- dist_subset(
  as.dist(read.table("Data/mamm_fb_turn.txt")),
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
write.csv(fb_df,"Data/mamm_mean_fb.csv")
