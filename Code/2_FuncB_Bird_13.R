#==================================================================#
# Calculate functional beta diversity for birds, cell groups 1 & 3 #
#==================================================================#

library(tidyverse)
library(mFD)

# Split cells into 4 groups, every fourth cell
cell_nums <- unique(read.csv(file="Data/bird_cell_df.csv",
                             stringsAsFactors=FALSE)$cell_id)
group_1 <- cell_nums[c(TRUE, FALSE, FALSE, FALSE)]
group_2 <- cell_nums[c(FALSE, TRUE, FALSE, FALSE)]
group_3 <- cell_nums[c(FALSE, FALSE, TRUE, FALSE)]
group_4 <- cell_nums[c(FALSE, FALSE, FALSE, TRUE)]

# Select set of cells to calculate functional beta diversity
cell_list <- sort(c(group_1, group_3))

# Load functional traits PCoA
bird_func_PCoA <- readRDS("Data/bird_func_PCoA.RDS")

# Create community matrix
bird_comm_mat <- read.csv(file="Data/bird_cell_df.csv",
                          stringsAsFactors=FALSE) %>%
  mutate(present = 1) %>%
  # Remove cells with <= 4 species
  group_by(cell_id) %>%
  mutate(SR = n()) %>%
  ungroup() %>%
  dplyr::filter(SR > 4,
                cell_id %in% cell_list) %>%
  dplyr::select(cell_id, names_IUCN, present) %>%
  pivot_wider(names_from = names_IUCN,
              values_from = present) %>%
  replace(.,is.na(.),0) %>%
  mutate(cell_id = paste("X",cell_id,sep="")) %>%
  column_to_rownames("cell_id") %>%
  as.matrix()
bird_comm_mat <- bird_comm_mat[,order(colnames(bird_comm_mat))]

print("Done loading data")

# Calculate functional beta diversity
start <- Sys.time()
bird_fb <- mFD::beta.fd.multidim(
  sp_faxes_coord   = bird_func_PCoA[ , c("PC1", "PC2",
                                         "PC3", "PC4")],
  asb_sp_occ       = bird_comm_mat,
  check_input      = TRUE,
  beta_family      = c("Sorensen"),
  details_returned = TRUE)
end <- Sys.time()
time_diff <- end - start

print("Done calculating functional beta diversity 13")
print(time_diff)

saveRDS(bird_fb, "Data/bird_fb_13.RDS")
#write.table(as.matrix(bird_fb$pairasb_fbd_indices$sor_diss),
#            file="Data/bird_fb_sor_13.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)
write.table(as.matrix(bird_fb$pairasb_fbd_indices$sor_turn),
            file="Data/bird_fb_turn_13.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)
#write.table(as.matrix(bird_fb$pairasb_fbd_indices$sor_nest),
#            file="Data/bird_fb_nest_13.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)

print("Done with bird_13")

# [Runtime 57 hrs]
# [62.2 GB memory used]
