#==========================================================#
# Calculate phylobetadiversity for bats, cell groups 2 & 3 #
#==========================================================#

library(tidyverse)
library(betapart)
library(picante)
library(ape)

# Split cells into 3 groups, every third cell
cell_nums <- unique(read.csv(file="Data/bat_cell_df.csv",
                             stringsAsFactors=FALSE)$cell_id)
group_1 <- cell_nums[c(TRUE, FALSE, FALSE)]
group_2 <- cell_nums[c(FALSE, TRUE, FALSE)]
group_3 <- cell_nums[c(FALSE, FALSE, TRUE)]

# Select set of cells to calculate phylobetadiversity
cell_list <- c(group_2, group_3)

# Create community matrix
bat_comm_mat <- read.csv(file="Data/bat_cell_df.csv",
                          stringsAsFactors=FALSE) %>%
  dplyr::filter(cell_id %in% cell_list) %>%
  group_by(cell_id, names_phylo) %>%
  slice(1) %>% ungroup() %>%
  mutate(present = 1) %>%
  dplyr::select(cell_id, names_phylo, present) %>%
  pivot_wider(names_from = names_phylo,
              values_from = present) %>%
  replace(.,is.na(.),0) %>%
  mutate(cell_id = paste("X",cell_id,sep="")) %>%
  column_to_rownames("cell_id") %>%
  as.matrix()
bat_comm_mat <- bat_comm_mat[,order(colnames(bat_comm_mat))]

# Load phylogeny, prune for species in community matrix
bat_phylo <- picante::prune.sample(
  bat_comm_mat,
  ape::read.tree("Data/Raw/Upham_mean_phylogeny_ultrametric.tree"))

print("Done loading data")
# Run phylogenetic beta diversity (Sorenson)
start_time <- Sys.time()
bat_pb <- betapart::phylo.beta.pair(
  bat_comm_mat,
  bat_phylo,
  index.family="sorensen")
end_time <- Sys.time()
print(end_time - start_time)

saveRDS(bat_pb, "Data/bat_pb_23.rds")
# sor is total phylobetdiveristy (turnover + nestedness)
# sim is phylobetdiveristy turnover
# sne is phylobetdiveristy nestedness
write.table(as.matrix(bat_pb$phylo.beta.sim),
            file="Data/bat_pb_sim_23.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)
#write.table(as.matrix(bat_pb$phylo.beta.sor),
#            file="Data/bat_pb_sor_23.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)
#write.table(as.matrix(bat_pb$phylo.beta.sne),
#            file="Data/bat_pb_sne_23.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)

print("Done with bat_23")
