#=============================================================#
# Calculate phylobetadiversity for mammals, cell groups 1 & 3 #
#=============================================================#

library(tidyverse)
library(betapart)
library(picante)
library(ape)

# Split cells into 4 groups, every fourth cell
cell_nums <- unique(read.csv(file="Data/mamm_cell_df.csv",
                             stringsAsFactors=FALSE)$cell_id)
group_1 <- cell_nums[c(TRUE, FALSE, FALSE, FALSE)]
group_2 <- cell_nums[c(FALSE, TRUE, FALSE, FALSE)]
group_3 <- cell_nums[c(FALSE, FALSE, TRUE, FALSE)]
group_4 <- cell_nums[c(FALSE, FALSE, FALSE, TRUE)]

# Select set of cells to calculate phylobetadiversity
cell_list <- c(group_1, group_3)

# Create community matrix
mamm_comm_mat <- read.csv(file="Data/mamm_cell_df.csv",
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
mamm_comm_mat <- mamm_comm_mat[,order(colnames(mamm_comm_mat))]

# Load phylogeny, prune for species in community matrix
mamm_phylo <- picante::prune.sample(
  mamm_comm_mat,
  ape::read.tree("Data/Raw/Upham_mean_phylogeny_ultrametric.tree"))

print("Done loading data")
# Run phylogenetic beta diversity (Sorenson)
start_time <- Sys.time()
mamm_pb <- betapart::phylo.beta.pair(
  mamm_comm_mat,
  mamm_phylo,
  index.family="sorensen")
end_time <- Sys.time()
print(end_time - start_time)

saveRDS(mamm_pb, "Data/mamm_pb_13.rds")
# sor is total phylobetdiveristy (turnover + nestedness)
# sim is phylobetdiveristy turnover
# sne is phylobetdiveristy nestedness
write.table(as.matrix(mamm_pb$phylo.beta.sim),
            file="Data/mamm_pb_sim_13.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)
#write.table(as.matrix(mamm_pb$phylo.beta.sor),
#            file="Data/mamm_pb_sor_13.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)
#write.table(as.matrix(mamm_pb$phylo.beta.sne),
#            file="Data/mamm_pb_sne_13.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)

print("Done with mamm_13")

# [Runtime 5 hr 57 min]
# [229 GB memory used]
