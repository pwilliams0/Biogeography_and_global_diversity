#===========================================================#
# Calculate phylobetadiversity for birds, cell groups 1 & 3 #
#===========================================================#

library(tidyverse)
library(betapart)
library(picante)
library(ape)

# Split cells into 5 groups, every fifth cell
cell_nums <- unique(read.csv(file="Data/bird_cell_df.csv",
                             stringsAsFactors=FALSE)$cell_id)
group_1 <- cell_nums[c(TRUE, FALSE, FALSE, FALSE, FALSE)]
group_2 <- cell_nums[c(FALSE, TRUE, FALSE, FALSE, FALSE)]
group_3 <- cell_nums[c(FALSE, FALSE, TRUE, FALSE, FALSE)]
group_4 <- cell_nums[c(FALSE, FALSE, FALSE, TRUE, FALSE)]
group_5 <- cell_nums[c(FALSE, FALSE, FALSE, FALSE, TRUE)]

# Select set of cells to calculate phylobetadiversity
cell_list <- c(group_1, group_3)

# Create community matrix
bird_comm_mat <- read.csv(file="Data/bird_cell_df.csv",
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
bird_comm_mat <- bird_comm_mat[,order(colnames(bird_comm_mat))]

# Load phylogeny, prune for species in community matrix
bird_phylo <- picante::prune.sample(
  bird_comm_mat,
  ape::read.tree("Data/Raw/Stage2_MayrParSho_Hackett_mean phylogeny_ultrametric.tre"))

print("Done loading data")
# Run phylogenetic beta diversity (Sorenson)
start_time <- Sys.time()
bird_pb <- betapart::phylo.beta.pair(
  bird_comm_mat,
  bird_phylo,
  index.family="sorensen")
end_time <- Sys.time()
print(end_time - start_time)

saveRDS(bird_pb, "Data/bird_pb_13.rds")
# sor is total phylobetdiveristy (turnover + nestedness)
# sim is phylobetdiveristy turnover
# sne is phylobetdiveristy nestedness
write.table(as.matrix(bird_pb$phylo.beta.sim),
            file="Data/bird_pb_sim_13.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)
#write.table(as.matrix(bird_pb$phylo.beta.sor),
#            file="Data/bird_pb_sor_13.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)
#write.table(as.matrix(bird_pb$phylo.beta.sne),
#            file="Data/bird_pb_sne_13.txt", col.names=TRUE, row.name=TRUE, sep="\t", quote=FALSE)

print("Done with bird_13")

# [Runtime 12 hr 3 min]
# [255 GB memory used]
