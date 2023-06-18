#==========================================================#
# Combine files to create full phylobetadiversity turnover #
#       distance matrix for mammals, save NMDS axes        #
#==========================================================#

# ----- FULL PHYLOBETADIVERSITY DISTANCE MATRIX ----------

library(progress)

# Load cell id numbers
cell_nums <- unique(read.csv(file="Data/mamm_cell_df.csv",
                             stringsAsFactors=FALSE)$cell_id)

# Load phylobetadiversity files
mamm_pb_sim_12 <- read.table("Data/mamm_pb_sim_12.txt")
mamm_pb_sim_13 <- read.table("Data/mamm_pb_sim_13.txt")
mamm_pb_sim_14 <- read.table("Data/mamm_pb_sim_14.txt")
mamm_pb_sim_23 <- read.table("Data/mamm_pb_sim_23.txt")
mamm_pb_sim_24 <- read.table("Data/mamm_pb_sim_24.txt")
mamm_pb_sim_34 <- read.table("Data/mamm_pb_sim_34.txt")

# Build phylobetadiversity dataframe for all cells
start <- Sys.time()
mamm_pb_all <- data.frame(row.names = cell_nums, stringsAsFactors = FALSE)
pb <- progress_bar$new(format=" running [:bar] :percent in :elapsed eta: :eta",
                       total = length(cell_nums), clear = FALSE, width = 80)
for(i in 1:length(cell_nums)){
  for(j in 1:length(cell_nums)){
    # Find the pairwise comparisons of cells i and j
    cell_1 <- paste("X",cell_nums[i],sep="")
    cell_2 <- paste("X",cell_nums[j],sep="")
    data_12 <- mamm_pb_sim_12[cell_1,cell_2]
    data_13 <- mamm_pb_sim_13[cell_1,cell_2]
    data_14 <- mamm_pb_sim_14[cell_1,cell_2]
    data_23 <- mamm_pb_sim_23[cell_1,cell_2]
    data_24 <- mamm_pb_sim_24[cell_1,cell_2]
    data_34 <- mamm_pb_sim_34[cell_1,cell_2]
    # Take average value
    # 3/4 of pairwise comparisons will have a single value
    # 1/4 of pairwise comparisons will have three values, should be the same
    mamm_pb_all[i,j] <- mean(na.omit(c(data_12,
                                      data_13,
                                      data_14,
                                      data_23,
                                      data_24,
                                      data_34)))
  }
  pb$tick()
}
colnames(mamm_pb_all) <- cell_nums
row.names(mamm_pb_all) <- cell_nums
Sys.time() - start # [1.9 hours]

# Save as a new distance matrix
write.table(mamm_pb_all, "Data/mamm_pb_all.txt")
mamm_pb_dist <- as.dist(mamm_pb_all)
saveRDS(mamm_pb_dist, "Data/mamm_pb_dist.RDS")

# ----- NMDS axes ----------

library(tidyverse)

# Load phylobetadiversity distance matrix
mamm_pb_dist <- readRDS("Data/mamm_pb_dist.RDS")

# Run NMDS, 3 axes
set.seed(123)
start <- Sys.time()
phylo_mds <- vegan::metaMDS(mamm_pb_dist, k=3)
Sys.time() - start
phylo_mds$stress
cell_nums <- unique(read.csv(file="Data/mamm_cell_df.csv",
                             stringsAsFactors=FALSE)$cell_id)
mamm_pb_NMDS <- as.data.frame(phylo_mds$points) %>%
  rownames_to_column("cell_id")

# Save NMDS axes
write.csv(mamm_pb_NMDS, "Data/mamm_pb_NMDS.csv")
