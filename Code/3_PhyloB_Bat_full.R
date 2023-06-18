#==========================================================#
# Combine files to create full phylobetadiversity turnover #
#         distance matrix for bats, save NMDS axes         #
#==========================================================#

# ----- FULL PHYLOBETADIVERSITY DISTANCE MATRIX ----------

library(progress)

# Load cell id numbers
cell_nums <- unique(read.csv(file="Data/bat_cell_df.csv",
                             stringsAsFactors=FALSE)$cell_id)

# Load phylobetadiversity files
bat_pb_sim_12 <- read.table("Data/bat_pb_sim_12.txt")
bat_pb_sim_13 <- read.table("Data/bat_pb_sim_13.txt")
bat_pb_sim_23 <- read.table("Data/bat_pb_sim_23.txt")

# Build phylobetadiversity dataframe for all cells
start <- Sys.time()
bat_pb_all <- data.frame(row.names = cell_nums, stringsAsFactors = FALSE)
pb <- progress_bar$new(format=" running [:bar] :percent in :elapsed eta: :eta",
                       total = length(cell_nums), clear = FALSE, width = 100)
for(i in 1:length(cell_nums)){
  for(j in 1:length(cell_nums)){
    # Find the pairwise comparisons of cells i and j
    cell_1 <- paste("X",cell_nums[i],sep="")
    cell_2 <- paste("X",cell_nums[j],sep="")
    data_12 <- bat_pb_sim_12[cell_1,cell_2]
    data_13 <- bat_pb_sim_13[cell_1,cell_2]
    data_23 <- bat_pb_sim_23[cell_1,cell_2]
    # Take average value
    # 2/3 of pairwise comparisons will have a single value
    # 1/3 of pairwise comparisons will have two value, should be the same
    bat_pb_all[i,j] <- mean(na.omit(c(data_12,
                                      data_13,
                                      data_23)))
  }
  pb$tick()
}
colnames(bat_pb_all) <- cell_nums
row.names(bat_pb_all) <- cell_nums
Sys.time() - start # [1.9 hours]

# Save as a new distance matrix
write.table(bat_pb_all, "Data/bat_pb_all.txt")
bat_pb_dist <- as.dist(bat_pb_all)
saveRDS(bat_pb_dist, "Data/bat_pb_dist.RDS")

# ----- NMDS axes ----------

library(tidyverse)

# Load phylobetadiversity distance matrix
bat_pb_dist <- readRDS("Data/bat_pb_dist.RDS")

# Run NMDS, 3 axes
set.seed(123)
start <- Sys.time()
phylo_mds <- vegan::metaMDS(bat_pb_dist, k=3)
Sys.time() - start
phylo_mds$stress
cell_nums <- unique(read.csv(file="Data/bat_cell_df.csv",
                             stringsAsFactors=FALSE)$cell_id)
bat_pb_NMDS <- as.data.frame(phylo_mds$points) %>%
  rownames_to_column("cell_id")

# Save NMDS axes
write.csv(bat_pb_NMDS, "Data/bat_pb_NMDS.csv")
