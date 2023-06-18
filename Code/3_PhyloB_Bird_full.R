#==========================================================#
# Combine files to create full phylobetadiversity turnover #
#        distance matrix for birds, save NMDS axes         #
#==========================================================#

# ----- FULL PHYLOBETADIVERSITY DISTANCE MATRIX ----------

library(progress)

# Load cell id numbers
cell_nums <- unique(read.csv(file="Data/bird_cell_df.csv",
                             stringsAsFactors=FALSE)$cell_id)

# Load phylobetadiversity files
bird_pb_sim_12 <- read.table("Data/bird_pb_sim_12.txt")
bird_pb_sim_13 <- read.table("Data/bird_pb_sim_13.txt")
bird_pb_sim_14 <- read.table("Data/bird_pb_sim_14.txt")
bird_pb_sim_15 <- read.table("Data/bird_pb_sim_15.txt")
bird_pb_sim_23 <- read.table("Data/bird_pb_sim_23.txt")
bird_pb_sim_24 <- read.table("Data/bird_pb_sim_24.txt")
bird_pb_sim_25 <- read.table("Data/bird_pb_sim_25.txt")
bird_pb_sim_34 <- read.table("Data/bird_pb_sim_34.txt")
bird_pb_sim_35 <- read.table("Data/bird_pb_sim_35.txt")
bird_pb_sim_45 <- read.table("Data/bird_pb_sim_45.txt")

# Build phylobetadiversity dataframe for all cells
start <- Sys.time()
bird_pb_all <- data.frame(row.names = cell_nums, stringsAsFactors = FALSE)
pb <- progress_bar$new(format=" running [:bar] :percent in :elapsed eta: :eta",
                       total = length(cell_nums), clear = FALSE, width = 80)
for(i in 1:length(cell_nums)){
  for(j in 1:length(cell_nums)){
    # Find the pairwise comparisons of cells i and j
    cell_1 <- paste("X",cell_nums[i],sep="")
    cell_2 <- paste("X",cell_nums[j],sep="")
    data_12 <- bird_pb_sim_12[cell_1,cell_2]
    data_13 <- bird_pb_sim_13[cell_1,cell_2]
    data_14 <- bird_pb_sim_14[cell_1,cell_2]
    data_15 <- bird_pb_sim_15[cell_1,cell_2]
    data_23 <- bird_pb_sim_23[cell_1,cell_2]
    data_24 <- bird_pb_sim_24[cell_1,cell_2]
    data_25 <- bird_pb_sim_25[cell_1,cell_2]
    data_34 <- bird_pb_sim_34[cell_1,cell_2]
    data_35 <- bird_pb_sim_35[cell_1,cell_2]
    data_45 <- bird_pb_sim_45[cell_1,cell_2]
    # Take average value
    # 4/5 of pairwise comparisons will have a single value
    # 1/5 of pairwise comparisons will have four values, should be the same
    bird_pb_all[i,j] <- mean(na.omit(c(data_12,
                                       data_13,
                                       data_14,
                                       data_15,
                                       data_23,
                                       data_24,
                                       data_25,
                                       data_34,
                                       data_35,
                                       data_45)))
  }
  pb$tick()
}
colnames(bird_pb_all) <- cell_nums
row.names(bird_pb_all) <- cell_nums
Sys.time() - start # [1.7 hours]

# Save as a new distance matrix
write.table(bird_pb_all, "Data/bird_pb_all.txt")
bird_pb_dist <- as.dist(bird_pb_all)
saveRDS(bird_pb_dist, "Data/bird_pb_dist.RDS")

# ----- NMDS axes ----------

library(tidyverse)

# Load phylobetadiversity distance matrix
bird_pb_dist <- readRDS("Data/bird_pb_dist.RDS")

# Run NMDS, 3 axes
set.seed(123)
start <- Sys.time()
phylo_mds <- vegan::metaMDS(bird_pb_dist, k=3)
Sys.time() - start
phylo_mds$stress
cell_nums <- unique(read.csv(file="Data/bird_cell_df.csv",
                             stringsAsFactors=FALSE)$cell_id)
bird_pb_NMDS <- as.data.frame(phylo_mds$points) %>%
  rownames_to_column("cell_id")

# Save NMDS axes
write.csv(bird_pb_NMDS, "Data/bird_pb_NMDS.csv")
