#====================================================#
# Calculate full functional beta diversity for birds #
#====================================================#

library(progress)
library(tidyverse)

# ----- FULL TURNOVER DISTANCE MATRIX ----------

# Load cell id numbers
cell_nums <- read.csv(file="Data/bird_FR_cells.csv",
                      stringsAsFactors=FALSE)$cell_id

# Load functional beta diversity files
bird_fb_turn_12 <- read.table("Data/bird_fb_turn_12.txt")
bird_fb_turn_13 <- read.table("Data/bird_fb_turn_13.txt")
bird_fb_turn_14 <- read.table("Data/bird_fb_turn_14.txt")
bird_fb_turn_23 <- read.table("Data/bird_fb_turn_23.txt")
bird_fb_turn_24 <- read.table("Data/bird_fb_turn_24.txt")
bird_fb_turn_34 <- read.table("Data/bird_fb_turn_34.txt")

# Build functional beta diversity dataframe for all cells
start <- Sys.time()
bird_fb_all <- data.frame(row.names = cell_nums, stringsAsFactors = FALSE)
pb <- progress_bar$new(format=" running [:bar] :percent in :elapsed eta: :eta",
                       total = length(cell_nums), clear = FALSE, width = 80)
for(i in 1:length(cell_nums)){
  for(j in 1:length(cell_nums)){
    # Find the pairwise comparisons of cells i and j
    cell_1 <- paste("X",cell_nums[i],sep="")
    cell_2 <- paste("X",cell_nums[j],sep="")
    data_12 <- bird_fb_turn_12[cell_1,cell_2]
    data_13 <- bird_fb_turn_13[cell_1,cell_2]
    data_14 <- bird_fb_turn_14[cell_1,cell_2]
    data_23 <- bird_fb_turn_23[cell_1,cell_2]
    data_24 <- bird_fb_turn_24[cell_1,cell_2]
    data_34 <- bird_fb_turn_34[cell_1,cell_2]
    # Take average value
    # 3/4 of pairwise comparisons will have a single value
    # 1/4 of pairwise comparisons will have three values, should be the same
    bird_fb_all[i,j] <- mean(na.omit(c(data_12,
                                      data_13,
                                      data_14,
                                      data_23,
                                      data_24,
                                      data_34)))
  }
  pb$tick()
}
colnames(bird_fb_all) <- paste("X", cell_nums, sep="")
row.names(bird_fb_all) <- paste("X", cell_nums, sep="")
Sys.time() - start # [1.41 hr]

# Save as a new distance matrix
write.table(bird_fb_all, "Data/bird_fb_turn.txt")
bird_fb_dist <- as.dist(bird_fb_all)
saveRDS(bird_fb_dist, "Data/bird_fb_turn_dist.RDS")

# ----- MEAN FUNCTIONAL BETA DIVERSITY TURNOVER ----------

library(usedist)

# Make list of cells with more species than functional axes
cell_list <- read.csv("Data/bird_SR_cells.csv") %>%
  filter(SR > 4) %>%
  mutate(cell_id = paste("X", cell_id, sep="")) 
# Load functional beta diversity turnover distance matrix
# Select only cells with more species than functional axes
fb_dist <- dist_subset(
  as.dist(read.table("Data/bird_fb_turn.txt")),
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
write.csv(fb_df,"Data/bird_mean_fb.csv")
