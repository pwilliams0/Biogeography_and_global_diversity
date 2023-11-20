#=========================================================#
# Create grid cells, create global environment variables, #
#      and create global bird/mammal/bat assemblages      #
#=========================================================#

library(tidyverse)
library(sf)
library(raster)
library(exactextractr)

# ----- CREATE 2x2 DEGREE GRID CELLS, SAVE REALM ----------
library(IceSat2R)
sf_use_s2(FALSE) # Avoid some issues of invalid geometries

# Load updated zoogeographic realms
# Downloaded from https://macroecology.ku.dk/resources/wallace
realms <- st_read("Data/Raw/newRealms/newRealms.shp") %>%
  st_make_valid()

# Create 2x2 degree global grid
grid <- IceSat2R::degrees_to_global_grid(minx = -180, maxx = 180,
                                         miny = -90, maxy = 90,
                                         degrees = 2) %>%
  dplyr::mutate(cell_id = row_number(),
                area = st_area(.)) %>%
  st_make_valid()

# Intersect with realms, keep realm name
# Keep cells >= 50% land
grid_realm_int <- grid %>%
  st_intersection(realms) %>%
  mutate(intersect_area = st_area(.),
         perc_land = as.numeric(intersect_area / area)) %>%
  dplyr::filter(perc_land >= 0.5) %>%
  dplyr::select(cell_id, Realm, perc_land) %>%
  st_drop_geometry()

# Create cells to be used in all analyses
cells <- grid %>% right_join(grid_realm_int, by = "cell_id") %>%
  dplyr::arrange(cell_id) %>%
  dplyr::mutate(cell_id = row_number()) # Convert to km2

# Save cells
st_write(cells, "Data/global_cells.shp")

# ----- ENVIRONMENT, LANDMASS AREA ----------

# Manually assign each cell to a landmass
# Manually enter area of each landmass using data from Global Islands dataset
# Count Eurasia as a single landmass
# Count North and South America as separate landmasses, divided by Isthmus of Panama
# Count Africa and Eurasia as separate landmasses, divided by Isthmus of Suez
# Save as "Landmass_area_cells.csv"

# ----- ENVIRONMENT, ELEVATION (MEAN AND RANGE) ----------
library(raster)
library(exactextractr)
cells <- st_read("Data/global_cells.shp")

# GMTED2010 from USGS
elev_rast <- raster("Data/Raw/mn30_grd/w001000.adf")

# Check that raster looks correct
plot(elev_rast)

elev_df <- cells
# Calculate mean elevation for each grid cell
elev_df$elev_mean <- exactextractr::exact_extract(elev_rast, cells, 'mean')
# Calculate max and min elevation to then calculate elevation range
elev_df$elev_max <- exactextractr::exact_extract(elev_rast, cells, 'max')
elev_df$elev_min <- exactextractr::exact_extract(elev_rast, cells, 'min')

elev_df <- elev_df %>% 
  mutate(elev_range = elev_max - elev_min) %>%
  dplyr::select(cell_id, elev_mean, elev_range) %>%
  st_drop_geometry()
write.csv(elev_df, "Data/elev_cells.csv")

# ----- ENVIRONMENT, PRESENT CLIMATE ----------
library(raster)
library(exactextractr)
cells <- st_read("Data/global_cells.shp")

# Data from worldclim.org, WorldClim 2.1
# Climate data from 1970-2000

# Load rasters of climate data, bioclimatic indicators
r01 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_1.tif")
r02 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_2.tif")
r03 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_3.tif")
r04 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_4.tif")
r05 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_5.tif")
r06 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_6.tif")
r07 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_7.tif")
r08 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_8.tif")
r09 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_9.tif")
r10 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_10.tif")
r11 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_11.tif")
r12 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_12.tif")
r13 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_13.tif")
r14 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_14.tif")
r15 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_15.tif")
r16 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_16.tif")
r17 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_17.tif")
r18 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_18.tif")
r19 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_19.tif")
# Create stack of all 19 bioclimatic variables
clim_stack <- raster::stack(r01, r02, r03, r04, r05, r06, r07, r08, r09, r10,
                            r11, r12, r13, r14, r15, r16, r17, r18, r19)

# Calculate mean climate value per polygon, repeat for each raster layer
mean_r01_val <- exactextractr::exact_extract(clim_stack[[1]], cells, 'mean')
mean_r02_val <- exactextractr::exact_extract(clim_stack[[2]], cells, 'mean')
mean_r03_val <- exactextractr::exact_extract(clim_stack[[3]], cells, 'mean')
mean_r04_val <- exactextractr::exact_extract(clim_stack[[4]], cells, 'mean')
mean_r05_val <- exactextractr::exact_extract(clim_stack[[5]], cells, 'mean')
mean_r06_val <- exactextractr::exact_extract(clim_stack[[6]], cells, 'mean')
mean_r07_val <- exactextractr::exact_extract(clim_stack[[7]], cells, 'mean')
mean_r08_val <- exactextractr::exact_extract(clim_stack[[8]], cells, 'mean')
mean_r09_val <- exactextractr::exact_extract(clim_stack[[9]], cells, 'mean')
mean_r10_val <- exactextractr::exact_extract(clim_stack[[10]], cells, 'mean')
mean_r11_val <- exactextractr::exact_extract(clim_stack[[11]], cells, 'mean')
mean_r12_val <- exactextractr::exact_extract(clim_stack[[12]], cells, 'mean')
mean_r13_val <- exactextractr::exact_extract(clim_stack[[13]], cells, 'mean')
mean_r14_val <- exactextractr::exact_extract(clim_stack[[14]], cells, 'mean')
mean_r15_val <- exactextractr::exact_extract(clim_stack[[15]], cells, 'mean')
mean_r16_val <- exactextractr::exact_extract(clim_stack[[16]], cells, 'mean')
mean_r17_val <- exactextractr::exact_extract(clim_stack[[17]], cells, 'mean')
mean_r18_val <- exactextractr::exact_extract(clim_stack[[18]], cells, 'mean')
mean_r19_val <- exactextractr::exact_extract(clim_stack[[19]], cells, 'mean')

# Save as dataframe
#   sqrt transform precipitation variables
bioclim_Present_cells <- as.data.frame(cells) %>%
  dplyr::select(cell_id) %>%
  mutate(BioClim01 = mean_r01_val,
         BioClim02 = mean_r02_val,
         BioClim03 = mean_r03_val,
         BioClim04 = mean_r04_val,
         BioClim05 = mean_r05_val,
         BioClim06 = mean_r06_val,
         BioClim07 = mean_r07_val,
         BioClim08 = mean_r08_val,
         BioClim09 = mean_r09_val,
         BioClim10 = mean_r10_val,
         BioClim11 = mean_r11_val,
         BioClim12 = sqrt(mean_r12_val),
         BioClim13 = sqrt(mean_r13_val),
         BioClim14 = sqrt(mean_r14_val),
         BioClim15 = mean_r15_val,
         BioClim16 = sqrt(mean_r16_val),
         BioClim17 = sqrt(mean_r17_val),
         BioClim18 = sqrt(mean_r18_val),
         BioClim19 = sqrt(mean_r19_val))

# Run PCA on all cells
clim_pca <- prcomp(bioclim_Present_cells[,c(2:20)],
                   center = TRUE,scale. = TRUE)
summary(clim_pca)

# Keep first 4 axes (explains >90% variance)
clim_Present_cells <- bioclim_Present_cells %>%
  mutate(clim_pca_1 = clim_pca[["x"]][,1],
         clim_pca_2 = clim_pca[["x"]][,2],
         clim_pca_3 = clim_pca[["x"]][,3],
         clim_pca_4 = clim_pca[["x"]][,4]) %>%
  dplyr::select(cell_id, clim_pca_1, clim_pca_2,
                clim_pca_3, clim_pca_4)

write.csv(clim_Present_cells, "Data/clim_Present_cells.csv")

# ----- ENVIRONMENT, *PAST CLIMATE (*not used in final analyses) ----------
library(raster)
library(exactextractr)
cells <- st_read("Data/global_cells.shp")

# ---------- Last Glacial Maximum (LGM) climate ----------

#https://www.worldclim.org/data/v1.4/paleo1.4.html

# Load rasters of climate data, bioclimatic indicators
r01 <- raster("Data/Raw/cclgmbi_10m/cclgmbi1.tif")
r02 <- raster("Data/Raw/cclgmbi_10m/cclgmbi2.tif")
r03 <- raster("Data/Raw/cclgmbi_10m/cclgmbi3.tif")
r04 <- raster("Data/Raw/cclgmbi_10m/cclgmbi4.tif")
r05 <- raster("Data/Raw/cclgmbi_10m/cclgmbi5.tif")
r06 <- raster("Data/Raw/cclgmbi_10m/cclgmbi6.tif")
r07 <- raster("Data/Raw/cclgmbi_10m/cclgmbi7.tif")
r08 <- raster("Data/Raw/cclgmbi_10m/cclgmbi8.tif")
r09 <- raster("Data/Raw/cclgmbi_10m/cclgmbi9.tif")
r10 <- raster("Data/Raw/cclgmbi_10m/cclgmbi10.tif")
r11 <- raster("Data/Raw/cclgmbi_10m/cclgmbi11.tif")
r12 <- raster("Data/Raw/cclgmbi_10m/cclgmbi12.tif")
r13 <- raster("Data/Raw/cclgmbi_10m/cclgmbi13.tif")
r14 <- raster("Data/Raw/cclgmbi_10m/cclgmbi14.tif")
r15 <- raster("Data/Raw/cclgmbi_10m/cclgmbi15.tif")
r16 <- raster("Data/Raw/cclgmbi_10m/cclgmbi16.tif")
r17 <- raster("Data/Raw/cclgmbi_10m/cclgmbi17.tif")
r18 <- raster("Data/Raw/cclgmbi_10m/cclgmbi18.tif")
r19 <- raster("Data/Raw/cclgmbi_10m/cclgmbi19.tif")
# Create stack of all 19 bioclimatic variables
clim_stack <- raster::stack(r01, r02, r03, r04, r05, r06, r07, r08, r09, r10,
                            r11, r12, r13, r14, r15, r16, r17, r18, r19)

# Calculate mean climate value per polygon, repeat for each raster layer
mean_r01_val <- exactextractr::exact_extract(clim_stack[[1]], cells, 'mean')
mean_r02_val <- exactextractr::exact_extract(clim_stack[[2]], cells, 'mean')
mean_r03_val <- exactextractr::exact_extract(clim_stack[[3]], cells, 'mean')
mean_r04_val <- exactextractr::exact_extract(clim_stack[[4]], cells, 'mean')
mean_r05_val <- exactextractr::exact_extract(clim_stack[[5]], cells, 'mean')
mean_r06_val <- exactextractr::exact_extract(clim_stack[[6]], cells, 'mean')
mean_r07_val <- exactextractr::exact_extract(clim_stack[[7]], cells, 'mean')
mean_r08_val <- exactextractr::exact_extract(clim_stack[[8]], cells, 'mean')
mean_r09_val <- exactextractr::exact_extract(clim_stack[[9]], cells, 'mean')
mean_r10_val <- exactextractr::exact_extract(clim_stack[[10]], cells, 'mean')
mean_r11_val <- exactextractr::exact_extract(clim_stack[[11]], cells, 'mean')
mean_r12_val <- exactextractr::exact_extract(clim_stack[[12]], cells, 'mean')
mean_r13_val <- exactextractr::exact_extract(clim_stack[[13]], cells, 'mean')
mean_r14_val <- exactextractr::exact_extract(clim_stack[[14]], cells, 'mean')
mean_r15_val <- exactextractr::exact_extract(clim_stack[[15]], cells, 'mean')
mean_r16_val <- exactextractr::exact_extract(clim_stack[[16]], cells, 'mean')
mean_r17_val <- exactextractr::exact_extract(clim_stack[[17]], cells, 'mean')
mean_r18_val <- exactextractr::exact_extract(clim_stack[[18]], cells, 'mean')
mean_r19_val <- exactextractr::exact_extract(clim_stack[[19]], cells, 'mean')

# Save as dataframe
#   Divide all temperature variables by 10 to get degrees C
#   sqrt transform precipitation variables
bioclim_LGM_cells <- as.data.frame(cells) %>%
  dplyr::select(cell_id) %>%
  mutate(BioClim01 = mean_r01_val/10,
         BioClim02 = mean_r02_val/10,
         BioClim03 = mean_r03_val,
         BioClim04 = mean_r04_val/10,
         BioClim05 = mean_r05_val/10,
         BioClim06 = mean_r06_val/10,
         BioClim07 = mean_r07_val/10,
         BioClim08 = mean_r08_val/10,
         BioClim09 = mean_r09_val/10,
         BioClim10 = mean_r10_val/10,
         BioClim11 = mean_r11_val/10,
         BioClim12 = sqrt(mean_r12_val),
         BioClim13 = sqrt(mean_r13_val),
         BioClim14 = sqrt(mean_r14_val),
         BioClim15 = mean_r15_val,
         BioClim16 = sqrt(mean_r16_val),
         BioClim17 = sqrt(mean_r17_val),
         BioClim18 = sqrt(mean_r18_val),
         BioClim19 = sqrt(mean_r19_val))

# Run PCA on all cells
clim_pca <- prcomp(bioclim_LGM_cells[,c(2:20)],
                   center = TRUE,scale. = TRUE)
summary(clim_pca)

# Keep first 4 axes (explains >90% variance)
clim_LGM_cells <- bioclim_LGM_cells %>%
  mutate(clim_pca_1 = clim_pca[["x"]][,1],
         clim_pca_2 = clim_pca[["x"]][,2],
         clim_pca_3 = clim_pca[["x"]][,3],
         clim_pca_4 = clim_pca[["x"]][,4]) %>%
  dplyr::select(cell_id, clim_pca_1, clim_pca_2,
                clim_pca_3, clim_pca_4)

write.csv(clim_LGM_cells, "Data/clim_LGM_cells.csv")

# ---------- Mid-Holocene (MidHol) climate ----------

#https://www.worldclim.org/data/v1.4/paleo1.4.html

# Load rasters of climate data, bioclimatic indicators
r01 <- raster("Data/Raw/ccmidbi_10m/ccmidbi1.tif")
r02 <- raster("Data/Raw/ccmidbi_10m/ccmidbi2.tif")
r03 <- raster("Data/Raw/ccmidbi_10m/ccmidbi3.tif")
r04 <- raster("Data/Raw/ccmidbi_10m/ccmidbi4.tif")
r05 <- raster("Data/Raw/ccmidbi_10m/ccmidbi5.tif")
r06 <- raster("Data/Raw/ccmidbi_10m/ccmidbi6.tif")
r07 <- raster("Data/Raw/ccmidbi_10m/ccmidbi7.tif")
r08 <- raster("Data/Raw/ccmidbi_10m/ccmidbi8.tif")
r09 <- raster("Data/Raw/ccmidbi_10m/ccmidbi9.tif")
r10 <- raster("Data/Raw/ccmidbi_10m/ccmidbi10.tif")
r11 <- raster("Data/Raw/ccmidbi_10m/ccmidbi11.tif")
r12 <- raster("Data/Raw/ccmidbi_10m/ccmidbi12.tif")
r13 <- raster("Data/Raw/ccmidbi_10m/ccmidbi13.tif")
r14 <- raster("Data/Raw/ccmidbi_10m/ccmidbi14.tif")
r15 <- raster("Data/Raw/ccmidbi_10m/ccmidbi15.tif")
r16 <- raster("Data/Raw/ccmidbi_10m/ccmidbi16.tif")
r17 <- raster("Data/Raw/ccmidbi_10m/ccmidbi17.tif")
r18 <- raster("Data/Raw/ccmidbi_10m/ccmidbi18.tif")
r19 <- raster("Data/Raw/ccmidbi_10m/ccmidbi19.tif")
# Create stack of all 19 bioclimatic variables
clim_stack <- raster::stack(r01, r02, r03, r04, r05, r06, r07, r08, r09, r10,
                            r11, r12, r13, r14, r15, r16, r17, r18, r19)

# Calculate mean climate value per polygon, repeat for each raster layer
mean_r01_val <- exactextractr::exact_extract(clim_stack[[1]], cells, 'mean')
mean_r02_val <- exactextractr::exact_extract(clim_stack[[2]], cells, 'mean')
mean_r03_val <- exactextractr::exact_extract(clim_stack[[3]], cells, 'mean')
mean_r04_val <- exactextractr::exact_extract(clim_stack[[4]], cells, 'mean')
mean_r05_val <- exactextractr::exact_extract(clim_stack[[5]], cells, 'mean')
mean_r06_val <- exactextractr::exact_extract(clim_stack[[6]], cells, 'mean')
mean_r07_val <- exactextractr::exact_extract(clim_stack[[7]], cells, 'mean')
mean_r08_val <- exactextractr::exact_extract(clim_stack[[8]], cells, 'mean')
mean_r09_val <- exactextractr::exact_extract(clim_stack[[9]], cells, 'mean')
mean_r10_val <- exactextractr::exact_extract(clim_stack[[10]], cells, 'mean')
mean_r11_val <- exactextractr::exact_extract(clim_stack[[11]], cells, 'mean')
mean_r12_val <- exactextractr::exact_extract(clim_stack[[12]], cells, 'mean')
mean_r13_val <- exactextractr::exact_extract(clim_stack[[13]], cells, 'mean')
mean_r14_val <- exactextractr::exact_extract(clim_stack[[14]], cells, 'mean')
mean_r15_val <- exactextractr::exact_extract(clim_stack[[15]], cells, 'mean')
mean_r16_val <- exactextractr::exact_extract(clim_stack[[16]], cells, 'mean')
mean_r17_val <- exactextractr::exact_extract(clim_stack[[17]], cells, 'mean')
mean_r18_val <- exactextractr::exact_extract(clim_stack[[18]], cells, 'mean')
mean_r19_val <- exactextractr::exact_extract(clim_stack[[19]], cells, 'mean')

# Save as dataframe
#   Divide all temperature variables by 10 to get degrees C
#   sqrt transform precipitation variables
bioclim_MidHol_cells <- as.data.frame(cells) %>%
  dplyr::select(cell_id) %>%
  mutate(BioClim01 = mean_r01_val/10,
         BioClim02 = mean_r02_val/10,
         BioClim03 = mean_r03_val,
         BioClim04 = mean_r04_val/10,
         BioClim05 = mean_r05_val/10,
         BioClim06 = mean_r06_val/10,
         BioClim07 = mean_r07_val/10,
         BioClim08 = mean_r08_val/10,
         BioClim09 = mean_r09_val/10,
         BioClim10 = mean_r10_val/10,
         BioClim11 = mean_r11_val/10,
         BioClim12 = sqrt(mean_r12_val),
         BioClim13 = sqrt(mean_r13_val),
         BioClim14 = sqrt(mean_r14_val),
         BioClim15 = mean_r15_val,
         BioClim16 = sqrt(mean_r16_val),
         BioClim17 = sqrt(mean_r17_val),
         BioClim18 = sqrt(mean_r18_val),
         BioClim19 = sqrt(mean_r19_val))

# Run PCA on all cells
clim_pca <- prcomp(bioclim_MidHol_cells[,c(2:20)],
                   center = TRUE,scale. = TRUE)
summary(clim_pca)

# Keep first 4 axes (explains >90% variance)
clim_MidHol_cells <- bioclim_MidHol_cells %>%
  mutate(clim_pca_1 = clim_pca[["x"]][,1],
         clim_pca_2 = clim_pca[["x"]][,2],
         clim_pca_3 = clim_pca[["x"]][,3],
         clim_pca_4 = clim_pca[["x"]][,4]) %>%
  dplyr::select(cell_id, clim_pca_1, clim_pca_2,
                clim_pca_3, clim_pca_4)

write.csv(clim_MidHol_cells, "Data/clim_MidHol_cells.csv")

# ---------- Climate distance from present ----------

# Combine climate data from all eras
bioclim_Present_cells <- bioclim_Present_cells %>%
  mutate(Era = "Present")
bioclim_LGM_cells <- bioclim_LGM_cells %>%
  mutate(Era = "LGM")
bioclim_MidHol_cells <- bioclim_MidHol_cells %>%
  mutate(Era = "MidHol")
bioclim_all <- bioclim_Present_cells %>%
  bind_rows(bioclim_LGM_cells) %>%
  bind_rows(bioclim_MidHol_cells)

# Run PCA on all cells/eras
clim_pca <- prcomp(bioclim_all[,c(2:20)],
                   center = TRUE,scale. = TRUE)
summary(clim_pca)

# Keep first 4 axes (explains >90% variance)
clim_All_cells <- bioclim_all %>%
  mutate(clim_pca_1 = clim_pca[["x"]][,1],
         clim_pca_2 = clim_pca[["x"]][,2],
         clim_pca_3 = clim_pca[["x"]][,3],
         clim_pca_4 = clim_pca[["x"]][,4])
# Create separate files for Present, LGM, MidHol
clim_Present_cells <- clim_All_cells %>%
  dplyr::filter(Era == "Present") %>%
  dplyr::select(cell_id, clim_pca_1, clim_pca_2,
                clim_pca_3, clim_pca_4)
clim_LGM_cells <- clim_All_cells %>%
  dplyr::filter(Era == "LGM") %>%
  dplyr::select(cell_id, clim_pca_1, clim_pca_2,
                clim_pca_3, clim_pca_4)
clim_MidHol_cells <- clim_All_cells %>%
  dplyr::filter(Era == "MidHol") %>%
  dplyr::select(cell_id, clim_pca_1, clim_pca_2,
                clim_pca_3, clim_pca_4)

# For each cell, calculate Euclidean distance in climate PCA space

# Present vs. LGM
Present_LGM_dist <- c()
for(i in 1:nrow(clim_Present_cells)){
  temp <- dplyr::bind_rows(clim_Present_cells[i,],
                           clim_LGM_cells[i,])
  Present_LGM_dist[i] <- as.numeric(dist(temp))
}

# Present vs. MidHol
Present_MidHol_dist <- c()
for(i in 1:nrow(clim_Present_cells)){
  temp <- dplyr::bind_rows(clim_Present_cells[i,],
                           clim_MidHol_cells[i,])
  Present_MidHol_dist[i] <- as.numeric(dist(temp))
}

# Save files
clim_diff_cells <- cells %>%
  mutate(clim_diff_LGM = Present_LGM_dist,
         clim_diff_MidHol = Present_MidHol_dist) %>%
  dplyr::select(cell_id, clim_diff_LGM,
                clim_diff_MidHol) %>%
  sf::st_drop_geometry()
write.csv(clim_diff_cells, "Data/clim_diff_cells.csv")

# ----- ENVIRONMENT, *PAST ICE COVER (*not used in final analyses) ----------
cells <- st_read("Data/global_cells.shp")

# pastclim package
# https://evolecolgroup.github.io/pastclim/index.html

# Install packages
install.packages("terra")
install.packages("curl")
devtools::install_github("EvolEcolGroup/pastclim")

library(pastclim)

# Download the data
pastclim::set_data_path()
pastclim::download_dataset(dataset="Beyer2020",
                           bio_variables = c("biome"))

# Create raster of ice cover for each time step
# 30,000 years ago to present

# Get time steps
time_steps <- pastclim::get_time_steps("Beyer2020")
time_steps <- time_steps[time_steps >= -30000]
# For each time step, create raster of ice cover
#   Then calculate % ice cover
ice_mat <- matrix(nrow=nrow(cells),
                  ncol=length(time_steps))
for(i in 1:length(time_steps)){
  ice_raster <- pastclim::region_slice(
    time_bp = time_steps[i],
    bio_variables = "biome",
    dataset = "Beyer2020")
  ice_raster$ice <- ice_raster$biome
  ice_raster$ice[ice_raster$ice != 28] <- FALSE
  ice_raster$ice[ice_raster$ice == 28] <- TRUE
  
  ice_mat[,i] <- exact_extract(ice_raster$ice, cells, 'mean')
}
# Calculate if cell ever >50% ice cover for any time step
ice_mat <- ice_mat > .5

ice_cells <- as.data.frame(cells) %>%
  dplyr::select(cell_id) %>%
  mutate(ice_cover = NA)
for(i in 1:nrow(ice_cells)){
  ice_cells$ice_cover[i] <- ifelse(sum(ice_mat[i,]) >= 1,
                                   1, 0)
}

# Save file
write.csv(ice_cells, "Data/ice_cells.csv")

# ----- ENVIRONMENT, *HISTORY OF LAND CONVERSION (*not used in final analyses) ----------
library(raster)
library(exactextractr)
library(ncdf4)
cells <- st_read("Data/global_cells.shp")

# KK10 dataset
# Kaplan, Jed O; Krumhardt, Kristen M (2011): The KK10 Anthropogenic Land Cover
#   Change scenario for the preindustrial Holocene, link to data in NetCDF
#   format. PANGAEA, https://doi.org/10.1594/PANGAEA.871369,
# Supplement to: Kaplan, Jed O; Krumhardt, Kristen M; Ellis, Erle C; Ruddiman,
#   William F; Lemmen, Carsten; Klein Goldewijk, Kees (2011): Holocene carbon
#   emissions as a result of anthropogenic land cover change. The Holocene,
#   21(5), 775-791, https://doi.org/10.1177/0959683610386983
# First significant land conversion (years before present)
#   Mean land conversion of cell > 20% (Ellis et al., 2013)
#   This threshold used in several studies

start <- Sys.time()
KK10 <- brick("Data/Raw/KK10.nc")
# Calculate mean land conversion per cell per year
# Split into groups of cells so computer can handle it
# Breaks based on where exact_extract has crashed
cells_1 <- cells[1:3241,]
cells_2 <- cells[3242:3375,]
cells_3 <- cells[3376:3718,]
KK10_df_1 <- exact_extract(KK10, cells_1, 'mean')
KK10_df_2 <- exact_extract(KK10, cells_2, 'mean')
KK10_df_3 <- exact_extract(KK10, cells_3, 'mean')
# Cells/years with > 20% land conversion
KK10_df <- dplyr::bind_rows(KK10_df_1, KK10_df_2) %>%
  dplyr::bind_rows(KK10_df_3)

end <- Sys.time()
end - start # 1.9 hours

KK10_df <- KK10_df > 0.2

# Find year where land conversion first exceeded 20%
# If never > 20%, set to zero
# Column name is years before 1950
KK10_list <- c()
for(i in 1:nrow(KK10_df)){
  index <- min(which(KK10_df[i,] == TRUE))
  KK10_list[i] <- ifelse(is.integer(index),
                         as.integer(
                           substring(
                             colnames(KK10_df)[index], 7)),
                    0)
}

KK10_cells <- cells %>%
  mutate(LandConvYears = KK10_list) %>%
  dplyr::select(cell_id, LandConvYears) %>%
  st_drop_geometry()
write.csv(KK10_cells, "Data/KK10_cells.csv")

# ----- ENVIRONMENT, *HUMAN IMPACT INDEX (*not used in final analyses) ----------
library(raster)
library(exactextractr)
cells <- st_read("Data/global_cells.shp")

# Human Impact Index, 2019 data
HII_2019 <- raster("Data/Raw/hii_2019-01-01.tif")

HII_df <- cells
HII_df$HII <- exactextractr::exact_extract(HII_2019, cells, 'mean')

HII_df <- HII_df %>% dplyr::select(cell_id, HII) %>%
  st_drop_geometry()
write.csv(HII_df, "Data/HII_cells.csv")

# ----- ASSEMBLAGES, BIRD ----------
library(progress)
sf_use_s2(FALSE) # Avoid some issues of invalid geometries
cells <- st_read("Data/global_cells.shp")

# ---------- List of birds in each cell ----------

# Original data for bird ranges from Birds Of The World dataset
#   bird <- st_read("Data/Raw/BOTW/BOTW.gdb", layer="All_species")
# Huge dataset, so loaded into QGIS, split into two files, and saved as shapefile
#   "BIRDS_resident" consists of resident ranges (Seasonal == 1)
#   "BIRDS_nonresident" consists of other ranges ((Seasonal != 1))

# Load both files for bird ranges
#   Exclude introduced or vagrant ranges (Origin is 3 or 4)
#   Exclude "passage", brief migratory stops (Seasonal is 4)
#   Exclude data errors (Origin is 0, one species)
bird_Nres <- st_read("Data/Raw/BOTW/BIRDS_nonresident.shp") %>%
  dplyr::filter(!ORIGIN %in% c(0,3,4),
                SEASONAL != 4) %>%
  dplyr::select(SCINAME) %>%
  tidyr::separate(SCINAME, c("genus", "spp"), " ") %>%
  dplyr::mutate(names_IUCN = paste(genus, spp, sep = "_")) %>%
  dplyr::select(names_IUCN)
bird_res <- st_read("Data/Raw/BOTW/BIRDS_resident.shp") %>%
  dplyr::filter(!ORIGIN %in% c(0,3,4),
                SEASONAL != 4) %>%
  dplyr::select(SCINAME) %>%
  tidyr::separate(SCINAME, c("genus", "spp"), " ") %>%
  dplyr::mutate(names_IUCN = paste(genus, spp, sep = "_")) %>%
  dplyr::select(names_IUCN)

# Save list of birds in dataset
write.csv(unique(c(bird_res$names_IUCN, bird_Nres$names_IUCN)),
          "Data/bird_names_IUCN.csv")

# Resident bird file too large, split in half
bird_res1 <- bird_res[1:round(nrow(bird_res)/2),]
bird_res2 <- bird_res[(1+round(nrow(bird_res)/2)):nrow(bird_res),]
st_write(bird_res1, "Data/Raw/BOTW/bird_res1.shp")
st_write(bird_res2, "Data/Raw/BOTW/bird_res2.shp")

# Create list of birds in each cell

# Non-residents [4 days 3 hours]
pb <- progress_bar$new(format=" running [:bar] :percent in :elapsed eta: :eta ",
                       total = nrow(cells), clear = FALSE, width = 100)
bird_list_cell_Nres <- c()
for(i in 1:nrow(cells)){
  intersect  <- suppressMessages(
    st_intersects(bird_Nres, cells[i,], sparse = FALSE))
  intersect_vec <- which(intersect == TRUE)
  spp_list <- bird_Nres[intersect_vec, 1]
  st_geometry(spp_list) <- NULL
  bird_list_cell_Nres[i] <- spp_list
  pb$tick()
}
# Save list of birds in each cell, non-residents
saveRDS(bird_list_cell_Nres, "bird_list_cell_Nres.RDS")

# Residents 1 [3 days 13 hours]
pb <- progress_bar$new(format=" running [:bar] :percent in :elapsed eta: :eta ",
                       total = nrow(cells), clear = FALSE, width = 100)
bird_list_cell_res1 <- c()
for(i in 1:nrow(cells)){
  intersect  <- suppressMessages(
    st_intersects(bird_res1, cells[i,], sparse = FALSE))
  intersect_vec <- which(intersect == TRUE)
  spp_list <- bird_res1[intersect_vec, 1]
  st_geometry(spp_list) <- NULL
  bird_list_cell_res1[i] <- spp_list
  pb$tick()
}
# Save list of birds in each cell, residents 1
saveRDS(bird_list_cell_res1, "bird_list_cell_res1.RDS")

# Residents 2 [3 days 8 hours]
pb <- progress_bar$new(format=" running [:bar] :percent in :elapsed eta: :eta ",
                       total = nrow(cells), clear = FALSE, width = 100)
bird_list_cell_res2 <- c()
for(i in 1:nrow(cells)){
  intersect  <- suppressMessages(
    st_intersects(bird_res2, cells[i,], sparse = FALSE))
  intersect_vec <- which(intersect == TRUE)
  spp_list <- bird_res2[intersect_vec, 1]
  st_geometry(spp_list) <- NULL
  bird_list_cell_res2[i] <- spp_list
  pb$tick()
}
# Save list of birds in each cell, residents 2
saveRDS(bird_list_cell_res, "bird_list_cell_res2.RDS")

# ---------- Birds in each cell, with trait data ----------

# Compile full list of birds in each cell
bird_list_cell_Nres <- readRDS("Data/bird_list_cell_Nres.RDS")
bird_list_cell_res1 <- readRDS("Data/bird_list_cell_res1.RDS")
bird_list_cell_res2 <- readRDS("Data/bird_list_cell_res2.RDS")
bird_list_cell <- list()
for(i in 1:length(bird_list_cell_Nres)){
  bird_list_cell[[i]] <- c(bird_list_cell_Nres[[i]],
                           bird_list_cell_res1[[i]],
                           bird_list_cell_res2[[i]])
}

# Import Elton traits data for birds
# Errors in reading original text file, read CSV instead
elton_bird <- read.csv(file="Data/Raw/BirdFuncDat.csv", header=TRUE,
                       stringsAsFactors = FALSE) %>%
  separate(Scientific, c("Genus_only", "Species_only"), sep = " ") %>%
  mutate(names_elton = paste(Genus_only, Species_only, sep = "_")) %>%
  dplyr::group_by(names_elton) %>%
  mutate(Diet.Vert = sum(Diet.Vend,
                         Diet.Vect,
                         Diet.Vunk),
         Diet.Fish = Diet.Vfish,
         Diet.Herb = Diet.PlantO,
         ForStrat.aquatic = sum(ForStrat.watbelowsurf,
                                ForStrat.wataroundsurf),
         ForStrat.canopy = sum(ForStrat.canopy,
                               ForStrat.midhigh),
         LogBodyMass = log10(BodyMass.Value)) %>%
  ungroup() %>%
  dplyr::select(names_elton, BLFamilyLatin, IOCOrder,
                Diet.Inv, Diet.Vert, Diet.Fish, Diet.Scav,
                Diet.Fruit, Diet.Nect, Diet.Seed, Diet.Herb,
                ForStrat.aquatic,
                ForStrat.ground,
                ForStrat.understory,
                ForStrat.canopy,
                ForStrat.aerial,
                PelagicSpecialist,
                LogBodyMass)

# Match species names
# bird_names_cells.csv created by manual matching species names
bird_traits <- read.csv(file="Data/Raw/bird_names_cells.csv", header=TRUE,
                        stringsAsFactors=FALSE) %>%
  left_join(elton_bird, by="names_elton") %>%
  dplyr::select(-diff_IUCN_elton, -diff_IUCN_phylo)

# Get traits of each species in each grid cell
bird_cell_df <- enframe(bird_list_cell,
                        name="cell_id",
                        value="names_IUCN") %>%
  unnest(cols=c(names_IUCN)) %>%
  left_join(bird_traits,by="names_IUCN") %>%
  dplyr::filter(PelagicSpecialist == 0) %>% # Exclude pelagic specialist
  dplyr::group_by(cell_id, names_IUCN) %>% slice(1)

# Get dataframe of species in each grid cell
write.csv(bird_cell_df, "Data/bird_cell_df.csv")

# ----- ASSEMBLAGES, MAMMAL & BAT ----------
library(progress)
sf_use_s2(FALSE) # Avoid some issues of invalid geometries
cells <- st_read("Data/global_cells.shp")

# ---------- List of mammals in each cell ----------

# Load terrestrial mammal ranges
#   Exclude non-terrestrial ranges
#   Exclude introduced ranges
mamm <- st_read("Data/Raw/MAMMALS/MAMMALS.shp") %>%
  dplyr::filter(terrestial == "t") %>%
  dplyr::filter(!legend %in% c("Extant & Introduced (resident)",
                               "Extinct & Introduced",
                               "Introduced",
                               "Probably Extant & Introduced (resident)")) %>%
  dplyr::select(binomial) %>%
  tidyr::separate(binomial, c("genus", "spp"), " ") %>%
  dplyr::mutate(species = paste(genus, spp, sep = "_")) %>%
  dplyr::select(species)

# Save list of mammals in dataset
write.csv(unique(mamm$species), "Data/mamm_names_IUCN.csv")

# Create list of mammals in each cell
pb <- progress_bar$new(format=" running [:bar] :percent in :elapsed eta: :eta ",
                       total = nrow(cells), clear = FALSE, width = 100)
mamm_list_cell <- c()
for(i in 1:nrow(cells)){
  intersect  <- suppressMessages(
    st_intersects(mamm, cells[i,], sparse = FALSE))
  intersect_vec <- which(intersect == TRUE)
  spp_list <- mamm[intersect_vec, 1]
  st_geometry(spp_list) <- NULL
  mamm_list_cell[i] <- spp_list
  pb$tick()
}

# Save list of mammals in each cell
saveRDS(mamm_list_cell, "Data/mamm_list_cell.RDS")

# ---------- Mammals in each cell, with trait data ----------

# Import Elton traits data for mammals
elton_mamm <- read.table(file="Data/Raw/MamFuncDat.txt", sep="\t", header=TRUE,
                         stringsAsFactors = FALSE) %>%
  dplyr::filter(MSW3_ID > 0) %>% # Remove empty rows
  separate(Scientific, c("Genus_only", "Species_only"), sep = " ") %>%
  mutate(names_elton = paste(Genus_only, Species_only, sep = "_")) %>%
  dplyr::group_by(names_elton) %>%
  mutate(Diet.Vert = sum(Diet.Vend,
                         Diet.Vect,
                         Diet.Vunk),
         Diet.Fish = Diet.Vfish,
         Diet.Herb = Diet.PlantO,
         ForStrat = ForStrat.Value,
         LogBodyMass = log10(BodyMass.Value)) %>%
  dplyr::select(names_elton, MSWFamilyLatin,
                Diet.Inv, Diet.Vert, Diet.Fish, Diet.Scav,
                Diet.Fruit, Diet.Nect, Diet.Seed, Diet.Herb,
                ForStrat, LogBodyMass)

# Match species names
# mamm_names_cells.csv creating by manual matching species names
mamm_traits <- read.csv(file="Data/Raw/mamm_names_cells.csv", header=TRUE,
                        stringsAsFactors=FALSE) %>%
  left_join(elton_mamm, by="names_elton") %>%
  dplyr::select(-diff_IUCN_elton, -diff_IUCN_phylo)

# Get traits of each species in each grid cell
mamm_cell_df <- enframe(readRDS("Data/mamm_list_cell.RDS"),
                        name="cell_id",
                        value="names_IUCN") %>%
  unnest(cols=c(names_IUCN)) %>%
  left_join(mamm_traits,by="names_IUCN") %>%
  dplyr::filter(ForStrat != "M") %>% # Exclude marine mammals
  dplyr::group_by(cell_id, names_IUCN) %>% slice(1)

# Get dataframe of species in each grid cell
write.csv(mamm_cell_df, "Data/mamm_cell_df.csv")

# ---------- Bat assemblages ----------

# Create list of bat families
bat_fams <- c("Pteropodidae", "Rhinolophidae",
              "Hipposideridae", "Megadermatidae",
              "Rhinopomatidae", "Craseonycteridae",
              "Emballonuridae", "Nycteridae",
              "Myzopodidae", "Mystacinidae",
              "Phyllostomidae", "Mormoopidae",
              "Noctilionidae", "Furipteridae",
              "Thyropteridae", "Natalidae",
              "Molossidae", "Vespertilionidae")

# Save dataframe of bat species in each grid cell
bat_cell_df <- read.csv("Data/mamm_cell_df.csv") %>%
  filter(MSWFamilyLatin %in% bat_fams)
write.csv(bat_cell_df, "Data/bat_cell_df.csv")
