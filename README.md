# Deep biogeographic barriers explain divergent global vertebrate communities

### Peter Jeffrey Williams, Elise F. Zipkin, Jedediah F. Brodie

---------------------------------

## Abstract

Regional and historical processes can lead to differences in biodiversity across biogeographic regions. However, the influence of biogeography on global biodiversity patterns remains unclear. Although attempts to understand biogeographic effects on biodiversity generally treat regions as discrete units, differences in biogeographic distances among species assemblages are actually continuous. Here, we partition variance in the species richness, functional richness, phylogenetic diversity, and mean functional beta diversity turnover of mammal and bird assemblages across the world into biogeographic and environmental components using phylobetadiversity, a metric that allows us to calculate continuous biogeographic distances among assemblages. We found that, on average, environment (climate, elevation, topography, and landmass area) explained global diversity patterns much better than phylobetadiversity, especially for birds. However, mammals in deeply isolated regions such as Australia were strongly influenced by biogeographic history, with lower taxonomic and functional richness than predicted by environment alone. Mammal assemblages in biogeographically isolated regions contained unique combinations of functional traits as well. Differences in Neotropical vs. Paleotropical bat functional richness also reflect the influence of biogeographic isolation. These examples show how deep, long-lasting biogeographic barriers can lead to divergent diversity patterns, against the backdrop of environmental determinism in diversity across most of the world.

## Code

1. [0_DataPrep_AllTaxa.R](Code/0_DataPrep_AllTaxa.R): Create global 2x2 degree grid cells, calculate environment variables for each grid cell, and assemblages for each grid cell. Environment variables include mean elevation, elevation range, present-day climate, past climate, distance between present-day climate and past climate, past ice cover, years since significant land conversion, and Human Impact Index. Assemblages are lists of all bird, mammal, or bat species with range maps that overlap a given cell.
2. 1_SR_PD_FR_XXXX.R: Calculate species richness (SR), phylogenetic diversity (PD), and functional richness (FR) for each grid cell. Separate scripts for birds, mammals, and bats.
    - [1_SR_PD_FR_Bird.R](Code/1_SR_PD_FR_Bird.R)
    - [1_SR_PD_FR_Mamm.R](Code/1_SR_PD_FR_Mamm.R)
    - [1_SR_PD_FR_Bat.R](Code/1_SR_PD_FR_Bat.R)
3. 2_PhyloB_Bird_XX.R: Calculate phylobetadiversity for birds on subsets of grid cells. Calculating all pairwise comparisions is computationally intractable, so grid cells are split into 5 groups.
    - [2_PhyloB_Bird_12.R](Code/2_PhyloB_Bird_12.R): Bird phylobetadiversity for cells in groups 1 and 2.
    - [2_PhyloB_Bird_13.R](Code/2_PhyloB_Bird_13.R): Bird phylobetadiversity for cells in groups 1 and 3.
    - [2_PhyloB_Bird_14.R](Code/2_PhyloB_Bird_14.R): Bird phylobetadiversity for cells in groups 1 and 4.
    - [2_PhyloB_Bird_15.R](Code/2_PhyloB_Bird_15.R): Bird phylobetadiversity for cells in groups 1 and 5.
    - [2_PhyloB_Bird_23.R](Code/2_PhyloB_Bird_23.R): Bird phylobetadiversity for cells in groups 2 and 3.
    - [2_PhyloB_Bird_24.R](Code/2_PhyloB_Bird_24.R): Bird phylobetadiversity for cells in groups 2 and 4.
    - [2_PhyloB_Bird_25.R](Code/2_PhyloB_Bird_25.R): Bird phylobetadiversity for cells in groups 2 and 5.
    - [2_PhyloB_Bird_34.R](Code/2_PhyloB_Bird_34.R): Bird phylobetadiversity for cells in groups 3 and 4.
    - [2_PhyloB_Bird_35.R](Code/2_PhyloB_Bird_35.R): Bird phylobetadiversity for cells in groups 3 and 5.
    - [2_PhyloB_Bird_45.R](Code/2_PhyloB_Bird_45.R): Bird phylobetadiversity for cells in groups 4 and 5.
4. 2_PhyloB_Mamm_XX.R: Calculate phylobetadiversity for mammals on subsets of grid cells. Calculating all pairwise comparisions is computationally intractable, so grid cells are split into 4 groups.
    - [2_PhyloB_Mamm_12.R](Code/2_PhyloB_Mamm_12.R): Mammal phylobetadiversity for cells in groups 1 and 2.
    - [2_PhyloB_Mamm_13.R](Code/2_PhyloB_Mamm_13.R): Mammal phylobetadiversity for cells in groups 1 and 3.
    - [2_PhyloB_Mamm_14.R](Code/2_PhyloB_Mamm_14.R): Mammal phylobetadiversity for cells in groups 1 and 4.
    - [2_PhyloB_Mamm_23.R](Code/2_PhyloB_Mamm_23.R): Mammal phylobetadiversity for cells in groups 2 and 3.
    - [2_PhyloB_Mamm_24.R](Code/2_PhyloB_Mamm_24.R): Mammal phylobetadiversity for cells in groups 2 and 4.
    - [2_PhyloB_Mamm_34.R](Code/2_PhyloB_Mamm_34.R): Mammal phylobetadiversity for cells in groups 3 and 4.
5. 2_PhyloB_Bat_XX.R: Calculate phylobetadiversity for bats on subsets of grid cells. Calculating all pairwise comparisions is computationally intractable, so grid cells are split into 3 groups.
    - [2_PhyloB_Bat_12.R](Code/2_PhyloB_Bat_12.R): Bat phylobetadiversity for cells in groups 1 and 2.
    - [2_PhyloB_Bat_13.R](Code/2_PhyloB_Bat_13.R): Bat phylobetadiversity for cells in groups 1 and 3.
    - [2_PhyloB_Bat_23.R](Code/2_PhyloB_Bat_23.R): Bat phylobetadiversity for cells in groups 2 and 3.
6. 2_FuncB_Bird_XX.R: Calculate functional beta diversity for birds on subsets of grid cells. Calculating all pairwise comparisions is computationally intractable, so grid cells are split into 4 groups.
    - [2_FuncB_Bird_12.R](Code/2_FuncB_Bird_12.R): Bird functional beta diversity for cells in groups 1 and 2.
    - [2_FuncB_Bird_13.R](Code/2_FuncB_Bird_13.R): Bird functional beta diversity for cells in groups 1 and 3.
    - [2_FuncB_Bird_14.R](Code/2_FuncB_Bird_14.R): Bird functional beta diversity for cells in groups 1 and 4.
    - [2_FuncB_Bird_23.R](Code/2_FuncB_Bird_23.R): Bird functional beta diversity for cells in groups 2 and 3.
    - [2_FuncB_Bird_24.R](Code/2_FuncB_Bird_24.R): Bird functional beta diversity for cells in groups 2 and 4.
    - [2_FuncB_Bird_34.R](Code/2_FuncB_Bird_34.R): Bird functional beta diversity for cells in groups 3 and 4.
7. 3_PhyloB_XXXX_full.R: Combine files to create full phylobetadiversity turnover distance matrix, and save NMDS axes. Separate scripts for birds, mammals, and bats.
    - [3_PhyloB_Bird_full.R](Code/3_PhyloB_Bird_full.R)
    - [3_PhyloB_Mamm_full.R](Code/3_PhyloB_Mamm_full.R)
    - [3_PhyloB_Bat_full.R](Code/3_PhyloB_Bat_full.R)
8. [3_FuncB_Bird_full.R](Code/3_FuncB_Bird_full.R): Combine files to create functional beta diversity turnover distance matrix for birds, then calculate mean functional beta diversity turnover for birds.
9. 3_FuncB_XXXX.R: Calculate functional beta diversity, then calculate mean functional beta diversity turnover. Separate scripts for mammals and bats.
    - [3_FuncB_Mamm.R](Code/3_FuncB_Mamm.R)
    - [3_FuncB_Bat.R](Code/3_FuncB_Bat.R)
10. 4_Analysis_XXXX.R: Run models for species richness, phylogenetic diversity, functional richness, and mean functional beta diversity turnover. Save variance partitioning results. Calculate differences in residuals and save figures (Fig. 3 and Extended Data Fig. 1). Separate scripts for birds, mammals, and bats.
    - [4_Analysis_Bird.R](Code/4_Analysis_Bird.R)
    - [4_Analysis_Mamm.R](Code/4_Analysis_Mamm.R)
    - [4_Analysis_Bat.R](Code/4_Analysis_Bat.R)
11. [5_BarPlots.R](Code/5_BarPlots.R): Create barplots of variance partitioning results (Fig. 2, Fig. 6a, and Extended Data Fig. 7).
12. [5_DiversityMaps.R](Code/5_DiversityMaps.R): Create maps of global diversity patterns of species richness, phylogenetic diversity, functional richness, and mean functional beta diversity turnover for birds, mammals, and bats (Fig. 4, Fig. 6b, and Extended Data Fig. 2).
13. [5_PhyloB_Figures.R](Code/5_PhyloB_Figures.R): Create phylobetadiversity figures, including NMDS plots using 2 axes (Fig. 1), NMDS plots using 3 axes (Extended Data Fig. 3), and scree plots (Extended Data Fig. 5).
14. [5_FuncSpace_Figures.R](Code/5_FuncSpace_Figures.R): Create figures showing overlap in functional space (Fig. 5 and Extended Data Fig. 4).
15. [5_PhyloSignal.R](Code/5_PhyloSignal.R): Calculate phylogenetic signal of continuous traits (Extended Data Table 1).
16. [5_AdditionalVariables.R](Code/5_AdditionalVariables.R): Calculate additional variance explained by environment variables besides climate (Extended Data Table 2).
17. [5_PresentVsPastClimate.R](Code/5_PresentVsPastClimate.R): Calculate variance explained using climate data from present day, Last Glacial Maximum, and mid-Holocene (Extended Data Table 3).

## Data

### Raw
1. [newRealms](Data/Raw/newRealms): Updated zoogeographic realms. Available from the Center for Macroecology, Evolution and Climate at the University of Copenhagen.
2. BOTW: BirdLife range map data for birds (not publicly available). Due to the large size, we split the original GDB file into two shapefiles prior to data processing. "BIRDS_resident.shp" consists of resident ranges (Seasonal == 1), and "BIRDS_nonresident.shp" consists of other ranges (Seasonal != 1). Data downloaded in 2018. Available from BirdLife International.
3. MAMMALS: IUCN range map data for mammals (not publicly available). Data downloaded in 2018. Available from the IUCN Red List.
4. [BirdFuncDat.csv](Data/Raw/BirdFuncDat.csv): EltonTraits 1.0 data for birds. Available from Figshare. Converted to CSV due to issues with original TXT file. 
5. [MamFuncDat.txt](Data/Raw/MamFuncDat.txt): EltonTraits 1.0 data for mammals. Available from Figshare.
6. [bird_names_cells.csv](Data/Raw/bird_names_cells.csv): Bird taxonomies for range map data (names_ICUN), trait data (names_elton), and phylogeny (names_phylo), manually aligned.
7. [mamm_names_cells.csv](Data/Raw/mamm_names_cells.csv): Mammal taxonomies for range map data (names_ICUN), trait data (names_elton), and phylogeny (names_phylo), manually aligned.
8. [bird_traits_cat.csv](Data/Raw/bird_traits_cat.csv): Trait type for each bird trait.
9. [mamm_traits_cat.csv](Data/Raw/mamm_traits_cat.csv): Trait type for each mammal trait.
10. [Stage2_MayrParSho_Hackett_mean phylogeny_ultrametric.tre](Data/Raw/Stage2_MayrParSho_Hackett_mean%20phylogeny_ultrametric.tre): Bird phylogeny, ultrametric consensus tree based on 1,000 credible phylogenies available from VertLife.
11. [Upham_mean_phylogeny_ultrametric.tree](Data/Raw/Upham_mean_phylogeny_ultrametric.tree): Mammal phylogeny, ultrametric consensus tree based on 1,000 credible phylogenies available from VertLife.
12. [Landmass_area_cells.csv](Data/Raw/Landmass_area_cells.csv): Landmass area for each cell. Manually entered using the Global Islands dataset available from the United States Geological Survey.
13. mn30_grd: Global Multi-resolution Terrain Elevation Data (GMTED2010). Files too large to push to GitHub. Available from the United States Geological Survey.
14. [wc2.1_10m_bio](Data/Raw/wc2.1_10m_bio): Bioclimatic variables, present day. Available from WorldClim.
15. [cclgmbi_10m](Data/Raw/cclgmbi_10m): WorldClim 1.4 downscaled paleo climate, bioclimatic variables, Last Glacial Maximum. Available from WorldClim.
16. [ccmidbi_10m](Data/Raw/ccmidbi_10m): WorldClim 1.4 downscaled paleo climate, bioclimatic variables, Mid Holocene. Available from WorldClim.
17. KK10.nc: Anthropogenic land cover change over time. File too large to push to GitHub. Available from PANGAEA.
18. hii_2019-01-01.tif: Human Impact Index, 2019 data. File too large to push to GitHub. Available from the Wildlife Conservation Society.

### Intermediate data products (for running analyses and making figures)

1. [global_cells.shp](Data/global_cells.shp): 2x2 degree grid cells used in all analyses.
2. [clim_Present_cells.csv](Data/clim_Present_cells.csv): Climate PCA coordinates for 4 axes for each grid cell.
3. [elev_cells.csv](Data/elev_cells.csv): Mean elevation and elevation range of each grid cell.
4. XXXX_cell_df.csv: Lists of species in each grid cell. Separate files for birds, mammals, and bats.
    - bird_cell_df.csv: File too large to push to GitHub.
    - [mamm_cell_df.csv](Data/mamm_cell_df.csv)
    - [bat_cell_df.csv](Data/bat_cell_df.csv)
5. XXXX_SR_cells.csv: Species richness of each grid cell. Separate files for birds, mammals, and bats.
    - [bird_SR_cells.csv](Data/bird_SR_cells.csv)
    - [mamm_SR_cells.csv](Data/mamm_SR_cells.csv)
    - [bat_SR_cells.csv](Data/bat_SR_cells.csv)
6. XXXX_PD_cells.csv: Phylogenetic diversity of each grid cell. Separate files for birds, mammals, and bats.
    - [bird_PD_cells.csv](Data/bird_PD_cells.csv)
    - [mamm_PD_cells.csv](Data/mamm_PD_cells.csv)
    - [bat_PD_cells.csv](Data/bat_PD_cells.csv)
7. XXXX_FR_cells.csv: Functional richness of each grid cell. Separate files for birds, mammals, and bats.
    - [bird_FR_cells.csv](Data/bird_FR_cells.csv)
    - [mamm_FR_cells.csv](Data/mamm_FR_cells.csv)
    - [bat_FR_cells.csv](Data/bat_FR_cells.csv)
8. XXXX_mean_fb.csv: Mean functional beta diversity turnover of each grid cell. Separate files for birds, mammals, and bats.
    - [bird_mean_fb.csv](Data/bird_mean_fb.csv)
    - [mamm_mean_fb.csv](Data/mamm_mean_fb.csv)
    - [bat_mean_fb.csv](Data/bat_mean_fb.csv)
9. XXXX_pb_NMDS.csv: Phylobetadiversity turnover NMDS coordinates for 3 axes for each grid cell. Separate files for birds, mammals, and bats.
    - [bird_pb_NMDS.csv](Data/bird_pb_NMDS.csv)
    - [mamm_pb_NMDS.csv](Data/mamm_pb_NMDS.csv)
    - [bat_pb_NMDS.csv](Data/bat_pb_NMDS.csv)
10. XXXX_pb_dist.RDS: Distance matrix of phylobetadiversity turnover. Separate files for birds, mammals, and bats.
    - [bird_pb_dist.RDS](Data/bird_pb_dist.RDS)
    - [mamm_pb_dist.RDS](Data/mamm_pb_dist.RDS)
    - [bat_pb_dist.RDS](Data/bat_pb_dist.RDS)
11. XXXX_fb.RDS: Output from functional beta diversity calculation, including distance matrices for total beta diversity, turnover, and nestedness. Separate files for mammals and bats. Files too large to push to GitHub.
    - mamm_fb.RDS
    - bat_fb.RDS
