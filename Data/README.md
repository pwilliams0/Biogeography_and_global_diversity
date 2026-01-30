## [Data](Data)

### Raw data
1. **newRealms**: Zoogeographic realms. Available from the Center for Macroecology, Evolution and Climate at the University of Copenhagen ([https://macroecology.ku.dk/resources/wallace](https://macroecology.ku.dk/resources/wallace)).
2. **BOTW**: BirdLife range map data for birds. Due to the large size, we split the original GDB file into two shapefiles prior to data processing. 'BIRDS_resident.shp' consisted of resident ranges (Seasonal == 1), and 'BIRDS_nonresident.shp' consisted of other ranges (Seasonal != 1). Data downloaded in 2018. Available upon request from BirdLife International ([http://datazone.birdlife.org/species/requestdis](http://datazone.birdlife.org/species/requestdis)).
3. **MAMMALS**: IUCN range map data for mammals. Data downloaded in 2018. Available from the IUCN Red List ([https://www.iucnredlist.org/resources/spatial-data-download](https://www.iucnredlist.org/resources/spatial-data-download)).
4. **BirdFuncDat.csv**: EltonTraits 1.0 data for birds. Converted to CSV due to issues with original TXT file. TXT file available from Figshare ([https://doi.org/10.6084/m9.figshare.c.3306933.v1](https://doi.org/10.6084/m9.figshare.c.3306933.v1))
5. **MamFuncDat.txt**: EltonTraits 1.0 data for mammals. Available from Figshare ([https://doi.org/10.6084/m9.figshare.c.3306933.v1](https://doi.org/10.6084/m9.figshare.c.3306933.v1)).
6. **[bird_names_cells.csv](Data/Raw/bird_names_cells.csv)**: Bird taxonomies for range map data (names_ICUN), trait data (names_elton), and phylogeny (names_phylo), manually aligned.
7. **[mamm_names_cells.csv](Data/Raw/mamm_names_cells.csv)**: Mammal taxonomies for range map data (names_ICUN), trait data (names_elton), and phylogeny (names_phylo), manually aligned.
8. **[bird_traits_cat.csv](Data/Raw/bird_traits_cat.csv)**: Trait type for each bird trait.
9. **[mamm_traits_cat.csv](Data/Raw/mamm_traits_cat.csv)**: Trait type for each mammal trait.
10. **Stage2_MayrParSho_Hackett_mean phylogeny_ultrametric.tre**: Bird phylogeny, ultrametric consensus tree created based on 1,000 credible phylogenies. Credible phylogenies available from BirdTree ([https://birdtree.org/subsets/](https://birdtree.org/subsets/)).
11. **Upham_mean_phylogeny_ultrametric.tree**: Mammal phylogeny, ultrametric consensus tree based on 1,000 credible phylogenies. Credible phylogenies available from VertLife ([http://vertlife.org/phylosubsets/](http://vertlife.org/phylosubsets/)).
12. **[Landmass_area_cells.csv](Data/Raw/Landmass_area_cells.csv)**: Landmass area for each cell. Manually entered using the Global Islands dataset available from the United States Geological Survey ([https://rmgsc.cr.usgs.gov/gie/gie.shtml](https://rmgsc.cr.usgs.gov/gie/gie.shtml)).
13. **mn30_grd**: Global Multi-resolution Terrain Elevation Data (GMTED2010). Available from the United States Geological Survey ([https://topotools.cr.usgs.gov/gmted_viewer/gmted2010_global_grids.php](https://topotools.cr.usgs.gov/gmted_viewer/gmted2010_global_grids.php)).
14. **wc2.1_10m_bio**: Bioclimatic variables, present day. Available from WorldClim ([https://www.worldclim.org/data/worldclim21.html](https://www.worldclim.org/data/worldclim21.html)).
15. **cclgmbi_10m**: WorldClim 1.4 downscaled paleo climate, bioclimatic variables, Last Glacial Maximum. Available from WorldClim ([https://www.worldclim.org/data/v1.4/paleo1.4.html](https://www.worldclim.org/data/v1.4/paleo1.4.html)).
16. **ccmidbi_10m**: WorldClim 1.4 downscaled paleo climate, bioclimatic variables, Mid Holocene. Available from WorldClim ([https://www.worldclim.org/data/v1.4/paleo1.4.html](https://www.worldclim.org/data/v1.4/paleo1.4.html)).
17. **KK10.nc**: Anthropogenic land cover change over time. Available from PANGAEA ([https://doi.pangaea.de/10.1594/PANGAEA.871369](https://doi.pangaea.de/10.1594/PANGAEA.871369)).
18. **hii_2019-01-01.tif**: Human Impact Index, 2019 data. Available from the Wildlife Conservation Society ([https://wcshumanfootprint.org/data-access](https://wcshumanfootprint.org/data-access)).

### Intermediate data products

1. **[global_cells.shp](Data/global_cells.shp)**: 2x2 degree latitude-longitude grid cells used in all analyses.
2. **[clim_Present_cells.csv](Data/clim_Present_cells.csv)**: Climate PCA coordinates for 4 axes for each grid cell.
3. **[elev_cells.csv](Data/elev_cells.csv)**: Mean elevation and elevation range of each grid cell.
4. **XXXX_cell_df.csv**: Lists of species in each grid cell. Separate files for birds, mammals, and bats.
    - bird_cell_df.csv: File too large to push to GitHub.
    - [mamm_cell_df.csv](Data/mamm_cell_df.csv)
    - [bat_cell_df.csv](Data/bat_cell_df.csv)
5. **XXXX_SR_cells.csv**: Species richness of each grid cell. Separate files for birds, mammals, and bats.
    - [bird_SR_cells.csv](Data/bird_SR_cells.csv)
    - [mamm_SR_cells.csv](Data/mamm_SR_cells.csv)
    - [bat_SR_cells.csv](Data/bat_SR_cells.csv)
6. **XXXX_PD_cells.csv**: Phylogenetic diversity of each grid cell. Separate files for birds, mammals, and bats.
    - [bird_PD_cells.csv](Data/bird_PD_cells.csv)
    - [mamm_PD_cells.csv](Data/mamm_PD_cells.csv)
    - [bat_PD_cells.csv](Data/bat_PD_cells.csv)
7. **XXXX_FR_cells.csv**: Functional richness of each grid cell. Separate files for birds, mammals, and bats.
    - [bird_FR_cells.csv](Data/bird_FR_cells.csv)
    - [mamm_FR_cells.csv](Data/mamm_FR_cells.csv)
    - [bat_FR_cells.csv](Data/bat_FR_cells.csv)
8. **XXXX_mean_fb.csv**: Mean functional beta diversity turnover of each grid cell. Separate files for birds, mammals, and bats.
    - [bird_mean_fb.csv](Data/bird_mean_fb.csv)
    - [mamm_mean_fb.csv](Data/mamm_mean_fb.csv)
    - [bat_mean_fb.csv](Data/bat_mean_fb.csv)
9. **XXXX_pb_NMDS.csv**: Phylobetadiversity turnover NMDS coordinates for 3 axes for each grid cell. Separate files for birds, mammals, and bats.
    - [bird_pb_NMDS.csv](Data/bird_pb_NMDS.csv)
    - [mamm_pb_NMDS.csv](Data/mamm_pb_NMDS.csv)
    - [bat_pb_NMDS.csv](Data/bat_pb_NMDS.csv)
10. **XXXX_pb_dist.RDS**: Distance matrix of phylobetadiversity turnover. Separate files for birds, mammals, and bats.
    - [bird_pb_dist.RDS](Data/bird_pb_dist.RDS)
    - [mamm_pb_dist.RDS](Data/mamm_pb_dist.RDS)
    - [bat_pb_dist.RDS](Data/bat_pb_dist.RDS)
11. **XXXX_func_PCoA.RDS**: Functional PCoA coordinates for each species. Separate files for birds, mammals, and bats.
    - [bird_func_PCoA.RDS](Data/bird_func_PCoA.RDS)
    - [mamm_func_PCoA.RDS](Data/mamm_func_PCoA.RDS)
    - [bat_func_PCoA.RDS](Data/bat_func_PCoA.RDS)
12. **XXXX_fb.RDS**: Output from functional beta diversity calculation, including distance matrices for total beta diversity, turnover, and nestedness. Separate files for mammals and bats.
    - mamm_fb.RDS: File too large to push to GitHub.
    - [bat_fb.RDS](Data/bat_fb.RDS)
