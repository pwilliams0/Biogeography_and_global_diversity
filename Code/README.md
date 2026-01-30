## [Code](Code)

1. **[0_DataPrep_AllTaxa.R](Code/0_DataPrep_AllTaxa.R)**: Create global 2x2 degree grid cells, calculate environment variables for each grid cell, and create assemblages for each grid cell. Environment variables include mean elevation, elevation range, present-day climate, past climate, distance between present-day climate and past climate, past ice cover, years since significant land conversion, and Human Impact Index. Assemblages are lists of all bird, mammal, or bat species with range maps that overlap a given cell.
2. **1_SR_PD_FR_XXXX.R**: Calculate species richness (SR), phylogenetic diversity (PD), and functional richness (FR) for each grid cell. Also calculate functional PCoA. Separate scripts for birds, mammals, and bats.
    - [1_SR_PD_FR_Bird.R](Code/1_SR_PD_FR_Bird.R)
    - [1_SR_PD_FR_Mamm.R](Code/1_SR_PD_FR_Mamm.R)
    - [1_SR_PD_FR_Bat.R](Code/1_SR_PD_FR_Bat.R)
3. **2_PhyloB_Bird_XX.R**: Calculate phylobetadiversity for birds on subsets of grid cells. Calculating all pairwise comparisions is computationally intractable, so grid cells are split into 5 groups.
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
4. **2_PhyloB_Mamm_XX.R**: Calculate phylobetadiversity for mammals on subsets of grid cells. Calculating all pairwise comparisions is computationally intractable, so grid cells are split into 4 groups.
    - [2_PhyloB_Mamm_12.R](Code/2_PhyloB_Mamm_12.R): Mammal phylobetadiversity for cells in groups 1 and 2.
    - [2_PhyloB_Mamm_13.R](Code/2_PhyloB_Mamm_13.R): Mammal phylobetadiversity for cells in groups 1 and 3.
    - [2_PhyloB_Mamm_14.R](Code/2_PhyloB_Mamm_14.R): Mammal phylobetadiversity for cells in groups 1 and 4.
    - [2_PhyloB_Mamm_23.R](Code/2_PhyloB_Mamm_23.R): Mammal phylobetadiversity for cells in groups 2 and 3.
    - [2_PhyloB_Mamm_24.R](Code/2_PhyloB_Mamm_24.R): Mammal phylobetadiversity for cells in groups 2 and 4.
    - [2_PhyloB_Mamm_34.R](Code/2_PhyloB_Mamm_34.R): Mammal phylobetadiversity for cells in groups 3 and 4.
5. **2_PhyloB_Bat_XX.R**: Calculate phylobetadiversity for bats on subsets of grid cells. Calculating all pairwise comparisions is computationally intractable, so grid cells are split into 3 groups.
    - [2_PhyloB_Bat_12.R](Code/2_PhyloB_Bat_12.R): Bat phylobetadiversity for cells in groups 1 and 2.
    - [2_PhyloB_Bat_13.R](Code/2_PhyloB_Bat_13.R): Bat phylobetadiversity for cells in groups 1 and 3.
    - [2_PhyloB_Bat_23.R](Code/2_PhyloB_Bat_23.R): Bat phylobetadiversity for cells in groups 2 and 3.
6. **2_FuncB_Bird_XX.R**: Calculate functional beta diversity for birds on subsets of grid cells. Calculating all pairwise comparisions is computationally intractable, so grid cells are split into 4 groups.
    - [2_FuncB_Bird_12.R](Code/2_FuncB_Bird_12.R): Bird functional beta diversity for cells in groups 1 and 2.
    - [2_FuncB_Bird_13.R](Code/2_FuncB_Bird_13.R): Bird functional beta diversity for cells in groups 1 and 3.
    - [2_FuncB_Bird_14.R](Code/2_FuncB_Bird_14.R): Bird functional beta diversity for cells in groups 1 and 4.
    - [2_FuncB_Bird_23.R](Code/2_FuncB_Bird_23.R): Bird functional beta diversity for cells in groups 2 and 3.
    - [2_FuncB_Bird_24.R](Code/2_FuncB_Bird_24.R): Bird functional beta diversity for cells in groups 2 and 4.
    - [2_FuncB_Bird_34.R](Code/2_FuncB_Bird_34.R): Bird functional beta diversity for cells in groups 3 and 4.
7. **3_PhyloB_XXXX_full.R**: Combine files to create full phylobetadiversity turnover distance matrix, and save NMDS coordinates for 3 axes. Separate scripts for birds, mammals, and bats.
    - [3_PhyloB_Bird_full.R](Code/3_PhyloB_Bird_full.R)
    - [3_PhyloB_Mamm_full.R](Code/3_PhyloB_Mamm_full.R)
    - [3_PhyloB_Bat_full.R](Code/3_PhyloB_Bat_full.R)
8. **[3_FuncB_Bird_full.R](Code/3_FuncB_Bird_full.R)**: Combine files to create functional beta diversity turnover distance matrix for birds, then calculate mean functional beta diversity turnover for birds.
9. **3_FuncB_XXXX.R**: Calculate functional beta diversity, then calculate mean functional beta diversity turnover. Separate scripts for mammals and bats.
    - [3_FuncB_Mamm.R](Code/3_FuncB_Mamm.R)
    - [3_FuncB_Bat.R](Code/3_FuncB_Bat.R)
10. **4_Analysis_XXXX.R**: Run models for species richness, phylogenetic diversity, functional richness, and mean functional beta diversity turnover. Save variance partitioning results. Calculate differences in residuals and save figures (Fig. 3 and Supplementary Fig. 1). Separate scripts for birds, mammals, and bats.
    - [4_Analysis_Bird.R](Code/4_Analysis_Bird.R)
    - [4_Analysis_Mamm.R](Code/4_Analysis_Mamm.R)
    - [4_Analysis_Bat.R](Code/4_Analysis_Bat.R)
11. **[5_BarPlots.R](Code/5_BarPlots.R)**: Create barplots of variance partitioning results for birds, mammals, and bats (Fig. 2, Fig. 6a, and Supplementary Fig. 5).
12. **[5_DiversityMaps.R](Code/5_DiversityMaps.R)**: Create maps of global diversity patterns of species richness, phylogenetic diversity, functional richness, and mean functional beta diversity turnover for birds, mammals, and bats (Fig. 4, Fig. 6b, and Supplementary Fig. 3).
13. **[5_PhyloB_Figures.R](Code/5_PhyloB_Figures.R)**: Create phylobetadiversity figures, including NMDS plots using 2 axes (Fig. 1), NMDS plots using 3 axes (Supplementary Fig. 4), and scree plots (Supplementary Fig. 7).
14. **[5_FuncSpace_Figures.R](Code/5_FuncSpace_Figures.R)**: Create figures showing overlap in functional space (Fig. 5 and Supplementary Fig. 6).
15. **5_Analysis_Realm_XXXX.R**: Run models for species richness, phylogenetic diversity, functional richness, and mean functional beta diversity turnover, using realm instead of phylobetadiversity.Calculate differences in residuals and save figures (Supplementary Fig. 2). Separate scripts for birds, mammals, and bats.
    - [5_Analysis_Realm_Bird.R](Code/5_Analysis_Realm_Bird.R)
    - [5_Analysis_Realm_Mamm.R](Code/5_Analysis_Realm_Mamm.R)
    - [5_Analysis_Realm_Bat.R](Code/5_Analysis_Realm_Bat.R)
16. **[5_PhyloSignal.R](Code/5_PhyloSignal.R)**: Calculate phylogenetic signal of traits (Supplementary Table 1).
17. **[5_AdditionalVariables.R](Code/5_AdditionalVariables.R)**: Calculate additional variance explained by environment variables besides climate (Supplementary Table 2).
18. **[5_PresentVsPastClimate.R](Code/5_PresentVsPastClimate.R)**: Calculate variance explained using climate data from present day, Last Glacial Maximum, and mid-Holocene (Supplementary Table 4).
