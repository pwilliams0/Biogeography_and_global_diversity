#=================================================================#
# Calculate phylogenetic signal of traits (Extended Data Table 1) #
#=================================================================#

library(tidyverse)
library(ape)
library(phytools)

# ----- MAMMALS ----------

# Load trait data
mamm_traits <- read.csv("Data/mamm_cell_df.csv") %>%
  group_by(names_phylo) %>% slice(1) %>% ungroup() %>%
  dplyr::select(names_phylo, Diet.Inv, Diet.Vert,
                Diet.Fish, Diet.Scav, Diet.Fruit,
                Diet.Nect, Diet.Seed, Diet.Herb,
                ForStrat, LogBodyMass) %>%
  mutate(ForStrat = as.factor(ForStrat)) %>%
  column_to_rownames("names_phylo")

# Load phylogeny
mamm_phylo <- ape::read.tree("Data/Raw/Upham_mean_phylogeny_ultrametric.tree")

# Extract each trait
Diet.Inv <- setNames(mamm_traits$Diet.Inv,
                      rownames(mamm_traits))
Diet.Vert <- setNames(mamm_traits$Diet.Vert,
                   rownames(mamm_traits))
Diet.Fish <- setNames(mamm_traits$Diet.Fish,
                    rownames(mamm_traits))
Diet.Scav <- setNames(mamm_traits$Diet.Scav,
                    rownames(mamm_traits))
Diet.Fruit <- setNames(mamm_traits$Diet.Fruit,
                    rownames(mamm_traits))
Diet.Nect <- setNames(mamm_traits$Diet.Nect,
                    rownames(mamm_traits))
Diet.Seed <- setNames(mamm_traits$Diet.Seed,
                    rownames(mamm_traits))
Diet.Herb <- setNames(mamm_traits$Diet.Herb,
                    rownames(mamm_traits))
LogBodyMass <- setNames(mamm_traits$LogBodyMass,
                      rownames(mamm_traits))
# Calculate phylogenetic signal K
K.Diet.Fish <- phylosig(mamm_phylo, Diet.Fish,
                   test=TRUE)
print(K.Diet.Fish)
K.Diet.Fruit <- phylosig(mamm_phylo, Diet.Fruit,
                   test=TRUE)
print(K.Diet.Fruit)
K.Diet.Herb <- phylosig(mamm_phylo, Diet.Herb,
                   test=TRUE)
print(K.Diet.Herb)
K.Diet.Inv <- phylosig(mamm_phylo, Diet.Inv,
                   test=TRUE)
print(K.Diet.Inv)
K.Diet.Nect <- phylosig(mamm_phylo, Diet.Nect,
                   test=TRUE)
print(K.Diet.Nect)
K.Diet.Scav <- phylosig(mamm_phylo, Diet.Scav,
                   test=TRUE)
print(K.Diet.Scav)
K.Diet.Seed <- phylosig(mamm_phylo, Diet.Seed,
                   test=TRUE)
print(K.Diet.Seed)
K.Diet.Vert <- phylosig(mamm_phylo, Diet.Vert,
                   test=TRUE)
print(K.Diet.Vert)
K.LogBodyMass <- phylosig(mamm_phylo,LogBodyMass,
                     test=TRUE)
print(K.LogBodyMass)

# ----- BIRDS ----------

# Load trait data
bird_traits <- read.csv("Data/bird_cell_df.csv") %>%
  group_by(names_phylo) %>% slice(1) %>% ungroup() %>%
  dplyr::select(names_phylo, Diet.Inv, Diet.Vert,
                Diet.Fish, Diet.Scav, Diet.Fruit,
                Diet.Nect, Diet.Seed, Diet.Herb,
                ForStrat.aquatic,
                ForStrat.ground,
                ForStrat.understory,
                ForStrat.canopy,
                ForStrat.aerial,
                LogBodyMass) %>%
  column_to_rownames("names_phylo")

# Load phylogeny
bird_phylo <- ape::read.tree("Data/Raw/Stage2_MayrParSho_Hackett_mean phylogeny_ultrametric.tre")

# Extract each trait
Diet.Inv <- setNames(bird_traits$Diet.Inv,
                     rownames(bird_traits))
Diet.Vert <- setNames(bird_traits$Diet.Vert,
                      rownames(bird_traits))
Diet.Fish <- setNames(bird_traits$Diet.Fish,
                      rownames(bird_traits))
Diet.Scav <- setNames(bird_traits$Diet.Scav,
                      rownames(bird_traits))
Diet.Fruit <- setNames(bird_traits$Diet.Fruit,
                       rownames(bird_traits))
Diet.Nect <- setNames(bird_traits$Diet.Nect,
                      rownames(bird_traits))
Diet.Seed <- setNames(bird_traits$Diet.Seed,
                      rownames(bird_traits))
Diet.Herb <- setNames(bird_traits$Diet.Herb,
                      rownames(bird_traits))
ForStrat.aquatic <- setNames(bird_traits$ForStrat.aquatic,
                            rownames(bird_traits))
ForStrat.ground <- setNames(bird_traits$ForStrat.ground,
                            rownames(bird_traits))    
ForStrat.understory <- setNames(bird_traits$ForStrat.understory,
                                rownames(bird_traits))
ForStrat.canopy <- setNames(bird_traits$ForStrat.canopy,
                            rownames(bird_traits))    
ForStrat.aerial <- setNames(bird_traits$ForStrat.aerial,
                            rownames(bird_traits))
LogBodyMass <- setNames(bird_traits$LogBodyMass,
                        rownames(bird_traits))

# Calculate phylogenetic signal K
K.Diet.Fish <- phylosig(bird_phylo, Diet.Fish,
                        test=TRUE)
print(K.Diet.Fish)
K.Diet.Fruit <- phylosig(bird_phylo, Diet.Fruit,
                         test=TRUE)
print(K.Diet.Fruit)
K.Diet.Herb <- phylosig(bird_phylo, Diet.Herb,
                        test=TRUE)
print(K.Diet.Herb)
K.Diet.Inv <- phylosig(bird_phylo, Diet.Inv,
                       test=TRUE)
print(K.Diet.Inv)
K.Diet.Nect <- phylosig(bird_phylo, Diet.Nect,
                        test=TRUE)
print(K.Diet.Nect)
K.Diet.Scav <- phylosig(bird_phylo, Diet.Scav,
                        test=TRUE)
print(K.Diet.Scav)
K.Diet.Seed <- phylosig(bird_phylo, Diet.Seed,
                        test=TRUE)
print(K.Diet.Seed)
K.Diet.Vert <- phylosig(bird_phylo, Diet.Vert,
                        test=TRUE)
print(K.Diet.Vert)
K.ForStrat.aerial <- phylosig(bird_phylo, ForStrat.aerial,
                              test=TRUE)
print(K.ForStrat.aerial)
K.ForStrat.aquatic <- phylosig(bird_phylo, ForStrat.aquatic,
                               test=TRUE)
print(K.ForStrat.aquatic)
K.ForStrat.canopy <- phylosig(bird_phylo, ForStrat.canopy,
                              test=TRUE)
print(K.ForStrat.canopy)
K.ForStrat.ground <- phylosig(bird_phylo, ForStrat.ground,
                              test=TRUE)
print(K.ForStrat.ground)
K.ForStrat.understory <- phylosig(bird_phylo, ForStrat.understory,
                                  test=TRUE)
print(K.ForStrat.understory)
K.LogBodyMass <- phylosig(bird_phylo,LogBodyMass,
                          test=TRUE)
print(K.LogBodyMass)
