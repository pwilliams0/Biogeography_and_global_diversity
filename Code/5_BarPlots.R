#============================================#
# Bar plots of variance partitioning results #
#============================================#

library(tidyverse)

# ----- MAMMALS AND BIRDS (Fig. 2) ----------

# ---------- Species richness ----------
# Load data
results <- read.csv("Results/bird_SR_results.csv") %>%
  bind_rows(read.csv("Results/mamm_SR_results.csv")) %>%
  mutate(X = as.character(X),
         R2 = as.character(R2)) %>%
  filter(varpart != "Shared") %>%
  mutate(R2 = as.numeric(R2),
         varpart = factor(varpart,
                          levels = c("Phylobetadiversity only",
                                     "Environment only"))) 
# Save figure
svg("Results/barplot_SR.svg",
    width = 1.4, height = 1.4)
ggplot(results) +
  aes(fill=varpart, y=R2, x=Taxon) +
  geom_col(position=position_dodge(.75), colour="white",
           width=.75, linewidth=0.2) +
  theme_classic() +
  scale_fill_manual(values = c("#8E063B","#023FA5")) +
  scale_y_continuous(limits = c(0,.34), breaks = seq(0,.3,len = 7),
                     labels = c("0.0","","0.1","","0.2","","0.3"),
                     expand = expand_scale(add=c(0,.01))) +
  labs(y = "Variance explained") +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        line = element_line(colour = "black", linewidth = 0.2),
        axis.text = element_text(colour = "black", size=7),
        axis.title = element_text(colour = "black", size=7),
        axis.text.x = element_blank(),
        legend.position = "none") +
  guides(fill=guide_legend(ncol=1))
dev.off()

# ---------- Phylogenetic diversity ----------
# Load data
results <- read.csv("Results/bird_PD_results.csv") %>%
  bind_rows(read.csv("Results/mamm_PD_results.csv")) %>%
  mutate(X = as.character(X),
         R2 = as.character(R2)) %>%
  filter(varpart != "Shared") %>%
  mutate(R2 = as.numeric(R2),
         varpart = factor(varpart,
                          levels = c("Phylobetadiversity only",
                                     "Environment only"))) 
# Save figure
svg("Results/barplot_PD.svg",
    width = 1.4, height = 1.4)
ggplot(results) +
  aes(fill=varpart, y=R2, x=Taxon) +
  geom_col(position=position_dodge(.75), colour="white",
           width=.75, linewidth=0.2) +
  theme_classic() +
  scale_fill_manual(values = c("#8E063B","#023FA5")) +
  scale_y_continuous(limits = c(0,.34), breaks = seq(0,.3,len = 7),
                     labels = c("0.0","","0.1","","0.2","","0.3"),
                     expand = expand_scale(add=c(0,.01))) +
  labs(y = "Variance explained") +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        line = element_line(colour = "black", linewidth = 0.2),
        axis.text = element_text(colour = "black", size=7),
        axis.title = element_text(colour = "black", size=7),
        axis.text.x = element_blank(),
        legend.position = "none") +
  guides(fill=guide_legend(ncol=1))
dev.off()

# ---------- Functional richness ----------
# Load data
results <- read.csv("Results/bird_FR_results.csv") %>%
  bind_rows(read.csv("Results/mamm_FR_results.csv")) %>%
  mutate(X = as.character(X),
         R2 = as.character(R2)) %>%
  filter(varpart != "Shared") %>%
  mutate(R2 = as.numeric(R2),
         varpart = factor(varpart,
                          levels = c("Phylobetadiversity only",
                                     "Environment only"))) 
# Save figure
svg("Results/barplot_FR.svg",
    width = 1.4, height = 1.4)
ggplot(results) +
  aes(fill=varpart, y=R2, x=Taxon) +
  geom_col(position=position_dodge(.75), colour="white",
           width=.75, linewidth=0.2) +
  theme_classic() +
  scale_fill_manual(values = c("#8E063B","#023FA5")) +
  scale_y_continuous(limits = c(0,.34), breaks = seq(0,.3,len = 7),
                     labels = c("0.0","","0.1","","0.2","","0.3"),
                     expand = expand_scale(add=c(0,.01))) +
  labs(y = "Variance explained") +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        line = element_line(colour = "black", linewidth = 0.2),
        axis.text = element_text(colour = "black", size=7),
        axis.title = element_text(colour = "black", size=7),
        axis.text.x = element_blank(),
        legend.position = "none") +
  guides(fill=guide_legend(ncol=1))
dev.off()

# ---------- Mean func beta turnover ----------
# Load data
results <- read.csv("Results/bird_FB_results.csv") %>%
  bind_rows(read.csv("Results/mamm_FB_results.csv")) %>%
  mutate(X = as.character(X),
         R2 = as.character(R2)) %>%
  filter(varpart != "Shared") %>%
  mutate(R2 = as.numeric(R2),
         varpart = factor(varpart,
                          levels = c("Phylobetadiversity only",
                                     "Environment only"))) 
# Save figure
svg("Results/barplot_FB.svg",
    width = 1.4, height = 1.4)
ggplot(results) +
  aes(fill=varpart, y=R2, x=Taxon) +
  geom_col(position=position_dodge(.75), colour="white",
           width=.75, linewidth=0.2) +
  theme_classic() +
  scale_fill_manual(values = c("#8E063B","#023FA5")) +
  scale_y_continuous(limits = c(0,.34), breaks = seq(0,.3,len = 7),
                     labels = c("0.0","","0.1","","0.2","","0.3"),
                     expand = expand_scale(add=c(0,.01))) +
  labs(y = "Variance explained") +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        line = element_line(colour = "black", linewidth = 0.2),
        axis.text = element_text(colour = "black", size=7),
        axis.title = element_text(colour = "black", size=7),
        axis.text.x = element_blank(),
        legend.position = "none") +
  guides(fill=guide_legend(ncol=1))
dev.off()

# ----- BAT FUNCTIONAL RICHNESS (Fig. 6a) ----------
# Load data
results <- read.csv("Results/bat_FR_results.csv") %>%
  mutate(X = as.character(X),
         R2 = as.character(R2)) %>%
  filter(varpart != "Shared") %>%
  mutate(R2 = as.numeric(R2),
         varpart = factor(varpart,
                          levels = c("Phylobetadiversity only",
                                     "Environment only"))) 
# Save figure
svg("Results/barplot_bat_FR.svg",
    width = 1.2, height = 1.3)
ggplot(results) +
  aes(fill=varpart, y=R2, x=Taxon) +
  geom_col(position=position_dodge(.75), colour="white",
           width=.75, linewidth=0.2) +
  theme_classic() +
  scale_fill_manual(values = c("#8E063B","#023FA5")) +
  scale_y_continuous(limits = c(0,.255), breaks = seq(0,.25,len = 6),
                     labels = c("0.0","","0.1","","0.2",""),
                     expand = expand_scale(add=c(0,.01))) +
  labs(y = "Variance explained") +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        line = element_line(colour = "black", size = .25),
        axis.text = element_text(colour = "black", size=7),
        axis.text.x=element_blank(),
        axis.title = element_text(colour = "black", size=7),
        legend.position = "none") +
  guides(fill=guide_legend(ncol=1))
dev.off()

# ----- BAT OTHER RESULTS (Extended Data Fig. 7) ----------

# ---------- Species richness ----------
# Load data
results <- read.csv("Results/bat_SR_results.csv") %>%
  mutate(X = as.character(X),
         R2 = as.character(R2)) %>%
  filter(varpart != "Shared") %>%
  mutate(R2 = as.numeric(R2),
         varpart = factor(varpart,
                          levels = c("Phylobetadiversity only",
                                     "Environment only"))) 
# Save figure
svg("Results/barplot_bat_SR.svg",
    width = 1.2, height = 1.3)
ggplot(results) +
  aes(fill=varpart, y=R2, x=Taxon) +
  geom_col(position=position_dodge(.75), colour="white",
           width=.75, linewidth=0.2) +
  theme_classic() +
  scale_fill_manual(values = c("#8E063B","#023FA5")) +
  scale_y_continuous(limits = c(0,.57), breaks = seq(0,.5,len = 6),
                     expand = expand_scale(add=c(0,.01))) +
  labs(y = "Variance explained") +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        line = element_line(colour = "black", size = .25),
        axis.text = element_text(colour = "black", size=7),
        axis.title = element_text(colour = "black", size=7),
        axis.text.x = element_blank(),
        legend.position = "none") +
  guides(fill=guide_legend(ncol=1))
dev.off()

# ---------- Phylogenetic diversity ----------
# Load data
results <- read.csv("Results/bat_PD_results.csv") %>%
  mutate(X = as.character(X),
         R2 = as.character(R2)) %>%
  filter(varpart != "Shared") %>%
  mutate(R2 = as.numeric(R2),
         varpart = factor(varpart,
                          levels = c("Phylobetadiversity only",
                                     "Environment only"))) 
# Save figure
svg("Results/barplot_bat_PD.svg",
    width = 1.2, height = 1.3)
ggplot(results) +
  aes(fill=varpart, y=R2, x=Taxon) +
  geom_col(position=position_dodge(.75), colour="white",
           width=.75, linewidth=0.2) +
  theme_classic() +
  scale_fill_manual(values = c("#8E063B","#023FA5")) +
  scale_y_continuous(limits = c(0,.57), breaks = seq(0,.5,len = 6),
                     expand = expand_scale(add=c(0,.01))) +
  labs(y = "Variance explained") +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        line = element_line(colour = "black", size = .25),
        axis.text = element_text(colour = "black", size=7),
        axis.title = element_text(colour = "black", size=7),
        axis.text.x = element_blank(),
        legend.position = "none") +
  guides(fill=guide_legend(ncol=1))
dev.off()

# ---------- Mean func beta turnover ----------
# Load data
results <- read.csv("Results/bat_FB_results.csv") %>%
  mutate(X = as.character(X),
         R2 = as.character(R2)) %>%
  filter(varpart != "Shared") %>%
  mutate(R2 = as.numeric(R2),
         varpart = factor(varpart,
                          levels = c("Phylobetadiversity only",
                                     "Environment only"))) 
# Save figure
svg("Results/barplot_bat_FB.svg",
    width = 1.2, height = 1.3)
ggplot(results) +
  aes(fill=varpart, y=R2, x=Taxon) +
  geom_col(position=position_dodge(.75), colour="white",
           width=.75, linewidth=0.2) +
  theme_classic() +
  scale_fill_manual(values = c("#8E063B","#023FA5")) +
  scale_y_continuous(limits = c(0,.57), breaks = seq(0,.5,len = 6),
                     expand = expand_scale(add=c(0,.01))) +
  labs(y = "Variance explained") +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        line = element_line(colour = "black", size = .25),
        axis.text = element_text(colour = "black", size=7),
        axis.title = element_text(colour = "black", size=7),
        axis.text.x = element_blank(),
        legend.position = "none") +
  guides(fill=guide_legend(ncol=1))
dev.off()
