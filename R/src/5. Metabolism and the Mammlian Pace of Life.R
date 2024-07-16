# Exploring Relationships between Life History Speed, 
# Fast-Slow Life History Strategies and Metabolism in Mammals

#Author(s): Reilly O'Connor
#Version: 2024-06-11

#Pkgs
library(tidyverse)
library(RColorBrewer)
library(Hmisc)

#load all data

#Hatton Data
df_metab <- read.csv("../Mammalian-Fast-Slow-Stability/Data/slow_fast_metabolism resolved taxonomy.csv", header = T)
df_growth <- read.csv("../Mammalian-Fast-Slow-Stability/Data/slow_fast_growth resolved taxonomy.csv", header = T)

#Beccari et al 2024
df_mammal_fs <- read.csv("../Mammalian-Fast-Slow-Stability/Data/Beccari Fast-Slow resolved taxonomy.csv", header = T)

##### Code #####

#Combine species across dataframes
#calculate residuals for rmax
lm_rmax <- lm(log10(rmax) ~ log10(Mass_g), data = df_growth)
summary(lm_rmax)
df_growth$residual_rmax_all <- residuals(lm_rmax)

df_rmax_fs <- merge(df_mammal_fs, df_growth, by = c("GBIF_ID"))

#Average mass across dataframes, center and scale PC1
df_rmax_fs_bs <- df_rmax_fs %>%
  mutate(Mass_g = (Mass_g + bm)/2,
         PC1 = as.numeric(scale(PC1))) %>%
  dplyr::select(GBIF_ID, Species.y, Mass_g, rmax, residual_rmax_all, PC1, PC2) %>%
  rename(Species = Species.y)

lm_rmax_fs <- lm(log10(rmax) ~ log10(Mass_g), data = df_rmax_fs_bs)
summary(lm_rmax_fs)
df_rmax_fs_bs$residual_rmax <- residuals(lm_rmax_fs)

rcorr(df_rmax_fs_bs$PC1, df_rmax_fs_bs$residual_rmax, type = "pearson")

lm_rmax_fs <- lm(log10(rmax) ~ log10(Mass_g) + PC1, data = df_rmax_fs_bs)
summary(lm_rmax_fs)

#Calculate residual metabolic rate for overlapping data
lm_mmammal_metab <- lm(log10(Mass_Spec_Meta_Watt) ~ log10(Mass_g), data = df_metab)
summary(lm_mmammal_metab)
df_metab$residual_met <- residuals(lm_mmammal_metab)

df_rmax_metab <- merge(df_metab, df_growth, by = c("GBIF_ID"))

#Average mass across dataframes, center and scale residual metabolism
df_rmax_metab_bs <- df_rmax_metab %>%
  mutate(Mass_g = (Mass_g.x + Mass_g.y)/2,
         residual_met = scale(residual_met)) %>%
  dplyr::select(GBIF_ID, Species.y, Mass_g, residual_met, Metabolism_W, Mass_Spec_Meta_Watt, rmax) %>%
  rename(Species = Species.y)

df_fs_metab <- merge(df_metab, df_mammal_fs, by = c("GBIF_ID"))

#Average mass across dataframes, center and scale PC1 & residual metabolism
df_fs_metab_bs <- df_fs_metab %>%
  mutate(Mass_g = (Mass_g + bm)/2,
         residual_met = as.numeric(scale(residual_met)),
         ls = as.numeric(scale(ls)),
         ly = as.numeric(scale(ly)),
         long = as.numeric(scale(long)),
         gest = as.numeric(scale(gest)),
         wea = as.numeric(scale(wea)),
         fmat = as.numeric(scale(fmat)),
         PC1 = as.numeric(scale(PC1))) %>%
  rename(Species = Species.y)

lm_rmax_bs <- lm(log10(rmax) ~ log10(Mass_g), data = df_rmax_metab_bs)
summary(lm_rmax_bs)

df_rmax_metab_bs$residual_rmax <- residuals(lm_rmax_bs)

lm_rmax_fs <- lm(log10(rmax) ~ log10(Mass_g) + residual_met, data = df_rmax_metab_bs)
summary(lm_rmax_fs)

rcorr(df_rmax_metab_bs$residual_rmax, df_rmax_metab_bs$residual_met, type = "pearson")

lm_rmax_fs <- lm(PC1 ~ residual_met, data = df_fs_metab_bs)
summary(lm_rmax_fs)

rcorr(df_fs_metab_bs$PC1, df_fs_metab_bs$residual_met, type = "pearson")

df_rmax_fs_met <- merge(df_rmax_fs, df_metab, by = c("GBIF_ID"))

df_rmax_fs_met_bs <- df_rmax_fs_met %>%
  mutate(Mass_g = (Mass_g.x + Mass_g.y + bm)/3,
         residual_met = as.numeric(scale(residual_met)),
         ls = as.numeric(scale(ls)),
         ly = as.numeric(scale(ly)),
         long = as.numeric(scale(long)),
         gest = as.numeric(scale(gest)),
         wea = as.numeric(scale(wea)),
         fmat = as.numeric(scale(fmat)),
         PC1 = as.numeric(scale(PC1)),
         residual_rmax_all = scale(residual_rmax_all))

lm_rmax_fs_met <- lm(log10(rmax) ~ log10(Mass_g) + PC1 + residual_met, data = df_rmax_fs_met_bs)
summary(lm_rmax_fs_met)

rcorr(df_rmax_fs_met_bs$PC1, df_rmax_fs_met_bs$residual_rmax_all, type = "pearson")
rcorr(df_rmax_fs_met_bs$PC1, df_rmax_fs_met_bs$residual_met, type = "pearson")
