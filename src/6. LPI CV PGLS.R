#Supplementary Material - Phylogenetic Least Squares Regression
#Author(s): Reilly O'Connor
#Version: 2023-11-22   

#Pkgs
library(tidyverse)
library(RColorBrewer)
library(beepr)
library(ape)
library(caper)
library(easystats)
library(ggpubr)
library(extrafont)
library(visreg)
library(MuMIn)
library(brms)
library(phytools)
library(phylolm)

#load data
#Complete Phylogenetic Tree - From Phylacine (https://zenodo.org/records/3690867)
df_tree <- read.nexus("../Large Data Files/Phylacine v1.2 Data/Phylogenies/Complete_phylogeny.nex")

#LPI data & Mammal Groups
df_lpi <- read.csv("../Mammalian-Fast-Slow-Stability/Data/lpi mammal population cv detrended resolved taxonomy.csv", header = T)
df_mammal_groups <- read.csv("../Mammalian-Fast-Slow-Stability/Data/Suggested Mammal Groups.csv", header = T)

#Hatton Data
df_metab <- read.csv("../Mammalian-Fast-Slow-Stability/Data/slow_fast_metabolism resolved taxonomy.csv", header = T)
df_growth <- read.csv("../Mammalian-Fast-Slow-Stability/Data/slow_fast_growth resolved taxonomy.csv", header = T)
df_mortality <- read.csv("../Mammalian-Fast-Slow-Stability/Data/slow_fast_mortality resolved taxonomy.csv", header = T)

#Beccari et al 2024
df_mammal_fs <- read.csv("../Mammalian-Fast-Slow-Stability/Data/Beccari Fast-Slow resolved taxonomy.csv", header = T)

#Body Size Data for CV - See Mean Body Size Calculation.R Script
df_body_mass <- read.csv("../Mammalian-Fast-Slow-Stability/Data/mean species body size.csv")

##### Code #####
#Subset final species list
df_lpi_species_list <- df_lpi %>% 
  dplyr::select(Binomial, Order, Family, Genus, Species)

df_lpi_species_list <- unique(df_lpi_species_list)

#write.csv(df_lpi_species_list, "../Mammalian-Fast-Slow-Stability/Data/LPI CV Final Species List.csv")

df_lpi$Latitude <- as.numeric(df_lpi$Latitude)
df_lpi$Longitude <- as.numeric(df_lpi$Longitude)

df_lpi_groups <- merge(df_lpi, df_mammal_groups, by = "Binomial")

df_lpi_mammal_bs <- merge(df_lpi_groups, df_body_mass, by = "GBIF_ID") 

df_lpi_mbs <- df_lpi_mammal_bs %>%
  dplyr::select(Group, GBIF_ID, species_names, Species.y, Binomial, Common_name, cv_window, cv_value, mean_body_mass, Latitude, Longitude) %>%
  rename(Species = Species.y)

unique <- unique(df_lpi_mbs$Species)

df_lpi_fs <- merge(df_lpi_mbs, df_mammal_fs, by = "GBIF_ID")

df_mammal_cv_bs <- df_lpi_mbs %>%
  mutate(log_cv = log10(cv_value)) %>%
  group_by(Group, Binomial, GBIF_ID, cv_window) %>%
  reframe(mean_cv = mean(cv_value, na.rm = T),
          mean_log_cv = mean(log_cv, na.rm = T),
          mean_body_mass = mean(mean_body_mass, na.rm = T),
          mean_latitude = mean(Latitude, na.rm = T))

df_mammal_cv_fs <- df_lpi_fs %>%
  mutate(log_cv = log10(cv_value)) %>%
  group_by(Group, Binomial.x, GBIF_ID, cv_window) %>%
  reframe(mean_cv = mean(cv_value, na.rm = T),
          mean_log_cv = mean(log_cv, na.rm = T),
          mean_body_mass = mean(mean_body_mass, na.rm = T),
          mean_body_mass_fs = mean(bm, na.rm = T),
          mean_PC1 = mean(PC1, na.rm = T),
          mean_PC2 = mean(PC2, na.rm = T)) %>%
  rename(Binomial = Binomial.x)


#Set up x-axis (log body mass limits) & global color palette (by order)
#Set x-axis limits
x_min <- 1e+0
x_max <- 1e+9

#Set breaks for log10 scale (these are just example values)
xlog_breaks <- c(1e+01, 1e+02, 1e+03, 1e+04, 1e+05, 1e+06, 1e+07, 1e+08)

options(scipen = 0)

unique(df_mammal_cv_bs$Group)

##### 5 Year Window - Body Mass vs CV #####
df_cv_10 <- df_mammal_cv_bs %>% filter(cv_window == 5) #%>%
  #filter(! (Group == "Marine Mammals" | Group == "Bats"))

#Set order of mammal groups
df_cv_10$Group <- factor(df_cv_10$Group, levels = c("Insectivores", "Bats", "Glires", "Marsupials", "Primates", 
                                                    "Carnivores", "Ungulates", "Elephants", "Marine Mammals"))

lm_cv_bs_10 <- lm(mean_log_cv ~ log10(mean_body_mass), data = df_cv_10)
summary(lm_cv_bs_10)
par(mfrow = c(2,2))
plot(lm_cv_bs_10)

##### Calculate Residuals for Overlapping Data #####
#Merge dataframes and investigate residual correlations
df_metcv <- merge(df_cv_10, df_metab, by = c("GBIF_ID"))
non_overlapping_metab <- anti_join(df_cv_10, df_metab, by = c("GBIF_ID"))

df_grwcv <- merge(df_cv_10, df_growth, by = c("GBIF_ID"))
non_overlapping_growth <- anti_join(df_cv_10, df_growth, by = c("GBIF_ID"))

df_mortcv <- merge(df_cv_10, df_mortality, by = c("GBIF_ID"))
non_overlapping_metab <- anti_join(df_cv_10, df_mortality, by = c("GBIF_ID"))

options(scipen = 0)

#Metabolism
lm_metcv_bs <- lm(log10(Mass_Spec_Meta_Watt) ~ log10(Mass_g), data = df_metcv)
summary(lm_metcv_bs)
par(mfrow = c(2,2))
plot(lm_metcv_bs)

#Calculate Mass-Independent Specific Metabolism
df_metcv$residual_met <- residuals(lm_metcv_bs)

lm_cvmet_bs <- lm(mean_log_cv ~ log10(Mass_g), data = df_metcv)
summary(lm_cvmet_bs)
par(mfrow = c(2,2))
plot(lm_cvmet_bs)

#calculate residual CV - just for quick look
df_metcv$residual_cv <- residuals(lm_cvmet_bs)

#Quick look at residual relationships
lm_resid_met_cv <- lm(data = df_metcv, residual_cv ~ residual_met)
summary(lm_resid_met_cv)
par(mfrow = c(2,2))
plot(lm_resid_met_cv)

#Growth
lm_grwcv <- lm(log10(rmax) ~ log10(Mass_g), data = df_grwcv)
summary(lm_grwcv)
par(mfrow = c(2,2))
plot(lm_grwcv)

#Calculate Mass-Independent rmax
df_grwcv$residual_rmax <- residuals(lm_grwcv)

lm_cvgrw_bs <- lm(mean_log_cv ~ log10(Mass_g), data = df_grwcv)
summary(lm_cvgrw_bs)
par(mfrow = c(2,2))
plot(lm_cvgrw_bs)

#calculate residual CV - just for quick look
df_grwcv$residual_cv <- residuals(lm_cvgrw_bs)

lm_resid_grw_cv <- lm(data = df_grwcv, residual_cv ~ residual_rmax)
summary(lm_resid_grw_cv)
par(mfrow = c(2,2))
plot(lm_resid_grw_cv)

#Mortality
lm_mortcv <- lm(log10(Mortality_per_yr) ~ log10(Mass_g), data = df_mortcv)
summary(lm_mortcv)
par(mfrow = c(2,2))
plot(lm_mortcv)

#Calculate Mass-Independent Mortality
df_mortcv$residual_mort <- residuals(lm_mortcv)

lm_cvmort_bs <- lm(mean_log_cv ~ log10(Mass_g), data = df_mortcv)
summary(lm_cvmort_bs)
par(mfrow = c(2,2))
plot(lm_cvmort_bs)

#calculate residual CV - just for quick look
df_mortcv$residual_cv <- residuals(lm_cvmort_bs)

lm_resid_mort_cv <- lm(data = df_mortcv, residual_cv ~ residual_mort)
summary(lm_resid_mort_cv)
check_model(lm_resid_mort_cv)


##### Phylogenetic Generalized Least Squares Regression #####
#Ensure the trait data is ordered to match the tree tip labels
#CV ~ body size
tree_tips <- df_tree[[1]]$tip.label
df_cv_10_tree <- df_cv_10[match(tree_tips, df_cv_10$Binomial), ]
df_cv_10_tree_na <- na.omit(df_cv_10_tree)
df_cv_10_tree_na <- as.data.frame(df_cv_10_tree_na)
rownames(df_cv_10_tree_na) <- df_cv_10_tree_na$Binomial

phylolm_cv_model <- phylolm(formula = mean_log_cv ~ log10(mean_body_mass),
                            data = df_cv_10_tree_na, 
                            phy = df_tree[[1]], 
                            model = "lambda")

summary(phylolm_cv_model)

#CV ~ body size + metabolism
tree_tips <- df_tree[[1]]$tip.label
df_cv_metab_10_tree <- df_metcv[match(tree_tips, df_metcv$Binomial), ]
df_cv_metab_10_tree_na <- na.omit(df_cv_metab_10_tree)
df_cv_metab_10_tree_na <- as.data.frame(df_cv_metab_10_tree_na)
rownames(df_cv_metab_10_tree_na) <- df_cv_metab_10_tree_na$Binomial

phylolm_metab_model <- phylolm(formula = mean_log_cv ~ log10(mean_body_mass) + residual_met,
                               data = df_cv_metab_10_tree_na, 
                               phy = df_tree[[1]], 
                               model = "lambda")

# iew the summary of the model
summary(phylolm_metab_model)

#CV ~ body size + rmax
tree_tips <- df_tree[[1]]$tip.label
df_cv_grw_10_tree <- df_grwcv[match(tree_tips, df_grwcv$Binomial), ]
df_cv_grw_10_tree_na <- na.omit(df_cv_grw_10_tree)
df_cv_grw_10_tree_na <- as.data.frame(df_cv_grw_10_tree_na)
rownames(df_cv_grw_10_tree_na) <- df_cv_grw_10_tree_na$Binomial

phylolm_grw_model <- phylolm(formula = mean_log_cv ~ log10(mean_body_mass) + residual_rmax,
                             data = df_cv_grw_10_tree_na, 
                             phy = df_tree[[1]], 
                             model = "lambda")

# View the summary of the model
summary(phylolm_grw_model)


#CV ~ body size + mortality
tree_tips <- df_tree[[1]]$tip.label
df_cv_mort_10_tree <- df_mortcv[match(tree_tips, df_mortcv$Binomial), ]
df_cv_mort_10_tree_na <- na.omit(df_cv_mort_10_tree)
df_cv_mort_10_tree_na <- as.data.frame(df_cv_mort_10_tree_na)
rownames(df_cv_mort_10_tree_na) <- df_cv_mort_10_tree_na$Binomial

phylolm_mort_model <- phylolm(formula = mean_log_cv ~ log10(mean_body_mass) + residual_mort,
                              data = df_cv_mort_10_tree_na, 
                              phy = df_tree[[1]], 
                              model = "lambda")

# View the summary of the model
summary(phylolm_mort_model)

#PC1
df_mammal_cv_fs_5 <- df_mammal_cv_fs %>% filter(cv_window == 5)# %>%
  #filter(! (Group == "Bats"| Group == "Marine Mammals"))

lm_cv_pc1 <- lm(mean_log_cv ~ log10(mean_body_mass) + mean_PC1, data = df_mammal_cv_fs_5)
summary(lm_cv_pc1)
par(mfrow = c(2,2))
plot(lm_cv_pc1)

#CV ~ body size + LH Traits Beccari et al 2024
tree_tips <- df_tree[[1]]$tip.label
df_cv_fs_10_tree <- df_mammal_cv_fs_5[match(tree_tips, df_mammal_cv_fs_5$Binomial), ]
df_cv_fs_10_tree_na <- na.omit(df_cv_fs_10_tree)
df_cv_fs_10_tree_na <- as.data.frame(df_cv_fs_10_tree_na)
rownames(df_cv_fs_10_tree_na) <- df_cv_fs_10_tree_na$Binomial

phylolm_fs_model <- phylolm(formula = mean_log_cv ~ log10(mean_body_mass) + mean_PC1,
                               data = df_cv_fs_10_tree_na, 
                               phy = df_tree[[1]], 
                               model = "lambda")

# iew the summary of the model
summary(phylolm_fs_model)


