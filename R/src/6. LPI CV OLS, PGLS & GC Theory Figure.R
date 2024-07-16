 #Main Analysis - Slow-Fast Body Size CV & Residual Relationships

#Author(s): Reilly O'Connor
#Version: 2024-01-24

#Pkgs
library(RColorBrewer)
library(beepr)
library(easystats)
library(ggpubr)
library(extrafont)
library(visreg)
library(MuMIn)
library(tidyverse)
library(ape)
library(caper)
library(phylolm)

#load data
df_lpi <- read.csv("../Mammalian-Fast-Slow-Stability/Data/lpi mammal population cv detrended resolved taxonomy.csv", header = T)

df_mammal_groups <- read.csv("../Mammalian-Fast-Slow-Stability/Data/Suggested Mammal Groups.csv", header = T)

#Hatton Data
df_metab <- read.csv("../Mammalian-Fast-Slow-Stability/Data/slow_fast_metabolism resolved taxonomy.csv", header = T)

df_growth <- read.csv("../Mammalian-Fast-Slow-Stability/Data/slow_fast_growth resolved taxonomy.csv", header = T)

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

unique(df_lpi_mammal_bs$ID)

df_lpi_mbs <- df_lpi_mammal_bs %>%
  dplyr::select(Group, GBIF_ID, species_names, Species.y, Binomial, Common_name, cv_window, cv_value, mean_body_mass, Latitude, Longitude) %>%
  rename(Species = Species.y)


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

##### Set colors for Groups #####
#Number of mammal groups
num_groups <- 9

#Generate color palette
palette <- c("#1E77B4", "#D5BD00", "#3CA373",
             "#62B400", "#8156B6", "#EB8F00",
             "#D13120", "#D90D8E", "#999999")

#Assign colors to mammal groups
groups <- c("Glires", "Carnivores", "Bats", 
            "Ungulates", "Primates", "Marine Mammals", 
            "Insectivores", "Marsupials", "Elephants")

color_mapping <- setNames(palette, groups)

##### 5 Year Window - Body Mass vs CV #####
df_cv_5 <- df_mammal_cv_bs %>% filter(cv_window == 5) #%>%
 #filter(! (Group == "Marine Mammals" | Group == "Bats"))

#Set order of mammal groups
df_cv_5$Group <- factor(df_cv_5$Group, levels = c("Insectivores", "Bats", "Glires", "Primates", "Marsupials", "Carnivores", "Ungulates", "Elephants", "Marine Mammals"))

lm_cv_bs_5 <- lm(mean_log_cv ~ log10(mean_body_mass), data = df_cv_5)
summary(lm_cv_bs_5)
par(mfrow = c(2,2))
plot(lm_cv_bs_5)


gg_cv_bs_5 <- ggplot(df_cv_5, aes(x = mean_body_mass, y = mean_cv, color = Group)) +
  geom_point(size = 1.5) +
  geom_smooth(se = T, method = "lm", aes(color = NA), color = "black", linetype = "solid", linewidth = 1.0) +
  ylab("Population Stability (CV)") +
  xlab("Body Mass (g)") +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks = element_line(linetype = "solid"),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y =element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        text = element_text(family = "Arial")) +
  scale_color_manual(values = color_mapping) +
  scale_x_log10(limits = c(x_min, x_max), breaks = xlog_breaks) +
  scale_y_log10()

gg_cv_bs_5

#ggsave("../Mammalian-Fast-Slow-Stability/Figures/Figure S2 - Body Mass vs CV Terrestrial Mammals.jpeg", plot = gg_cv_bs_5, width = 8, height = 6)

gg_mean_cv_groups <- ggplot(df_cv_5, aes(x = Group, y = mean_cv, fill = Group)) +
  geom_boxplot(outliers = F) +
  geom_jitter(aes(color = Group), size = 2, shape = 21, stroke = 0.25, color = "black") +
  ylab("Population Stability (CV)") +
  xlab("Group") +
theme(axis.line = element_line(linetype = "solid"),
      axis.ticks = element_line(linetype = "blank"),
      panel.background = element_rect(fill = NA),
      legend.key = element_rect(fill = NA),
      legend.background = element_rect(fill = NA),
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y =element_text(size = 12), 
      axis.title.x = element_text(size = 12), 
      text = element_text(family = "Arial")) +
  scale_fill_manual(values = color_mapping) +
  scale_y_log10()

gg_mean_cv_groups

aov_cv <- aov(mean_log_cv ~ Group, data = df_cv_5)
summary(aov_cv)
TukeyHSD(aov_cv)

options(scipen = 9)
gg_cv_bs_order_mean <- df_cv_5 %>%
  group_by(Group) %>%
  reframe(ave_cv = mean(mean_cv),
          sd_cv = sd(mean_cv),
          ave_body_mass = mean(mean_body_mass),
          sd_body_mass = sd(mean_body_mass),
          count = n(),
          se_cv = (sd_cv/(sqrt(count))),
          se_body_mass = (sd_body_mass/(sqrt(count)))) %>%
  arrange(ave_body_mass)
options(scipen = 0)

#ggsave("../Mammalian-Fast-Slow-Stability/Figures/Figure S3 - Mean CV Among Groups.jpeg", plot = gg_mean_cv_groups, width = 12, height = 8)


##### Calculate Residuals for Overlapping Data #####

#Merge dataframes and investigate residual correlations
df_grwcv <- merge(df_cv_5, df_growth, by = c("GBIF_ID"))
df_metcv <- merge(df_cv_5, df_metab, by = c("GBIF_ID"))

options(scipen = 0)

#rmax
lm_grwcv <- lm(log10(rmax) ~ log10(Mass_g), data = df_grwcv)
summary(lm_grwcv)
par(mfrow = c(2,2))
plot(lm_grwcv)

#Calculate Mass-Independent rmax
df_grwcv$residual_rmax <- residuals(lm_grwcv)

gg_grwcv_bs <- ggplot(df_grwcv, aes(x = Mass_g, y = rmax, color = Group.x)) +
  geom_point(size = 1.5) +
  geom_smooth(se = T, method = "lm", aes(color = NA), color = "black", linetype = "solid", linewidth = 0.75) +
  ylab("Growth Rate (per year)") +
  xlab("Body Mass (g)") +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks = element_line(linetype = "solid"),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y =element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        text = element_text(family = "Arial")) +
  scale_x_log10(limits = c(x_min, x_max), breaks = xlog_breaks) +
  scale_color_manual(values = color_mapping) +
  scale_y_log10()

gg_grwcv_bs

lm_cvgrw_bs <- lm(mean_log_cv ~ log10(Mass_g), data = df_grwcv)
summary(lm_cvgrw_bs)
par(mfrow = c(2,2))
plot(lm_cvgrw_bs)

#calculate residual CV - just for quick look
df_grwcv$residual_cv <- residuals(lm_cvgrw_bs)

gg_cvgrw_bs <- ggplot(df_grwcv, aes(x = Mass_g, y = mean_log_cv, color = Group.x)) +
  geom_point(size = 1.5) +
  geom_smooth(se = T, method = "lm", aes(color = NA), color = "black", linetype = "solid", linewidth = 0.75) +
  ylab("Population Variability (CV)") +
  xlab("Body Mass (g)") +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks = element_line(linetype = "solid"),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y =element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        text = element_text(family = "Arial")) +
  scale_color_manual(values = color_mapping) +
  scale_x_log10(limits = c(x_min, x_max), breaks = xlog_breaks)

gg_cvgrw_bs


lm_resid_grw_cv <- lm(data = df_grwcv, scale(residual_cv) ~ scale(residual_rmax))
summary(lm_resid_grw_cv)
par(mfrow = c(2,2))
plot(lm_resid_grw_cv)

gg_grwcv_resid <- ggplot(df_grwcv, aes(x = scale(residual_rmax), y = scale(residual_cv), color = Group.x)) +
  geom_point(size = 1.5) +
  geom_smooth(se = T, method = "lm", aes(x = scale(residual_rmax), y = scale(residual_cv)), color = "black",  linetype = "solid", linewidth = 0.8) +
  ylab("Residual CV") +
  xlab("Residual rmax") +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks = element_line(linetype = "solid"),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y =element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        text = element_text(family = "Arial")) +
  scale_color_manual(values = color_mapping)

gg_grwcv_resid

#Metabolism
lm_metcv_bs <- lm(log10(Mass_Spec_Meta_Watt) ~ log10(Mass_g), data = df_metcv)
summary(lm_metcv_bs)
par(mfrow = c(2,2))
plot(lm_metcv_bs)

#Calculate Mass-Independent Specific Metabolism
df_metcv$residual_met <- residuals(lm_metcv_bs)

gg_metabcv_bs <- ggplot(df_metcv, aes(x = Mass_g, y = Mass_Spec_Meta_Watt, color = Group.x)) +
  geom_point(size = 1.5) +
  geom_smooth(se = T, method = "lm", aes(color = NA), color = "black", linetype = "solid", linewidth = 0.75) +
  ylab("Mass Specific Metabolic Rate (Watts/g)") +
  xlab("Body Mass (g)") +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks = element_line(linetype = "solid"),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y =element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        text = element_text(family = "Arial")) +
  scale_color_manual(values = color_mapping) +
  scale_x_log10(limits = c(x_min, x_max), breaks = xlog_breaks) +
  scale_y_log10()

gg_metabcv_bs

lm_cvmet_bs <- lm(mean_log_cv ~ log10(Mass_g), data = df_metcv)
summary(lm_cvmet_bs)
par(mfrow = c(2,2))
plot(lm_cvmet_bs)

#calculate residual CV - just for quick look
df_metcv$residual_cv <- residuals(lm_cvmet_bs)

gg_cvmet_bs <- ggplot(df_metcv, aes(x = Mass_g, y = mean_log_cv, color = Group.x)) +
  geom_point(size = 1.5) +
  geom_smooth(se = T, method = "lm", aes(color = NA), color = "black", linetype = "solid", linewidth = 0.75) +
  ylab("Population Variability (CV)") +
  xlab("Body Mass (g)") +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks = element_line(linetype = "solid"),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y =element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        text = element_text(family = "Arial")) +
  scale_color_manual(values = color_mapping) +
  scale_x_log10(limits = c(x_min, x_max), breaks = xlog_breaks)

gg_cvmet_bs

#Quick look at residual relationships
lm_resid_met_cv <- lm(data = df_metcv, scale(residual_cv) ~ scale(residual_met))
summary(lm_resid_met_cv)
par(mfrow = c(2,2))
plot(lm_resid_met_cv)

gg_metcv_resid <- ggplot(df_metcv, aes(x = scale(residual_met), y = scale(residual_cv), color = Group.x)) +
  geom_point(size = 1.5) +
  geom_smooth(se = T, method = "lm", aes(x = scale(residual_met), y = scale(residual_cv)), color = "black",  linetype = "solid", linewidth = 0.8) +
  ylab("Residual CV") +
  xlab("Residual Specific Metabolism") +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks = element_line(linetype = "solid"),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y =element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        text = element_text(family = "Arial")) +
  scale_color_manual(values = color_mapping)

gg_metcv_resid

##### Life History Trait Relationships - Overlapping Data Panel #####
gg_met_grw_mort_cv_ol_3 <- ggarrange(gg_grwcv_bs, gg_metabcv_bs, 
                                     nrow = 2, ncol = 1,
                                     labels = c("A", "B"))

gg_met_grw_mort_cv_ol_3

#ggsave("../Mammalian-Fast-Slow-Stability/Figures/Figure S1 - LH Trait BS Overlap.jpeg", plot = gg_met_grw_mort_cv_ol_3, width = 6, height = 8)


##### Life History Trait - CV Body Size for Overlapping Data Panel #####
gg_cv_met_grw_mort_3 <- ggarrange(gg_cvmet_bs, gg_cvgrw_bs,
                                  nrow = 1, ncol = 2,
                                  labels = c("A", "B"))

gg_cv_met_grw_mort_3

#ggsave("../Mammalian-Fast-Slow-Stability/Figures/Figure SX - CV BS Overlap.jpeg", plot = gg_cv_met_grw_mort_3, width = 16, height = 4)


##### Multiple Linear Regressions with Residual Life History Traits #####\

##### MLM rmax #####
#scale residuals and ensure numeric vector for plotting
df_grwcv <- df_grwcv %>%
  mutate(residual_rmax = scale(residual_rmax))
df_grwcv$residual_rmax <- as.numeric(df_grwcv$residual_rmax)

mlm_grwcv <- lm(mean_log_cv ~ log10(mean_body_mass) + residual_rmax, data = df_grwcv)
summary(mlm_grwcv)
par(mfrow = c(2,2))
plot(mlm_grwcv)

gg_partial_grwcvbs <- visreg(mlm_grwcv, "mean_body_mass", overlay=FALSE, partial=TRUE, gg=TRUE, points=list(col="NA")) + 
  geom_point(data = df_grwcv, size = 1.5, aes(x = mean_body_mass, y = mean_log_cv, color = Group.x)) +
  geom_smooth(se = F, color = "black", method = "lm") +
  theme_classic() +
  ylab("Population Stability (CV)") +
  xlab("Body Size (g)") +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks = element_line(linetype = "solid"),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y =element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        text = element_text(family = "Arial")) +
  scale_x_log10(limits = c(x_min, x_max), breaks = xlog_breaks)   +
  scale_color_manual(values = color_mapping) +
  ylim(-2.0, 0.5)

gg_partial_grwcvbs

gg_partial_grwcv <- visreg(mlm_grwcv, "residual_rmax", overlay=FALSE, partial=TRUE, gg=TRUE, points=list(col="NA")) + 
  geom_point(data = df_grwcv, size = 1.5, aes(x = residual_rmax, y = mean_log_cv, color = Group.x)) +
  geom_smooth(se = F, color = "black", method = "lm") +
  theme_classic() +
  ylab("Population Stability (CV)") +
  xlab("Mass Independent rmax") +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks = element_line(linetype = "solid"),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y =element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        text = element_text(family = "Arial"))  +
  scale_color_manual(values = color_mapping) +
  ylim(-2.0, 0.5)

gg_partial_grwcv


##### MLM PC1 - Beccari et al 2024 #####
df_mammal_cv_fs_5 <- df_mammal_cv_fs %>% filter(cv_window == 5) #%>%
  #filter(! (Group == "Bats"| Group == "Marine Mammals"))

df_mammal_cv_fs_5 <- df_mammal_cv_fs_5 %>%
  mutate(mean_PC1 = scale(mean_PC1))
df_mammal_cv_fs_5$mean_PC1 <- as.numeric(df_mammal_cv_fs_5$mean_PC1)

lm_cv_bs_pc1 <- lm(mean_log_cv ~ log10(mean_body_mass), data = df_mammal_cv_fs_5)
summary(lm_cv_bs_pc1)
par(mfrow = c(2,2))
plot(lm_cv_bs_pc1)

lm_cv_pc1 <- lm(mean_log_cv ~ log10(mean_body_mass) + mean_PC1, data = df_mammal_cv_fs_5)
summary(lm_cv_pc1)
par(mfrow = c(2,2))
plot(lm_cv_pc1)

gg_partial_cvbspc1 <- visreg(lm_cv_pc1, "mean_body_mass", overlay=FALSE, partial=TRUE, gg=TRUE, points=list(col="NA")) + 
  geom_point(data = df_mammal_cv_fs_5, size = 1.5, aes(x = mean_body_mass, y = mean_log_cv, color = Group)) +
  geom_smooth(se = F, color = "black", method = "lm") +
  theme_classic() +
  ylab("Population Stability (CV)") +
  xlab("Body Size (g)") +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks = element_line(linetype = "solid"),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y =element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        text = element_text(family = "Arial")) +
  scale_color_manual(values = color_mapping) +
  scale_x_log10()

gg_partial_cvbspc1

gg_partial_cvpc1 <- visreg(lm_cv_pc1, "mean_PC1", overlay=FALSE, partial=TRUE, gg=TRUE, points=list(col="NA")) + 
  geom_point(data = df_mammal_cv_fs_5, size = 1.5, aes(x = mean_PC1, y = mean_log_cv, color = Group)) +
  geom_smooth(se = F, color = "black", method = "lm") +
  theme_classic() +
  ylab("Population Stability (CV)") +
  xlab("PC1 (Beccari et al 2024)") +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks = element_line(linetype = "solid"),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y =element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        text = element_text(family = "Arial")) +
  scale_color_manual(values = color_mapping)

gg_partial_cvpc1

##### MLM Metabolism #####

#scale residuals and ensure numeric vector for plotting
df_metcv <- df_metcv %>%
  mutate(residual_met = scale(residual_met))
df_metcv$residual_met <- as.numeric(df_metcv$residual_met)

mlm_metcv <- lm(mean_log_cv ~ log10(mean_body_mass) + residual_met, data = df_metcv)
summary(mlm_metcv)
par(mfrow = c(2,2))
plot(mlm_metcv)

gg_partial_metcvbs <- visreg(mlm_metcv, "mean_body_mass", overlay=FALSE, partial=TRUE, gg=TRUE, points=list(col="NA")) + 
  geom_point(data = df_metcv, size = 1.5, aes(x = mean_body_mass, y = mean_log_cv, color = Group.x)) +
  geom_smooth(se = F, color = "black", method = "lm") +
  theme_classic() +
  ylab("Population Stability (CV)") +
  xlab("Body Size (g)") +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks = element_line(linetype = "solid"),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y =element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        text = element_text(family = "Arial")) +
  scale_color_manual(values = color_mapping) +
  scale_x_log10(limits = c(x_min, 1e+07), breaks = xlog_breaks)  +
  ylim(-2.0, 0.5)

gg_partial_metcvbs

gg_partial_metcv <- visreg(mlm_metcv, "residual_met", overlay=FALSE, partial=TRUE, gg=TRUE, points=list(col="NA")) + 
  geom_point(data = df_metcv, size = 1.5, aes(x = residual_met, y = mean_log_cv, color = Group.x)) +
  geom_smooth(se = F, color = "black", method = "lm") +
  theme_classic() +
  ylab("Population Stability (CV)") +
  xlab("Mass Independent Specific Metabolism") +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks = element_line(linetype = "solid"),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y =element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        text = element_text(family = "Arial"))   +
  scale_color_manual(values = color_mapping) +
  ylim(-2.0, 0.5)

gg_partial_metcv


##### Partial Plots - Residual Life History Trait Body Size Panel #####
gg_residual_trait_cv <- ggarrange(gg_partial_grwcv, gg_partial_cvpc1,
                                  labels = c("B","C"),
                                  nrow = 1, ncol = 2)

gg_residual_trait_cv

#ggsave("../Mammalian-Fast-Slow-Stability/Figures/Figure 3 - Residual Traits + PC1 vs Log CV.jpeg", plot = gg_residual_trait_cv, width = 10, height = 8)

gg_big_panel <- ggarrange(gg_cv_bs_5, gg_residual_trait_cv,
                          ncol = 1, nrow = 2,
                          labels = c("A"))
gg_big_panel

#ggsave("../Mammalian-Fast-Slow-Stability/Figures/Figure 3 - CV body size + Residual Traits + PC1 vs Log CV.jpeg", plot = gg_big_panel, width =10, height = 10)


##### Phylogenetic Generalized Least Squares Regression #####
#Complete Phylogenetic Tree - From Phylacine (https://zenodo.org/records/3690867)
df_tree <- read.nexus("../Large Data Files/Phylacine v1.2 Data/Phylogenies/Complete_phylogeny.nex")

#Match Species names to tree tip labels
#CV ~ body size
tree_tips <- df_tree[[1]]$tip.label
df_cv_5_tree <- df_cv_5[match(tree_tips, df_cv_5$Binomial), ]
df_cv_5_tree_na <- na.omit(df_cv_5_tree)
df_cv_5_tree_na <- as.data.frame(df_cv_5_tree_na)
rownames(df_cv_5_tree_na) <- df_cv_5_tree_na$Binomial

phylolm_cv_model <- phylolm(formula = mean_log_cv ~ log10(mean_body_mass),
                            data = df_cv_5_tree_na, 
                            phy = df_tree[[1]], 
                            model = "lambda")

summary(phylolm_cv_model)

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

summary(phylolm_grw_model)

#PC1
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

summary(phylolm_fs_model)

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


summary(phylolm_metab_model)


##### Global Change Body Size Change Estimates - Implications for CV #####
#Calculate body mass changes (4%) with degrees Celsius warming Sensu Searing et al 2023...
options(scipen = 9)

df_cv_gc <- df_cv_5 %>%
  mutate(log_mean_body_mass = log10(mean_body_mass),
         Mass_g_2d = (mean_body_mass * 0.92),
         Mass_g_3d = (mean_body_mass  * 0.88),
         Mass_g_4d = (mean_body_mass  * 0.84),
         log_Mass_g_2d = log10(Mass_g_2d),
         log_Mass_g_3d = log10(Mass_g_3d),
         log_Mass_g_4d = log10(Mass_g_4d),
  Body_Mass_Group = case_when(
      mean_body_mass >= 2.0 & mean_body_mass <= 1000.0 ~ "Small",
      #mean_body_mass > 100.0 & mean_body_mass <= 1000.0 ~ "Small",
      mean_body_mass > 1000.0 & mean_body_mass <= 10000.0 ~ "Medium",
      #mean_body_mass > 10000.0 & mean_body_mass <= 100000.0 ~ "Large",
      mean_body_mass > 10000.0 ~ "Large",
      TRUE ~ "Other"
    ))


lm_gc_bscv <- lm(mean_log_cv ~ log_mean_body_mass, data = df_cv_gc)
summary(lm_gc_bscv)

#Predict mean_cv for projected changes in body mass
df_pred_cv_gc <- df_cv_gc %>%
  mutate(
    pred_logcv_baseline = predict(lm_gc_bscv),
    pred_logcv_2d = predict(lm_gc_bscv, newdata = data.frame(log_mean_body_mass = log_Mass_g_2d)),
    pred_logcv_3d = predict(lm_gc_bscv, newdata = data.frame(log_mean_body_mass = log_Mass_g_3d)),
    pred_logcv_4d = predict(lm_gc_bscv, newdata = data.frame(log_mean_body_mass = log_Mass_g_4d)),
    
    pred_cv_baseline = 10^(pred_logcv_baseline),
    pred_cv_2d_res = 10^(pred_logcv_2d),
    pred_cv_3d_res = 10^(pred_logcv_3d),
    pred_cv_4d_res = 10^(pred_logcv_4d))

#Create Body Size Groupings
df_delta_cv_gc <- df_pred_cv_gc %>%
  mutate(
    delta_cv2d = pred_cv_2d_res - pred_cv_baseline,
    delta_cv3d = pred_cv_3d_res - pred_cv_baseline,
    delta_cv4d = pred_cv_4d_res - pred_cv_baseline,
    
    percent_cv2d = (delta_cv2d/pred_cv_2d_res)*100,
    percent_cv3d = (delta_cv3d/pred_cv_3d_res)*100,
    percent_cv4d = (delta_cv4d/pred_cv_4d_res)*100)


df_percent_cv_gc <- df_delta_cv_gc %>%
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "Projection",
    values_to = "Predicted_mean_cv"
  ) %>%
  group_by(Projection) %>%
  reframe(Percent_Delta_CV = mean(Predicted_mean_cv))


#Set colors for Degree/Body Size Groupings
myPal <- c("#FEE391","#F7834D", "#DA464C", "#9E0142")


#GC Plot of Estiamted Percent Shift in CV with Projected Body Size Changes
gg_cv_gc <- ggplot(df_percent_cv_gc, aes(x = Projection, y = Percent_Delta_CV, fill = Projection)) +
  geom_col(width = 0.5) +
  ylim(0, 5) +
  theme_classic() +
  scale_fill_manual(values = myPal) +
  xlab("Projected Change in Body Size (per ˚C)") +
  ylab("Percent Change in Population Stability (CV)") +
  scale_x_discrete(labels=c("8% (2˚C)", "12% (3˚C)", "16% (4˚C)")) +
  theme(axis.line = element_line(linetype = "solid"),
      axis.ticks = element_line(linetype = "solid"),
      panel.background = element_rect(fill = NA),
      legend.key = element_rect(fill = NA),
      legend.background = element_rect(fill = NA),
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y =element_text(size = 12), 
      axis.title.x = element_text(size = 12), 
      text = element_text(family = "Arial"))

gg_cv_gc

#Estimate change in rmax with a shift in body size...
df_grw_gc <- df_grwcv %>%
  mutate(log_mean_body_mass = log10(Mass_g),
         Mass_g_2d = (Mass_g * 0.92),
         Mass_g_3d = (Mass_g  * 0.88),
         Mass_g_4d = (Mass_g  * 0.84),
         log_Mass_g_2d = log10(Mass_g_2d),
         log_Mass_g_3d = log10(Mass_g_3d),
         log_Mass_g_4d = log10(Mass_g_4d),
         Body_Mass_Group = case_when(
           Mass_g >= 2.0 & Mass_g <= 1000.0 ~ "Small",
           #mean_body_mass > 100.0 & mean_body_mass <= 1000.0 ~ "Small",
           Mass_g > 1000.0 & Mass_g <= 10000.0 ~ "Medium",
           #mean_body_mass > 10000.0 & mean_body_mass <= 100000.0 ~ "Large",
           Mass_g > 10000.0 ~ "Large",
           TRUE ~ "Other"
         ))

lm_grw_bs <- lm(log10(rmax) ~ log_mean_body_mass, data = df_grw_gc)
summary(lm_grw_bs)

#Predict mean_grw for projected changes in body mass
df_pred_grw_gc <- df_grw_gc %>%
  mutate(
    pred_loggrw_baseline = predict(lm_grw_bs),
    pred_loggrw_2d = predict(lm_grw_bs, newdata = data.frame(log_mean_body_mass = log_Mass_g_2d)),
    pred_loggrw_3d = predict(lm_grw_bs, newdata = data.frame(log_mean_body_mass = log_Mass_g_3d)),
    pred_loggrw_4d = predict(lm_grw_bs, newdata = data.frame(log_mean_body_mass = log_Mass_g_4d)),
    
    pred_grw_baseline = 10^(pred_loggrw_baseline),
    pred_grw_2d_res = 10^(pred_loggrw_2d),
    pred_grw_3d_res = 10^(pred_loggrw_3d),
    pred_grw_4d_res = 10^(pred_loggrw_4d))

#Create Body Size Groupings
df_delta_grw_gc <- df_pred_grw_gc %>%
  mutate(
    delta_grw2d = pred_grw_2d_res - pred_grw_baseline,
    delta_grw3d = pred_grw_3d_res - pred_grw_baseline,
    delta_grw4d = pred_grw_4d_res - pred_grw_baseline,
    
    percent_grw2d = (delta_grw2d/pred_grw_2d_res)*100,
    percent_grw3d = (delta_grw3d/pred_grw_3d_res)*100,
    percent_grw4d = (delta_grw4d/pred_grw_4d_res)*100)


df_percent_grw_gc <- df_delta_grw_gc %>%
  pivot_longer(
    cols = starts_with("percent"),
    names_to = "Projection",
    values_to = "Predicted_mean_grw"
  ) %>%
  group_by(Projection) %>%
  reframe(Percent_Delta_rmax = mean(Predicted_mean_grw))


#Set colors for Degree/Body Size Groupings
myPal <- c("#FEE391","#F7834D", "#DA464C", "#9E0142")


#GC Plot of Estiamted Percent Shift in CV with Projected Body Size Changes
gg_rmax_gc <- ggplot(df_percent_grw_gc, aes(x = Projection, y = Percent_Delta_rmax, fill = Projection)) +
  geom_col(width = 0.5) +
  theme_classic() +
  scale_fill_manual(values = myPal) +
  ylim(0, 5) +
  xlab("Projected Change in Body Size (per ˚C)") +
  ylab("Percent Change in rmax") +
  scale_x_discrete(labels=c("8% (2˚C)", "12% (3˚C)", "16% (4˚C)")) +
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks = element_line(linetype = "solid"),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y =element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        text = element_text(family = "Arial"))

gg_rmax_gc


gg_globchange <- ggarrange(gg_cv_gc, gg_rmax_gc,
                             nrow = 2, ncol = 1)
gg_globchange

#ggsave("../Mammalian-Fast-Slow-Stability/Figures/Figure 4 - Percent Change CV Rmax - Global Change.jpeg", plot = gg_globchange, width = 6, height = 8)


df_gc_theory <- read.csv("../Mammalian-Fast-Slow-Stability/Data/Theory/Theory GC Percent.csv")

gg_theory_gc <- ggplot(df_gc_theory, aes(x = Percent_Increase, y = Percent_Change_CV_High)) +
  theme_classic() +
  #scale_fill_manual(values = myPal) +
  xlab("Projected Change in rmax") +
  ylab("Theoretical Percent Change in CV") +
  ylim(0, 45) +
  geom_rect(aes(xmin = 0, xmax = 0.02047822, ymin = -Inf, ymax = Inf), fill = myPal[1], alpha = 0.5) +
  geom_rect(aes(xmin = 0.02047822, xmax = 0.03122348, ymin = -Inf, ymax = Inf), fill = myPal[2], alpha = 0.5) +
  geom_rect(aes(xmin = 0.03122348, xmax = 0.05, ymin = -Inf, ymax = Inf), fill = myPal[3], alpha = 0.5) +
  geom_smooth(se = F, linewidth = 3, color = "black") +
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks = element_line(linetype = "solid"),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y =element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        text = element_text(family = "Arial"))

gg_theory_gc

#ggsave("../Mammalian-Fast-Slow-Stability/Figures/Figure 4 - Percent Change CV Rmax - Global Change Theory.jpeg", plot = gg_theory_gc, width = 8, height = 6)


gg_globchange_all <- ggarrange(gg_globchange, gg_theory_gc,
                           nrow = 1, ncol = 2)
gg_globchange_all

#ggsave("../Mammalian-Fast-Slow-Stability/Figures/Figure 4 - All.jpeg", plot = gg_globchange_all, width = 10, height = 6)

