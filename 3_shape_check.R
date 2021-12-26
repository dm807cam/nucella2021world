library("tidyverse")
library("cmocean")
library("Momocs")
library("ggrepel")
library("patchwork")

rm(list=ls())

# Raster package masks dplyr select
select <- dplyr::select

# Import helper file. This file contain proprietary code and
# will not be supplied with the rest of the code.
# All functions used from this file will be highlighted in the script.
source("imp_func.R", local = TRUE)

# Import data -------------------------------------------------------------
# Read in files
df_caliper <- read_csv("data/shell_dimension_caliper.csv")
df_pca <- read_csv("data/shell_outline_PCA.csv")
df_env <- read_csv("data/environmental_comp.csv")
df_exp <- read_csv("data/exposure_levels.csv")

# Check files
glimpse(df_caliper)
glimpse(df_pca)
glimpse(df_env)
glimpse(df_exp)

# Restructure data for analysis -------------------------------------------
# Caliper measurements
df_caliper_long <- df_caliper %>%
# Create unique row ID  
  tibble::rowid_to_column("ID") %>%
# Pivot from wide to long format  
  pivot_longer(
    cols = c(
      -ID,-location,-lat,-lon,-year,-label
    ),
    names_to = 'assigned',
    values_to = 'value'
  )  %>%
# Extract strings from labels  
  mutate(
    shape = str_extract(assigned, "([^\\s]+\\s+[^\\s]+)"),
    repitition = str_extract(assigned, "[0-9]")
  ) %>%
# Remove redundant variable "assigned"
  dplyr::select(-assigned,-label,-year) %>%
# Group unique data to summarise duplicate measurements
# Note: Every shell descriptor on every shell was measured four times!  
    group_by(ID,
           location,
           lat,
           lon,
           shape) %>%
# Calculate average and standard deviation of all descriptors  
  summarise(
    sd = sd(value),
    value = mean(value)
  ) %>%
# Don't forget to un-group!  
  ungroup() %>%
# Make sure variables are "factors"  
  mutate(
    flocation = as.factor(location),
    fshape = as.factor(shape)
  ) %>% 
  select(-location, -shape)

# Turn to wide again to calculate meaningful shape ratios
df_caliper <- df_caliper_long %>% 
  select(-sd) %>% 
  pivot_wider(id_cols = c(flocation,lat,lon,ID), 
              names_from = fshape,
              values_from = value) %>% 
  mutate(shell_height = `shell height`,
         shell_width = `shell width`,
         aperture_size = pi * (0.5 * `aperture width`) * (0.5 * `aperture height`)) %>% 
  select(-c(`aperture height`:`shell width`))
  
# Pivot from wide to long format  
df_caliper_long <- df_caliper %>% 
  pivot_longer(
    cols = c(
      -ID,-flocation,-lat,-lon, -shell_height
    ),
    names_to = 'fshape',
    values_to = 'value') %>% 
  select (-ID)

# PCA data
df_pca_long <- df_pca %>% 
# Create unique row ID  
  tibble::rowid_to_column("ID") %>%
# Pivot from wide to long format  
  pivot_longer(
    cols = c(-ID, -sampleID, -sampleNo, -shell_length, -location, -lat, -lon),
    names_to = 'shape',
    values_to = 'value'
  ) %>%
# Make sure variables are in the correct format  
  mutate(
    flocation = as.factor(location),
    fshape = as.factor(shape)
  ) %>% 
  select(-location, -shape,-ID)

# Merge dataframes ------------------------
# Merge with environmental data frame
df_pca_long <- df_env %>%
  mutate(flocation = as.factor(location)) %>% 
  # un-select duplicate or unwanted variables
  select(-flag, -dist, -lat, -lon) %>% 
  full_join(df_pca_long, by="flocation")

# Merge with environmental data frame
df_caliper_long <- df_env %>%
  mutate(flocation = as.factor(location)) %>% 
  # un-select duplicate or unwanted variables
  select(-flag, -dist, -lat, -lon) %>% 
  full_join(df_caliper_long, by="flocation")

# Export data
write_csv(df_pca_long, "data/df_pca_long.csv")
write_csv(df_caliper_long, "data/df_caliper_long.csv")

# Data exploration --------------------------------------------------------
# Follow data exploration routine as in Zuur 2010

# 1. Outliers Y & X
# 2. Homogeneity Y
# 3. Normality Y
# 4. Zero inflation Y
# 5. Collinearity X
# 6. Relationships + Interactions Y & X
# 7. Independence Y

# 1. Outliers Y & X (cleveland dotplot and boxplots) -----------------------------------
dp_pca <- df_pca %>% 
  select(PC1:PC5) %>% 
  as.matrix()

lattice::dotplot(dp_pca, group =F, 
                 scales = list(x = list(relation ="free"), 
                               y = list(relation = "free"), draw = F))

dp_cal <- df_caliper %>% 
  select(5:7) %>% 
  as.matrix()

lattice::dotplot(dp_cal, group =F, 
                       scales = list(x = list(relation ="free"), 
                                     y = list(relation = "free"), draw = F))
# Larger variability between caliper groups than PCA groups

# Boxplots
n1 <- ggplot(df_pca_long, aes(flocation, value)) + 
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~fshape, scales="free_y")

n2 <- ggplot(df_caliper_long, aes(flocation, value)) + 
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~fshape, scales="free_y")

n1 + n2

# 2. Homogeneity ----------------------------------------------------------
n1 <- ggplot(df_pca_long, aes(flocation, value)) + 
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~fshape)

n2 <- ggplot(df_caliper_long, aes(flocation, value)) + 
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~fshape)

n1 + n2
# Homogeneity among observations per group (fshape) can be assumed. 
# BUT variance between groups vastly different.

# 3. Normality ------------------------------------------------------------
n1 <- df_pca_long %>% 
    ggplot(aes(sample = value)) +
    stat_qq(alpha=0.3, colour="black") +
    stat_qq_line() +
    theme_bw() +
    theme(text = element_text(family="Times"),
          title = element_text(family="Times"),
          legend.position = "none",
          aspect.ratio = 1) +
  xlab("Theoretical Quantiles") +
  facet_wrap(~fshape, scales = "free_y")

n2 <- df_caliper_long %>% 
  ggplot(aes(sample = value)) +
  stat_qq(alpha=0.3, colour="black") +
  stat_qq_line() +
  theme_bw() +
  theme(text = element_text(family="Times"),
        title = element_text(family="Times"),
        legend.position = "none",
        aspect.ratio = 1) +
  xlab("Theoretical Quantiles") +
  facet_wrap(~fshape, scales = "free_y")

n1 + n2 

# 4. Zero inflation -------------------------------------------------------
# NAs in the data frame could cause problems with the models, check if NAs are present
colSums(is.na(df_pca_long))
colSums(is.na(df_caliper_long))
#No NA's

# 5. Multicollinearity ----------------------------------------------------
# Check for (multi-) colinearity among descriptive variables
# and step wise exclusion of less relevant variables 
# based on a priori knowledge.

df_pca_long %>% 
  ungroup() %>% 
  select(Temp,
         Psal,
         #pH,
         #lat,
         SIR,
         #DIC,
         #HCO3,
         #ALK,
         #OmegaCalcite,
         #OmegaAragonite
         ) %>% 
  corvif() %>% 
  round(2)

df_caliper_long %>% 
  ungroup() %>% 
  select(Temp,
         Psal,
         #pH,
         #lat,
         SIR,
         #DIC,
         #HCO3,
         ALK,
         #OmegaCalcite,
         #OmegaAragonite
  ) %>% 
  corvif() %>% 
  round(2)

# 6. Relationships + Interactions Y & X -----------------------------------
# 6.1 Temperature ---------------------------------------------------------
n1 <- df_pca_long %>%
  ggplot(aes(Temp, value)) +
  geom_point(fill="grey60",pch=21, alpha=0.5) +
  geom_smooth(method ="loess", span = 1, size= 0.5, colour="orange") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  facet_wrap(~fshape,scales = "free_y",ncol = 2)

n2 <- df_caliper_long %>% 
  ggplot(aes(Temp, value)) +
  geom_point(fill="grey60",pch=21, alpha=0.5) +
  geom_smooth(method ="loess", span = 1, size= 0.5, colour="orange") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  facet_wrap(~fshape,scales = "free_y",ncol = 2)

n1 + n2

# 6.2 Salinity ------------------------------------------------------------
n1 <- df_pca_long %>% 
  ggplot(aes(Psal, value)) +
  geom_point(fill="grey60",pch=21, alpha=0.5) +
  geom_smooth(method ="loess", span = 1, size= 0.5, colour="orange") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  facet_wrap(~fshape,scales = "free_y",ncol = 2)

n2 <- df_caliper_long %>% 
  ggplot(aes(Psal, value)) +
  geom_point(fill="grey60",pch=21, alpha=0.5) +
  geom_smooth(method ="loess", span = 1, size= 0.5, colour="orange") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  facet_wrap(~fshape,scales = "free_y",ncol = 2)

n1 + n2

# 6.3 SIR ----------------------------------------------------------------
n1 <- df_pca_long %>% 
  ggplot(aes(SIR, value)) +
  geom_point(fill="grey60",pch=21, alpha=0.5) +
  geom_smooth(method ="loess", span = 1, size= 0.5, colour="orange") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  facet_wrap(~fshape,scales = "free_y",ncol = 2)

n2 <- df_caliper_long %>% 
  ggplot(aes(SIR, value)) +
  geom_point(fill="grey60",pch=21, alpha=0.5) +
  geom_smooth(method ="loess", span = 1, size= 0.5, colour="orange") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  facet_wrap(~fshape,scales = "free_y",ncol = 2)

n1 + n2

# 6.4 ALK ----------------------------------------------------------------
n1 <- df_pca_long %>% 
  ggplot(aes(ALK, value)) +
  geom_point(fill="grey60",pch=21, alpha=0.5) +
  geom_smooth(method ="loess", span = 1, size= 0.5, colour="orange") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  facet_wrap(~fshape,scales = "free_y",ncol = 2)

n2 <- df_caliper_long %>% 
  ggplot(aes(ALK, value)) +
  geom_point(fill="grey60",pch=21, alpha=0.5) +
  geom_smooth(method ="loess", span = 1, size= 0.5, colour="orange") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  facet_wrap(~fshape,scales = "free_y",ncol = 2)

n1 + n2

# 6.5 pH ------------------------------------------------------------------
n1 <- df_pca_long %>% 
  ggplot(aes(pH, value)) +
  geom_point(fill="grey60",pch=21, alpha=0.5) +
  geom_smooth(method ="loess", span = 1, size= 0.5, colour="orange") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  facet_wrap(~fshape,scales = "free_y",ncol = 2)

n2 <- df_caliper_long %>% 
  ggplot(aes(pH, value)) +
  geom_point(fill="grey60",pch=21, alpha=0.5) +
  geom_smooth(method ="loess", span = 1, size= 0.5, colour="orange") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  facet_wrap(~fshape,scales = "free_y",ncol = 2)

n1 + n2

# 6.6 Shell height --------------------------------------------------------
df_pca_long %>% 
  ggplot(aes(shell_length, value)) +
  geom_point(fill="grey60",pch=21, alpha=0.5) +
  geom_smooth(method ="loess", span = 1, size= 0.5, colour="orange") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  facet_wrap(~fshape,scales = "free_y",ncol = 2)

# 6.7 Latitude ------------------------------------------------------------
n1 <- df_pca_long %>% 
  ggplot(aes(lat, value)) +
  geom_point(fill="grey60",pch=21, alpha=0.5) +
  geom_smooth(method ="loess", span = 1, size= 0.5, colour="orange") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  facet_wrap(~fshape,scales = "free_y",ncol = 2)

n2 <- df_caliper_long %>% 
  ggplot(aes(lat, value)) +
  geom_point(fill="grey60",pch=21, alpha=0.5) +
  geom_smooth(method ="loess", span = 1, size= 0.5, colour="orange") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  facet_wrap(~fshape,scales = "free_y",ncol = 2)

n1 + n2

# Check for relationship of shape-PCs to latitude -------------------------

# Linear regression
par(mfrow = c(2, 2))

# shape-PC1
lm_pc1 <- lm(value ~ lat, data = filter(df_pca_long, fshape=="PC1"))
anova(lm_pc1)
plot(lm_pc1)

# shape-PC2
lm_pc2 <- lm(value ~ lat, data = filter(df_pca_long, fshape=="PC2"))
anova(lm_pc2)
plot(lm_pc2)

# shape-PC3
lm_pc3 <- lm(value ~ lat, data = filter(df_pca_long, fshape=="PC3"))
anova(lm_pc3)
plot(lm_pc3)

# shape-PC4
lm_pc4 <- lm(value ~ lat, data = filter(df_pca_long, fshape=="PC4"))
anova(lm_pc4)
plot(lm_pc4)

# shape-PC5
lm_pc5 <- lm(value ~ lat, data = filter(df_pca_long, fshape=="PC5"))
anova(lm_pc5)
plot(lm_pc5)

# absolute aperture size
# Log-transform descriptor
df_caliper_as <- filter(df_caliper_long, fshape=="aperture_size") %>% 
  mutate(log_value = log10(value))

df_caliper_as %>% 
  ggplot(aes(sample = log_value)) +
  stat_qq(alpha=0.3, colour="black") +
  stat_qq_line() +
  theme_bw() +
  theme(text = element_text(family="Times"),
        title = element_text(family="Times"),
        legend.position = "none",
        aspect.ratio = 1) +
  xlab("Theoretical Quantiles") 

lm_as <- lm(value ~ lat, data = df_caliper_as)
anova(lm_as)
plot(lm_as)
