library("tidyverse")
library("RColorBrewer")
library("Momocs")
library("patchwork")

# Import helper file. This file contain proprietary code and
# will not be supplied with the rest of the code.
# All functions used from this file will be highlighted in the script.
source("../imp_func.R", local = TRUE)

# Import shell shape file
nucella_PCs <- read_csv("../data/shell_outline_PCA.csv")
nucella_coe <- readRDS("../data/shell_coe_object.rds")

# Import env data
df_env <- read_csv("../data/environmental_comp.csv")

# Merge Temp with Pcs data frame
nucella_PCs <- df_env %>% 
  select(Temp, location) %>% 
  full_join(nucella_PCs)

# Calculate PCA from shapes
nucella_pca <- PCA(nucella_coe)

# Plot PC contribution 
PC_contrib <- PCcontrib(nucella_pca, nax=1:3, sd.r = c(-3,0,3))

# Scale PC1 and PC2 
pcdf <- nucella_PCs %>%
  mutate(
    PC1 = range01(PC1),
    PC2 = range01(PC2)
  ) %>% 
  mutate(zone = case_when(
    Temp < 9 ~ "cold",
    Temp >= 9 & Temp <= 12 ~ "intermediate",
    Temp > 12 ~ "warm"
  ))

# Extract Shell shapes from PCA
PC1_1 <- PC_contrib$shp$`1.1`[1:2]
PC1_3 <- PC_contrib$shp$`3.1`[1:2]
PC2_1 <- PC_contrib$shp$`1.2`[1:2]
PC2_3 <- PC_contrib$shp$`3.2`[1:2]

# Plot subs
sub1 <- ggplot(PC1_1 , aes(x,y)) + 
  geom_polygon(fill="lightgrey", size=0.35, colour="grey20") + 
  theme_void() + coord_fixed()
sub2 <- ggplot(PC1_3 , aes(x,y)) + 
  geom_polygon(fill="lightgrey", size=0.35, colour="grey20") + 
  theme_void() + coord_fixed()
sub3 <- ggplot(PC2_1 , aes(x,y)) + 
  geom_polygon(fill="lightgrey", size=0.35, colour="grey20") +
  theme_void() + coord_fixed()
sub4 <- ggplot(PC2_3 , aes(x,y)) + 
  geom_polygon(fill="lightgrey", size=0.35, colour="grey20") + 
  theme_void()+ coord_fixed()

pdf("PCA_plot.pdf",
    onefile = F, height = 3)

# Plot shape-PCs
ggplot(pcdf, aes(x=PC2, y=PC1, group=location)) +
  geom_hline(yintercept = 0.5, colour="grey", size=0.2, alpha=0.5) +
  geom_vline(xintercept = 0.5, colour="grey", size=0.2, alpha=0.5) +
  annotation_custom(grob=ggplotGrob(sub3), 
                    ymin = 0.35, ymax = 0.65, xmax=0) +
  annotation_custom(grob=ggplotGrob(sub4), 
                    ymin = 0.35, ymax = 0.65, xmin=1.05) + 
  annotation_custom(grob=ggplotGrob(sub1), 
                    ymin = -0.25, ymax=0.05, xmin=-0.2) +
  annotation_custom(grob=ggplotGrob(sub2), 
                    ymin = 0.95, ymax=1.25, xmin=-0.2) + 
  geom_point(pch=21, fill="white", colour="grey", alpha = 0.7) +
  stat_ellipse(type = "norm", linetype = 2, size = 0.4, level = 0.65, alpha=0.5) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    title = element_text(family="Times"),
    text = element_text(size = 10,family="Times")
  ) + 
  ylab(paste("shape-PC1"," - ",round(nucella_pca$eig[1]*100,2),"%",sep="")) +
  xlab(paste("shape-PC2"," - ",round(nucella_pca$eig[2]*100,2),"%",sep="")) +
  xlim(c(-0.2,1.2)) +
  ylim(c(-0.2,1.2)) +
  coord_fixed() 

dev.off()
