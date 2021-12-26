library("tidyverse")
library("RColorBrewer")
library("ggrepel")
library("Momocs")
library("patchwork")
library("grid")
library("mgcv")

rm(list=ls())

# Import helper file. This file contain proprietary code and
# will not be supplied with the rest of the code.
# All functions used from this file will be highlighted in the script.
source("../imp_func.R", local = TRUE)

# Import model data
ci_mody_norm <- read_csv("../data/predframe_eff_size.csv")

# Order variables for ggplot
ci_mody_norm$vars<- factor(ci_mody_norm$vars, levels=c("T", "Alk", "SIR", "SH"))

# Import shell shape file
nucella_coe <- readRDS("../data/shell_coe_object.rds")

# Calculate PCA from shapes
nucella_pca <- PCA(nucella_coe)

# Plot PC contribution 
PC_contrib <- PCcontrib(nucella_pca, nax=1:3, sd.r = c(-3,0,3))

# Extract Shell shapes from PCA
PC1_1 <- as.matrix(PC_contrib$shp$`1.1`[1:2])
PC1_3 <- as.matrix(PC_contrib$shp$`3.1`[1:2])
PC2_1 <- as.matrix(PC_contrib$shp$`1.2`[1:2])
PC2_3 <- as.matrix(PC_contrib$shp$`3.2`[1:2])
PC3_1 <- as.matrix(PC_contrib$shp$`1.3`[1:2])
PC3_3 <- as.matrix(PC_contrib$shp$`3.3`[1:2])

isopalette <- colorRampPalette(rev(brewer.pal(10, "RdBu")))
maxp <- rev(brewer.pal(11, "RdBu"))

PC1 <- wrap_elements(full = ~ tps_iso(PC1_1, PC1_3, shp.lwd = c(1.5, 1.5), shp.border = c(maxp[2], maxp[9]),
               amp = 0, iso.nb = 10000, iso.levels = 10, legend = FALSE, palette = isopalette,
               grid = F))
PC2 <- wrap_elements(full = ~ tps_iso(PC2_1, PC2_3, shp.lwd = c(1.5, 1.5), shp.border = c(maxp[2], maxp[9]),
               amp = 0, iso.nb = 10000, iso.levels = 10, legend = FALSE, palette = isopalette,
               grid = F))
PC3 <- wrap_elements(full = ~ tps_iso(PC3_1, PC3_3, shp.lwd = c(1.5, 1.5), shp.border = c(maxp[2], maxp[9]),
               amp = 0, iso.nb = 10000, iso.levels = 10, legend = FALSE, palette = isopalette,
               grid = F))


# Plot panels -------------------------------------------------------------
eff_pc1 <- ci_mody_norm %>% 
  filter(PC == "shape-PC1") %>% 
  ggplot(aes(vars, ymin = lower, ymax = upper, geom = "pointrange")) +
  geom_hline(yintercept = 0, alpha = I(5/12), lty = 2) + 
  geom_errorbar(position = position_dodge(width = 0.3), width = 0.2, colour="grey80") + 
  geom_point(aes(vars, est., fill=est.), position = position_dodge(width = 0.3), size = 1.5, shape = 21) + 
  scale_fill_gradient(low = maxp[2], high = maxp[9], limits = c(-0.4,0.4)) +
  theme_science() +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        text = element_text(family ="Times"),
        title =element_text(family ="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none") +
  ylim(-0.4,0.4)

eff_pc2 <- ci_mody_norm %>% 
  filter(PC == "shape-PC2") %>% 
  ggplot(aes(vars, ymin = lower, ymax = upper, geom = "pointrange")) +
  geom_hline(yintercept = 0, alpha = I(5/12), lty = 2) + 
  geom_errorbar(position = position_dodge(width = 0.3), width = 0.2, colour="grey80") + 
  geom_point(aes(vars, est., fill=est.), position = position_dodge(width = 0.3), size = 1.5, shape = 21) + 
  scale_fill_gradient(low = maxp[2], high = maxp[9], limits = c(-0.4,0.4)) +
  theme_science() +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        text = element_text(family ="Times"),
        title =element_text(family ="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")  +
  ylim(-0.4,0.4)

eff_pc3 <- ci_mody_norm %>% 
  filter(PC == "shape-PC3") %>% 
  ggplot(aes(vars, ymin = lower, ymax = upper, geom = "pointrange")) +
  geom_hline(yintercept = 0, alpha = I(5/12), lty = 2) + 
  geom_errorbar(position = position_dodge(width = 0.3), width = 0.2, colour="grey80") + 
  geom_point(aes(vars, est., fill=est.), position = position_dodge(width = 0.3), size = 1.5, shape = 21) +   
  scale_fill_gradient(low = maxp[2], high = maxp[9], limits = c(-0.4,0.4)) +
  theme_science() +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        text = element_text(family ="Times"),
        title =element_text(family ="Times"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  ylim(-0.4,0.4)

# Combine plots -----------------------------------------------------------

pdf("eff_size.pdf", width = 3, height = 4, onefile = F)

pc_panel_l <- (PC1 / PC2 / PC3) + plot_layout(heights = c(1,1,0.9))

eff_panel_r <- (eff_pc1 / eff_pc2 / eff_pc3) + plot_layout(heights = c(1,1,1))

(pc_panel_l | eff_panel_r) + plot_layout(widths = c(1,0.9))

grid.text("Standardised regression estimates", rot=90,
          x = 0.45, y = 0.5, gp=gpar(fontsize=10, fontfamily = "Times"))

grid.text("A", x = 0.1, y = 0.95, gp=gpar(fontsize=10, fontfamily = "Times"))
grid.text("B", x = 0.1, y = 0.65, gp=gpar(fontsize=10, fontfamily = "Times"))
grid.text("C", x = 0.1, y = 0.35, gp=gpar(fontsize=10, fontfamily = "Times"))

grid.text("*", x = 0.805, y = 0.81, gp=gpar(fontsize=10, fontfamily = "Times"))
grid.text("*", x = 0.888, y = 0.88, gp=gpar(fontsize=10, fontfamily = "Times"))

grid.text("*", x = 0.725, y = 0.515, gp=gpar(fontsize=10, fontfamily = "Times"))
grid.text("*", x = 0.805, y = 0.595, gp=gpar(fontsize=10, fontfamily = "Times"))
grid.text("*", x = 0.888, y = 0.57, gp=gpar(fontsize=10, fontfamily = "Times"))

grid.text("*", x = 0.805, y = 0.262, gp=gpar(fontsize=10, fontfamily = "Times"))
grid.text("*", x = 0.888, y = 0.28, gp=gpar(fontsize=10, fontfamily = "Times"))

dev.off()
