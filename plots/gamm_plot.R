library("tidyverse")
library("RColorBrewer")
library("latex2exp")
library("ggrepel")
library("Momocs")
library("patchwork")
library("grid")
library("mgcv")

# Import helper file. This file contain proprietary code and
# will not be supplied with the rest of the code.
# All functions used from this file will be highlighted in the script.
source("../imp_func.R", local = TRUE)

# Import model data
predframe_temp <- read_csv("../data/predframe_temp.csv")
predframe_SIR <- read_csv("../data/predframe_SIR.csv")
predframe_alk <- read_csv("../data/predframe_alk.csv")
predframe_ash <- read_csv("../data/predframe_ash.csv")
mody <- readRDS("../data/shell_shape_model.rds")

summary(mody)

# Import shell shape file
nucella_coe <- readRDS("../data/shell_coe_object.rds")

# Calculate PCA from shapes
nucella_pca <- PCA(nucella_coe)

# Plot PC contribution 
PC_contrib <- PCcontrib(nucella_pca, nax=1:3, sd.r = c(-2,0,2))

# Extract Shell shapes from PCA
PC1_1 <- as.matrix(PC_contrib$shp$`1.1`[1:2])
PC1_3 <- as.matrix(PC_contrib$shp$`3.1`[1:2])
PC2_1 <- as.matrix(PC_contrib$shp$`1.2`[1:2])
PC2_3 <- as.matrix(PC_contrib$shp$`3.2`[1:2])
PC3_1 <- as.matrix(PC_contrib$shp$`1.3`[1:2])
PC3_3 <- as.matrix(PC_contrib$shp$`3.3`[1:2])

isopalette <- colorRampPalette(rev(brewer.pal(10, "RdBu")))
maxp <- rev(brewer.pal(11, "RdBu"))

PC1 <- wrap_elements(full = ~ tps_iso(PC1_1, PC1_3, shp.lwd = c(2, 2), shp.border = c(maxp[2], maxp[9]),
                                      amp = 0.3, iso.nb = 10000, iso.levels = 10, legend = FALSE, palette = isopalette,
                                      grid = F))
PC2 <- wrap_elements(full = ~ tps_iso(PC2_1, PC2_3, shp.lwd = c(2, 2), shp.border = c(maxp[2], maxp[9]),
                                      amp = 0.3, iso.nb = 10000, iso.levels = 10, legend = FALSE, palette = isopalette,
                                      grid = F))
PC3 <- wrap_elements(full = ~ tps_iso(PC3_1, PC3_3, shp.lwd = c(2, 2), shp.border = c(maxp[2], maxp[9]),
                                      amp = 0.3, iso.nb = 10000, iso.levels = 10, legend = FALSE, palette = isopalette,
                                      grid = F))

# Plot panels -------------------------------------------------------------
scaleFUN <- function(x) sprintf("%.1f", x)
tpc1 <- predframe_temp %>% 
  filter(fshape == "PC1") %>% 
  ggplot() + geom_smooth(mapping=aes(Temp, preds, color=..y..), lwd = 0.5, se = F) +
  geom_smooth(aes(Temp, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(Temp, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  theme_science() +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_colour_gradient(low = maxp[2], high = maxp[9]) +
  scale_y_continuous(labels = scaleFUN, expand = c(0.1,0.1)) +
  xlab("Temperature (°C)") + ylab("shape-PC1") +
  annotate("text", x=Inf, y = Inf, label = "***", vjust=2, hjust=1.5, family="Times")

tpc2 <- predframe_temp %>% 
  filter(fshape == "PC2") %>% 
  ggplot() + geom_smooth(aes(Temp, preds), col="grey60", lwd = 0.5, se = F) +
  geom_smooth(aes(Temp, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(Temp, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  theme_science() +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_colour_gradient(low = maxp[2], high = maxp[9]) +
  scale_y_continuous(labels = scaleFUN, expand = c(0.1,0.1)) +
  xlab("Temperature (°C)") + ylab("shape-PC2") +
  annotate("text", x=Inf, y = Inf, label = "n.s.", vjust=2, hjust=1.5, family="Times")

tpc3 <- predframe_temp %>% 
  filter(fshape == "PC3") %>% 
  ggplot() + geom_smooth(aes(Temp, preds), col ="grey60", lwd = 0.5, se = F) +
  geom_smooth(aes(Temp, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(Temp, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  theme_science() +
  theme(legend.position = "none") +
  scale_colour_gradient(low = maxp[2], high = maxp[9]) +
  scale_y_continuous(labels = scaleFUN, expand = c(0.1,0.1)) +
  xlab("Temperature (°C)") + ylab("shape-PC3") +
  annotate("text", x=Inf, y = Inf, label = "n.s.", vjust=2, hjust=1.5, family="Times")

# Alk 
apc1 <- predframe_alk %>% 
  filter(fshape == "PC1") %>% 
  mutate(ALK = round(ALK * 10 ^ 6,0)) %>% 
  ggplot() + geom_smooth(aes(ALK, preds, color=..y..), lwd = 0.5, se = F) +
  geom_smooth(aes(ALK, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(ALK, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  theme_science() +
  theme(axis.title = element_blank(),
        legend.position = "none") +
  scale_colour_gradient(low = maxp[2], high = maxp[9]) +
  scale_y_continuous(labels = scaleFUN, expand = c(0.1,0.1)) +
  xlab(expression(paste("Alkalinity (µmol ", kg^{-1},")"))) + ylab("shape-PC1") +
  annotate("text", x=Inf, y = Inf, label = "*", vjust=2, hjust=1.5, family="Times")

apc2 <- predframe_alk %>% 
  filter(fshape == "PC2") %>% 
  mutate(ALK = round(ALK * 10 ^ 6,0)) %>% 
  ggplot() + geom_smooth(aes(ALK, preds, color=..y..), lwd = 0.5, se = F) +
  geom_smooth(aes(ALK, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(ALK, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  theme_science() +
  theme(axis.title = element_blank(),
        legend.position = "none") +
  scale_colour_gradient(low = maxp[2], high = maxp[9]) +
  scale_y_continuous(labels = scaleFUN, expand = c(0.1,0.1)) +
  xlab(expression(paste("Alkalinity (µmol ", kg^{-1},")"))) + ylab("shape-PC2") +
  annotate("text", x=Inf, y = Inf, label = "***", vjust=2, hjust=1.5, family="Times")

apc3 <- predframe_alk %>% 
  filter(fshape == "PC3") %>% 
  mutate(ALK = round(ALK * 10 ^ 6,0)) %>% 
  ggplot() + geom_smooth(aes(ALK, preds, color=..y..), lwd = 0.5, se = F) +
  geom_smooth(aes(ALK, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(ALK, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  theme_science() +
  theme(axis.title.y = element_blank(),
        legend.position = "none") +
  scale_colour_gradient(low = maxp[2], high = maxp[9]) +
  scale_y_continuous(labels = scaleFUN, expand = c(0.1,0.1)) +
  xlab(expression(paste("Alkalinity (µmol ", kg^{-1},")"))) + ylab("shape-PC3") +
  annotate("text", x=Inf, y = Inf, label = "***", vjust=2, hjust=1.5, family="Times")

# SIR
spc1 <- predframe_SIR %>% 
  filter(fshape == "PC1") %>% 
  ggplot() + geom_smooth(aes(SIR, preds, color=..y..), lwd = 0.5, se = F) +
  geom_smooth(aes(SIR, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(SIR, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  theme_science() +
  theme(axis.title = element_blank(),
        legend.position = "none") +
  scale_colour_gradient(low = maxp[2], high = maxp[9]) +
  scale_y_continuous(labels = scaleFUN, expand = c(0.1,0.1)) +
  xlab("SIR") + ylab("shape-PC1") +
  annotate("text", x=Inf, y = Inf, label = "***", vjust=2, hjust=2, family="Times")

spc2 <- predframe_SIR %>% 
  filter(fshape == "PC2") %>% 
  ggplot() + geom_smooth(aes(SIR, preds, color=..y..), lwd = 0.5, se = F) +
  geom_smooth(aes(SIR, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(SIR, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  theme_science() +
  theme(axis.title = element_blank(),
        legend.position = "none") +
  scale_colour_gradient(low = maxp[2], high = maxp[9]) +
  scale_y_continuous(labels = scaleFUN, expand = c(0.1,0.1)) +
  xlab("SIR") + ylab("shape-PC2") +
  annotate("text", x=Inf, y = Inf, label = "***", vjust=2, hjust=2, family="Times")

spc3 <- predframe_SIR %>% 
  filter(fshape == "PC3") %>% 
  ggplot() + geom_smooth(aes(SIR, preds, color=..y..), lwd = 0.5, se = F) +
  geom_smooth(aes(SIR, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(SIR, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  theme_science() +
  theme(axis.title.y = element_blank(),
        legend.position = "none") +
  scale_colour_gradient(low = maxp[2], high = maxp[9]) +
  scale_y_continuous(labels = scaleFUN, expand = c(0.1,0.1)) +
  xlab("SIR") + ylab("shape-PC3") +
  annotate("text", x=Inf, y = Inf, label = "*", vjust=2, hjust=2, family="Times")

# Shell height
slpc1 <- predframe_ash %>% 
  filter(fshape == "PC1") %>% 
  mutate(nSH = shell_length/mean(shell_length)) %>% 
  ggplot() + geom_smooth(aes(nSH, preds, color=..y..), lwd = 0.5, se = F) +
  geom_smooth(aes(nSH, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(nSH, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  theme_science() +
  theme(axis.title = element_blank(),
        legend.position = "none") +
  scale_colour_gradient(low = maxp[2], high = maxp[9]) +
  scale_y_continuous(labels = scaleFUN, expand = c(0.1,0.1)) +
  xlab("nSH") + ylab("shape-PC1") +
  annotate("text", x=Inf, y = Inf, label = "***", vjust=2, hjust=2, family="Times")

slpc2 <- predframe_ash %>% 
  filter(fshape == "PC2") %>% 
  mutate(nSH = shell_length/mean(shell_length)) %>% 
  ggplot() + geom_smooth(aes(nSH, preds, color=..y..), lwd = 0.5, se = F) +
  geom_smooth(aes(nSH, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(nSH, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  theme_science() +
  theme(axis.title = element_blank(),
        legend.position = "none") +
  scale_colour_gradient(low = maxp[2], high = maxp[9]) +
  scale_y_continuous(labels = scaleFUN, expand = c(0.1,0.1)) +
  xlab("nSH") + ylab("shape-PC2") +
  annotate("text", x=Inf, y = Inf, label = "*", vjust=2, hjust=2, family="Times")

slpc3 <- predframe_ash %>% 
  filter(fshape == "PC3") %>% 
  mutate(nSH = shell_length/mean(shell_length)) %>% 
  ggplot() + geom_smooth(aes(nSH, preds, color=..y..), lwd = 0.5, se = F) +
  geom_smooth(aes(nSH, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(nSH, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  theme_science() +
  theme(axis.title.y = element_blank(),
        legend.position = "none") +
  scale_colour_gradient(low = maxp[2], high = maxp[9]) +
  scale_y_continuous(labels = scaleFUN, expand = c(0.1,0.1)) +
  xlab("nSH") + ylab("shape-PC3") +
  annotate("text", x=Inf, y = Inf, label = "***", vjust=2, hjust=2, family="Times")

# Combine plots -----------------------------------------------------------

pdf("gamm_plot.pdf", height = 5, onefile = F)

# ((PC1 | tpc1 | apc1 | spc1) /
# (PC2 | tpc2 | apc2 | spc2) /
# (PC3 | tpc3 | apc3 | spc3) /
# (PC4 | tpc4 | apc4 | spc4) /
# (PC5 | tpc5 | apc5 | spc5)) +
# plot_layout(widths = c(.5,3,3,3,
#                        .5,3,3,3,
#                        .5,3,3,3,
#                        .5,3,3,3,
#                        .5,3,3,3),
#             heights = c(1,1,1,1,1.05))

((tpc1 | apc1 | spc1 | slpc1) /
    (tpc2 | apc2 | spc2 | slpc2) /
    (tpc3 | apc3 | spc3 | slpc3)) +
  plot_layout(widths = c(.5,3,3,3,
                         .5,3,3,3,
                         .5,3,3,3),
              heights = c(1,1,1.05))

grid.text("A", x = 0.015, y = 0.98, gp=gpar(fontsize=10, fontfamily = "Times"))
grid.text("B", x = 0.015, y = 0.68, gp=gpar(fontsize=10, fontfamily = "Times"))
grid.text("C", x = 0.015, y = 0.34, gp=gpar(fontsize=10, fontfamily = "Times"))

dev.off()
