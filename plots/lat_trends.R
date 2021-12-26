library("tidyverse")
library("patchwork")
library("latex2exp")
library("cmocean")
library("Momocs")
library("patchwork")
library("grid")

rm(list=ls())

# Import helper file. This file contain proprietary code and
# will not be supplied with the rest of the code.
# All functions used from this file will be highlighted in the script.
source("../imp_func.R", local = TRUE)

# Import caliper data
df_pca_long <- read_csv("../data/df_pca_long.csv")
df_caliper_long <- read_csv("../data/df_caliper_long.csv")


# Import shell shape file
nucella_coe <- readRDS("../data/shell_coe_object.rds")

# Calculate PCA from shapes
nucella_pca <- PCA(nucella_coe)

# Plot PC contribution 
PC_contrib <- PCcontrib(nucella_pca, nax=1:5, sd.r = c(-3,0,3))

# Extract Shell shapes from PCA
PC1_1 <- as.matrix(PC_contrib$shp$`1.1`[1:2])
PC1_3 <- as.matrix(PC_contrib$shp$`3.1`[1:2])
PC2_1 <- as.matrix(PC_contrib$shp$`1.2`[1:2])
PC2_3 <- as.matrix(PC_contrib$shp$`3.2`[1:2])
PC3_1 <- as.matrix(PC_contrib$shp$`1.3`[1:2])
PC3_3 <- as.matrix(PC_contrib$shp$`3.3`[1:2])
PC4_1 <- as.matrix(PC_contrib$shp$`1.4`[1:2])
PC4_3 <- as.matrix(PC_contrib$shp$`3.4`[1:2])
PC5_1 <- as.matrix(PC_contrib$shp$`1.5`[1:2])
PC5_3 <- as.matrix(PC_contrib$shp$`3.5`[1:2])

isopalette <- colorRampPalette(rev(brewer.pal(10, "RdBu")))
maxp <- rev(brewer.pal(11, "RdBu"))

PC1 <- wrap_elements(full = ~ tps_grid(PC1_1, PC1_3, shp.lwd = c(1.5, 1.5), 
                                       shp.border = c(maxp[2], maxp[9]),
                                       amp = 0,
                                       legend = F))
PC2 <- wrap_elements(full = ~ tps_grid(PC2_1, PC2_3, shp.lwd = c(1.5, 1.5), 
                                       shp.border = c(maxp[2], maxp[9]),
                                       amp = 0, 
                                       legend = F))
PC3 <- wrap_elements(full = ~ tps_grid(PC3_1, PC3_3, shp.lwd = c(1.5, 1.5), 
                                       shp.border = c(maxp[2], maxp[9]),
                                       amp = 0,
                                       legend = F))
PC4 <- wrap_elements(full = ~ tps_grid(PC4_1, PC4_3, shp.lwd = c(1.5, 1.5), 
                                       shp.border = c(maxp[2], maxp[9]),
                                       amp = 0.,
                                       legend = F))

PC5 <- wrap_elements(full = ~ tps_grid(PC5_1, PC5_3, shp.lwd = c(1.5, 1.5), 
                                       shp.border = c(maxp[2], maxp[9]),
                                       amp = 0,
                                       legend = F))

p1 <- df_pca_long %>% 
  filter(fshape == "PC1") %>% 
  group_by(flocation) %>% 
  summarise(avg = mean(value),
            stdev = sd(value),
            lat = lat[1]) %>% 
  arrange(lat) %>% 
  ggplot(aes(lat, avg)) +
  geom_errorbar(aes(ymin=avg-stdev, ymax=avg+stdev), colour="grey80") +
    geom_path(linetype=2, colour="grey60") +
    geom_point(pch=21, size=2, colour="grey20", fill="white") +
    geom_smooth(method = "lm", colour="orange", alpha=.5, size=0.5, se=F) +
    scale_x_continuous(limits = c(36, 72),
                       breaks=seq(36,72,by=5),
                       expand = c(0,0),
                       labels=function(x) paste0(x,"°N"), 
                       name = "latitude (degree)") +
    scale_y_continuous(expand = c(0.1,0)) +
    ylab("shape-PC1") + 
    theme_set(theme_classic(base_size = 10, base_family = 'Times')) +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"),
          axis.line = element_line(size=0.3),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
  coord_flip() + 
  annotate("text",  x=Inf, y = Inf, label = "***", vjust=1, hjust=1, family='Times')

p2 <- df_pca_long %>% 
  filter(fshape == "PC2") %>% 
  group_by(flocation) %>% 
  summarise(avg = mean(value),
            stdev = sd(value),
            lat = lat[1]) %>% 
  arrange(lat) %>% 
  ggplot(aes(lat, avg)) +
  geom_errorbar(aes(ymin=avg-stdev, ymax=avg+stdev), colour="grey80") +
  geom_path(linetype=2, colour="grey60") +
  geom_point(pch=21, size=2, colour="grey20", fill="white") +
  geom_smooth(method = "lm", colour="orange", alpha=.5, size=0.5, se=F) +
  scale_x_continuous(limits = c(36, 72),
                     breaks=seq(36,72,by=5),
                     expand = c(0,0),
                     labels=function(x) paste0(x,"°N"), 
                     name = "latitude (degree)") +
  scale_y_continuous(expand = c(0.1,0)) +
  ylab("shape-PC2") + 
  theme_set(theme_classic(base_size = 10, base_family = 'Times')) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"),
        axis.line = element_line(size=0.3),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() + 
  annotate("text",  x=Inf, y = Inf, label = "***", vjust=1, hjust=1, family='Times')

p3 <- df_pca_long %>% 
  filter(fshape == "PC3") %>% 
  group_by(flocation) %>% 
  summarise(avg = mean(value),
            stdev = sd(value),
            lat = lat[1]) %>% 
  arrange(lat) %>% 
  ggplot(aes(lat, avg)) +
  geom_errorbar(aes(ymin=avg-stdev, ymax=avg+stdev), colour="grey80") +
  geom_path(linetype=2, colour="grey60") +
  geom_point(pch=21, size=2, colour="grey20", fill="white") +
  geom_smooth(method = "lm", colour="orange", alpha=.5, size=0.5, se=F) +
  scale_x_continuous(limits = c(36, 72),
                     breaks=seq(36,72,by=5),
                     expand = c(0,0),
                     labels=function(x) paste0(x,"°N"), 
                     name = "latitude (degree)") +
  scale_y_continuous(expand = c(0.1,0)) +
  ylab("shape-PC3") + 
  theme_set(theme_classic(base_size = 10, base_family = 'Times')) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"),
        axis.line = element_line(size=0.3),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() + 
  annotate("text",  x=Inf, y = Inf, label = "**", vjust=1, hjust=1, family='Times')

p4 <- df_pca_long %>% 
  filter(fshape == "PC4") %>% 
  group_by(flocation) %>% 
  summarise(avg = mean(value),
            stdev = sd(value),
            lat = lat[1]) %>% 
  arrange(lat) %>% 
  ggplot(aes(lat, avg)) +
  geom_errorbar(aes(ymin=avg-stdev, ymax=avg+stdev), colour="grey80") +
  geom_path(linetype=2, colour="grey60") +
  geom_point(pch=21, size=2, colour="grey20", fill="white") +
  geom_smooth(method = "lm", colour="orange", alpha=.5, size=0.5, se=F) +
  scale_x_continuous(limits = c(36, 72),
                     breaks=seq(36,72,by=5),
                     expand = c(0,0),
                     labels=function(x) paste0(x,"°N"), 
                     name = "latitude (degree)") +
  scale_y_continuous(expand = c(0.1,0)) +
  ylab("shape-PC4") + 
  theme_set(theme_classic(base_size = 10, base_family = 'Times')) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"),
        axis.line = element_line(size=0.3),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() + 
  annotate("text",  x=Inf, y = Inf, label = "n.s.", vjust=1, hjust=1, family='Times')

p5 <- df_pca_long %>% 
  filter(fshape == "PC5") %>% 
  group_by(flocation) %>% 
  summarise(avg = mean(value),
            stdev = sd(value),
            lat = lat[1]) %>% 
  arrange(lat) %>% 
  ggplot(aes(lat, avg)) +
  geom_errorbar(aes(ymin=avg-stdev, ymax=avg+stdev), colour="grey80") +
  geom_path(linetype=2, colour="grey60") +
  geom_point(pch=21, size=2, colour="grey20", fill="white") +
  geom_smooth(method = "lm", colour="orange", alpha=.5, size=0.5, se=F) +
  scale_x_continuous(limits = c(36, 72),
                     breaks=seq(36,72,by=5),
                     expand = c(0,0),
                     labels=function(x) paste0(x,"°N"), 
                     name = "latitude (degree)") +
  scale_y_continuous(expand = c(0.1,0)) +
  ylab("shape-PC5") + 
  theme_set(theme_classic(base_size = 10, base_family = 'Times')) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"),
        axis.line = element_line(size=0.3),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() + 
  annotate("text",  x=Inf, y = Inf, label = "n.s.", vjust=1, hjust=1, family='Times')

PC1_perc <- round(scree(nucella_pca)$cumsum[1] * 100, 2)
PC2_perc <- round((scree(nucella_pca)$cumsum[2] - scree(nucella_pca)$cumsum[1]) * 100,2)
PC3_perc <- round((scree(nucella_pca)$cumsum[3] - scree(nucella_pca)$cumsum[2]) * 100,2)
PC4_perc <- round((scree(nucella_pca)$cumsum[4] - scree(nucella_pca)$cumsum[3]) * 100,2)
PC5_perc <- round((scree(nucella_pca)$cumsum[5] - scree(nucella_pca)$cumsum[4]) * 100,2)

# patchwork
# Export plots to PDF
pdf("lat_trends.pdf",height = 5,
    onefile = F)

(p1 | p2 | p3 | p4 | p5) /
  (PC1 | PC2 | PC3 | PC4 | PC5) + plot_layout(heights = c(2,1))

grid.text("A", x = 0.015, y = 0.98, gp=gpar(fontsize=10, fontfamily = "Times"))
grid.text("B", x = 0.015, y = 0.3, gp=gpar(fontsize=10, fontfamily = "Times"))

grid.text(paste(PC1_perc,"%", sep=""), x = 0.135, y = 0.17, gp=gpar(fontsize=9, fontfamily = "Times"))
grid.text(paste(PC2_perc,"%", sep=""), x = 0.355, y = 0.17, gp=gpar(fontsize=9, fontfamily = "Times"))
grid.text(paste(PC3_perc,"%", sep=""), x = 0.54, y = 0.17, gp=gpar(fontsize=9, fontfamily = "Times"))
grid.text(paste(PC4_perc,"%", sep=""), x = 0.72, y = 0.17, gp=gpar(fontsize=9, fontfamily = "Times"))
grid.text(paste(PC5_perc,"%", sep=""), x = 0.9, y = 0.17, gp=gpar(fontsize=9, fontfamily = "Times"))

dev.off()
