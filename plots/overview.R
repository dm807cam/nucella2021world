library("tidyverse")
library("rnaturalearth")
#devtools::install_github("ropensci/rnaturalearthhires")
library("cmocean")
library("RColorBrewer")
library("latex2exp")
library("ggrepel")
library("Momocs")
library("patchwork")
library("grid")

rm(list=ls())

# Import helper file. This file contain proprietary code and
# will not be supplied with the rest of the code.
# All functions used from this file will be highlighted in the script.
source("../imp_func.R", local = TRUE)

# Import world map as sf
world <- ne_countries(scale = 10, returnclass = "sf")

# Import environmental data
df_env <- read_csv("../data/environmental_comp.csv")
df_carb <- read_csv("../data/environmental_all.csv")

# Plot thermal sampling map
p1 <- ggplot(data = world) +
  geom_raster(data = df_carb,
              aes(lon, lat, fill = T),
              show.legend = F, alpha=0.75) +
  geom_sf(fill = "white", size = 0.2) +
  scale_fill_cmocean(name = "thermal") +
  geom_text_repel(
    df_env,
    mapping = aes(lon, lat, label = nlocation),
    force = 30, nudge_x = -6, fontface = 1, 
    min.segment.length = unit(0, 'lines'),
    family = "Times", size = 4, segment.size= 0.1) +
  theme_science() +
  theme(text = element_text(family ="Times"),
        title =element_text(family ="Times")) +
  coord_sf(xlim = c(-35, 24),
           ylim = c(36, 72),
           expand = F) +
  xlab("longitude (degree)") +
  ylab("latitude (degree)")

# Select variables to plot
dat <- c("Temp", "Psal", "pH", "OmegaCalcite")

# Create variable names in Latex
datname <-
  c("T $\\degree$C",
    "S (psu) ",
    "pH",
    "$\\Omega_{Ca}$")

# Define color scheme
col <- c("thermal", "haline", "speed", "dense")

# Create line plots
datalist <- list()
for (ii in 1:length(dat)) {
  
  is_even <- function(x) x %% 2 == 0  
  ij <- unique(dat)[ii]
  
  tmp <- ggplot(df_env, aes_string(ij, "lat", fill = ij)) +
    geom_path(linetype=2, colour="grey60") +
    geom_point(pch=21, size=2, colour="grey20") +
    scale_fill_cmocean(name = col[ii]) +
    scale_y_continuous(limits = c(36, 72),
                       breaks=seq(36,72,by=5),
                       expand = c(0,0)) +
    {if(!(is_even(ii))){scale_x_continuous(n.breaks = 3, expand = c(0.1,0))}} +
    {if(is_even(ii)){scale_x_continuous(n.breaks = 3, expand = c(0.1,0), position = "top")}}+
    xlab(unname(TeX(datname[ii])))
  
  datalist[[ii]]<- tmp + theme_set(theme_classic(base_size = 10, base_family = 'Times')) +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_line(size=0.3),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) 
  
}

# patchwork
# Export plots to PDF
pdf("overview.pdf", height = 4,
    onefile = F)

p1 | (datalist[[1]] | datalist[[2]] | datalist[[3]] | datalist[[4]])

dev.off()
