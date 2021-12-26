library("tidyverse")
library("patchwork")
library("cmocean")

rm(list=ls())

# Import helper file. This file contain proprietary code and
# will not be supplied with the rest of the code.
# All functions used from this file will be highlighted in the script.
source("../imp_func.R", local = TRUE)

# Import caliper data
df_pca_long <- read_csv("../data/df_pca_long.csv")
df_caliper_long <- read_csv("../data/df_caliper_long.csv")


df_apt <- df_caliper_long %>% 
  filter(fshape =="aperture_size") %>%
  mutate(svalue = range01(value)) %>% 
  group_by(flocation) %>% 
  summarise(avg_apt = mean(svalue),
            stdev_apt = sd(svalue))

df_pc1 <- df_pca_long %>% 
  filter(fshape =="PC1") %>% 
  mutate(svalue = range01(value)) %>% 
  group_by(flocation) %>% 
  summarise(avg_pc = mean(svalue),
            stdev_pc = sd(svalue))

df_pc2 <- df_pca_long %>% 
  filter(fshape =="PC2") %>% 
  mutate(svalue = range01(value)) %>% 
  group_by(flocation) %>% 
  summarise(avg_pc = mean(svalue),
            stdev_pc = sd(svalue))

df_pc3 <- df_pca_long %>% 
  filter(fshape =="PC3") %>% 
  mutate(svalue = range01(value)) %>% 
  group_by(flocation) %>% 
  summarise(avg_pc = mean(svalue),
            stdev_pc = sd(svalue))

df_pc4 <- df_pca_long %>% 
  filter(fshape =="PC4") %>% 
  mutate(svalue = range01(value)) %>% 
  group_by(flocation) %>% 
  summarise(avg_pc = mean(svalue),
            stdev_pc = sd(svalue))

df_pc5 <- df_pca_long %>% 
  filter(fshape =="PC5") %>% 
  mutate(svalue = range01(value)) %>% 
  group_by(flocation) %>% 
  summarise(avg_pc = mean(svalue),
            stdev_pc = sd(svalue))


p1 <- df_apt %>% 
full_join(df_pc1) %>% 
  ggplot(aes(avg_pc, avg_apt)) +
  geom_abline(slope=1, col="orange", linetype=2) +
  geom_errorbar(aes(ymin = avg_apt - stdev_apt, ymax = avg_apt + stdev_apt), colour="grey80") +
  geom_errorbarh(aes(xmin = avg_pc - stdev_pc, xmax = avg_pc + stdev_pc), colour="grey80") +
    geom_point(pch=21, size=2, colour="grey20", fill="white") +
    scale_x_continuous(expand = c(.01,.01),
                       labels=function(x) paste0(x,"°N"), 
                       name = "shape-PC1") +
    scale_y_continuous(expand = c(.1,.1)) +
    ylab("Aperture size") + 
    theme_science() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(), 
          aspect.ratio = 1) 

p2 <- df_apt %>% 
  full_join(df_pc2) %>% 
  ggplot(aes(avg_pc, avg_apt)) +
  geom_abline(slope=1, col="orange", linetype=2) +
  geom_errorbar(aes(ymin = avg_apt - stdev_apt, ymax = avg_apt + stdev_apt), colour="grey80") +
  geom_errorbarh(aes(xmin = avg_pc - stdev_pc, xmax = avg_pc + stdev_pc), colour="grey80") +
  geom_point(pch=21, size=2, colour="grey20", fill="white") +
  scale_x_continuous(expand = c(.01,.01),
                     labels=function(x) paste0(x,"°N"), 
                     name = "shape-PC2") +
  scale_y_continuous(expand = c(.1,.1)) +
  ylab("Aperture size") + 
  theme_science() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title.y = element_blank(),
        aspect.ratio = 1) 


p3 <- df_apt %>% 
  full_join(df_pc3) %>% 
  ggplot(aes(avg_pc, avg_apt)) +
  geom_abline(slope=1, col="orange", linetype=2) +
  geom_errorbar(aes(ymin = avg_apt - stdev_apt, ymax = avg_apt + stdev_apt), colour="grey80") +
  geom_errorbarh(aes(xmin = avg_pc - stdev_pc, xmax = avg_pc + stdev_pc), colour="grey80") +
  geom_point(pch=21, size=2, colour="grey20", fill="white") +
  scale_x_continuous(expand = c(.01,.01),
                     labels=function(x) paste0(x,"°N"), 
                     name = "shape-PC3") +
  scale_y_continuous(expand = c(.1,.1)) +
  ylab("Aperture size") + 
  theme_science() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title.y = element_blank(),
        aspect.ratio = 1) 

# patchwork
# Export plots to PDF
pdf("PCA_apt_size.pdf", height = 2,
    onefile = F)

(p1 | p2 | p3)

dev.off()
