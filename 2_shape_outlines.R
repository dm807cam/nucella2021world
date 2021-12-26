library("tidyverse")
library("Momocs")
library("cowplot")
library("viridis")

rm(list=ls())

# Import and prepare shape objects ----------------------------------------
# Load names of images into R
# First we have to sort the file list because R thinks that 10 is closer to 1 than 2. Silly R :)
# filelist <- tibble(path = list.files(path = "data/shell_images",
#                                      pattern = "\\.jpg",
#                                      full.names=TRUE)) %>%
#   separate(path,sep="-", into=c("id1", "id2")) %>%
#   mutate(sampleID = as.numeric(str_extract_all(id1, "[[:digit:]]+",
#                                                simplify = TRUE))) %>%
#   mutate(sampleNo = as.numeric(str_extract_all(id2, "[[:digit:]]+",
#                                                simplify = TRUE))) %>%
#   arrange(sampleID, sampleNo) %>%
#   mutate(path = paste(id1,id2, sep="-")) %>%
#   pull(path)
# 
# # Import images into R object
# nucella_xy <- import_jpg(filelist,
#                          auto.notcentered = TRUE,
#                          fun.notcentered = NULL,
#                          threshold = 0.5)
# 
# saveRDS(nucella_xy, "data/shell_xy_object.rds")

# Read xy-file
nucella_xy <- readRDS("data/shell_xy_object.rds")

# Auto align shapes along longitudinal axis.
rotate_shape <- function(x) {
  
  x <- coo_smooth(x, 100)
  
  dobj <- raster::pointDistance(x, lonlat=F)
  max_dist <- which(dobj==max(dobj), arr.ind=T)
  
  max_line <- rbind(x[max_dist[1],], x[max_dist[2],])
  
  mod <- lm(max_line[,2]~max_line[,1])
  
  xcorr <- coo_rotate(x, theta = -coef(mod)[2])
  
  return(xcorr)
}

nucella_xy_adj <- list()
pb <- txtProgressBar(0, length(nucella_xy), style = 3)
for(ii in 1:length(nucella_xy)) {
  
  setTxtProgressBar(pb, ii)
  
  ij <- names(nucella_xy)[ii]
  x <- nucella_xy[[ii]]
  tmp_out <- rotate_shape(x)

  nucella_xy_adj[[ij]] <- tmp_out
}

# Save nucella_xy_adj
saveRDS(nucella_xy_adj, "data/shell_xy_adj_object.rds")

# Read xy_adj-file
# nucella_xy_adj <- readRDS("data/shell_xy_adj_object.rds")

nucella_coo <- Out(nucella_xy_adj) 

# Save nucella_coo
saveRDS(nucella_coo, "data/shell_coo_object.rds")

# Read coo-file
# nucella_coo <- readRDS("data/shell_coo_object.rds")

# Clean shapes and align --------------------------------------------------
# Outline smoothing with 5 smoothing iterations
nucella_coo_pr <- nucella_coo %>% 
        # Center pseudo-landmarks
        coo_center() %>% 
        # Scale pseudo-landmarks
        coo_scale() %>%  
        # Sample pseudo-landmarks with 3000 points
        coo_sample(1000) %>% 
        # Normalization of the starting point
        coo_slidedirection("right") %>% 
        # Rotate to apex top
        coo_rotatecenter(-1.5708) # 90 degrees

# Save nucella_coo_pr
saveRDS(nucella_coo_pr, "data/shell_coo_pr_object.rds")

# Export smoothed forms to check the accuracy of outlines tracing
for(ii in 1:length(nucella_coo_pr)) {
   img_name <- rlist::list.names(nucella_coo_pr)[ii]
   
   png(paste("data/shell_images/outlines/",img_name,".png",sep=""))
   Momocs::coo_plot(nucella_coo_pr[ii])
   dev.off()
}

# Plot normalized outlines
stack(nucella_coo_pr)

# Elliptic Fourier Analysis -----------------------------------------------
# Identify number of harmonics required for the efourier reconstruction
# Shape reconstruction using an increasing number of harmonics 1-12
cal_shape <- calibrate_reconstructions_efourier(nucella_coo_pr, 
                                                range = 1:12)
plot(cal_shape) + 
        theme_bw() + 
        scale_fill_gradientn()

# Deviation between the best possible fit and a given number of harmonics
cal_dev <- calibrate_deviations_efourier(nucella_coo_pr, 
                                         id=1:20, 
                                         range=c(6,8,10,12))$gg
plot(cal_dev) + 
        theme_classic()

# Harmonic power
cal_harm <- calibrate_harmonicpower_efourier(nucella_coo_pr, nb.h = 13)

# Show harmonic power
cal_harm$minh # 10 harmonics gather 99% of the total harmonic power
cal_harm$gg$data$harm <- factor(cal_harm$gg$data$harm, 
                                levels=c("h1","h2","h3","h4","h5",
                                         "h6","h7","h8","h9","h10","h12"))

mo_pan_0 <- plot(cal_harm$gg) + 
        theme_bw() + 
        geom_hline(yintercept=100, colour="orange") + 
        theme(text = element_text(size = 11))
title_mo_pan_0 <- ggdraw() + 
        draw_label("Fourier power spectrum", fontface='italic', size = 10)
plot_grid(title_mo_pan_0, mo_pan_0, ncol=1, rel_heights=c(0.1, 1))

# Elliptic Fourier analysis of shell shapes 
nucella_coe <- efourier(nucella_coo_pr, nb.h = 10, norm=F)

# Save nucella_coe
saveRDS(nucella_coe, "data/shell_coe_object.rds")

# Read nucella_coe
# nucella_coe <- readRDS("data/shell_coe_object.rds")

# Overview of harmonic contributions
hc <- hcontrib(nucella_coe, main=NULL, harm.r=1:10, amp.r = c(0,1,2,3), id=1, col="grey20")

# Effect of all 9 harmonics on the shape reconstruction. Every harmonic is 
# represented when their corresponding coefficients are multiplied by an
# amplification factor that illustrates their removal [0], the normal shape [1]
# and exaggerated shapes [2+].

# PCA analysis of shape objects -------------------------------------------
# Calculate PCA from shapes
nucella_pca <- PCA(nucella_coe)

# The first 3 PCA describe ~81% and the first 5 PCA describe ~90% of the shape 
# variability among specimens

# Plot variability covered/exlained by each PC
# Boxplot function broken due to an error with ggplot - extract and plot data manually
df_nucella_pca <- gather(as.data.frame(nucella_pca$x), PC, values, PC1:PC12, factor_key=TRUE)

p1 <- ggplot(subset(df_nucella_pca), aes(PC, values)) + geom_boxplot() + theme_classic() 
p2 <- scree_plot(nucella_pca, nax = 1:12) + theme_classic() 
plot_grid(p1,p2, ncol=2)

# PCA contributions --------------------------------------------------------
PC_contrib <- PCcontrib(nucella_pca, nax=1:6, sd.r = c(-3,0,3))

# Plot PC contribution 
PC_contrib$gg + geom_polygon(aes(fill=shp), colour="black", size=1) + 
        theme_bw() + 
        scale_fill_viridis() +
        ylab("PCs") +
        xlab(expression(paste("Mean +- 3",sigma,sep="")))

# Export -------------------------------------
# Import explanatory variables
nucella_fac <- tibble(read_delim("data/sample_list.csv",delim=",")) %>% 
  arrange(as.numeric(sampleID), 
          as.numeric(sampleNo))

# Export first 6 PCs which cover 90% of possible shell shape variance.
nucella_PCs <- data.frame(PC1 = nucella_pca$x[, 1],
                          PC2 = nucella_pca$x[, 2],
                          PC3 = nucella_pca$x[, 3],
                          PC4 = nucella_pca$x[, 4],
                          PC5 = nucella_pca$x[, 5]) %>% 
                rownames_to_column(var = "picname") %>% 
                left_join(nucella_fac) %>% 
                select(-year, -folder, -picname, -country)

# Add shape length from images a PCA data frame
nucella_PCs$shell_length <- coo_length(nucella_coo)

# Export to csv
write.csv(nucella_PCs, "data/shell_outline_PCA.csv", row.names = FALSE)

