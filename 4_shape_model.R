library("tidyverse")
library("mgcv")
library("gratia")
library("gstat")
library("sp")
library("MuMIn")

rm(list=ls())

# Import helper file. This file contain proprietary code and
# will not be supplied with the rest of the code.
# All functions used from this file will be highlighted in the script.
source("imp_func.R", local = TRUE)

# Import data -------------------------------------------------------------
# Read in files
df_pca_long <- read_csv("data/df_pca_long.csv")

# Check files
glimpse(df_pca_long)

# Make sure location is treated as a factor
# Select shape descriptors that show significant latitudinal trends
df_pca_long <- df_pca_long %>% 
  mutate(fshape = as.factor(fshape),
         flocation = as.factor(flocation)) %>% 
  group_by(fshape) %>% 
  mutate(svalue = range01(value)+1) %>% 
  ungroup() %>% 
  filter(fshape != "PC4") %>% 
  filter(fshape != "PC5") %>% 
  glimpse()

# Simple function to calculate maximal VIF observed between predictors
# This function is based on the CORVIF function by Zuur et al.
max_corr <- function(x) {
  
  # Check if model is of class gam
  if(class(x)[1] != "gam"){
    print("Error: Invalid model class")}
  
  # Select model table as dataframe and drop dependent variable
  corr <- data.frame(x$model) %>% 
    select(-1)
  
  # Select only numeric terms and calculate VIF
  corr <-  corr %>% 
    select_if(is.numeric) %>% 
    # Use CROVIF function from the helper file
    corvif() %>% 
    round(2)
  
  # Vif table as matrix
  corrm <- as.matrix(corr)
  
  # Find mxx VIF value 
  if (length(corrm)<=1){
    corrm <- 0
    max(abs(corrm))
  } else {
    max(abs(corrm))
  }  
}

# Fit full shell shape-PCs Model ------------------------------
m_gam <- gam(svalue ~ fshape + 
                      s(Temp, bs="cr", k = 3, by = fshape) + 
                      s(Psal, bs="cr", k = 3, by = fshape) + 
                      s(pH, bs="cr", k = 3, by = fshape) + 
                      s(ALK, bs="cr", k = 3, by = fshape) + 
                      s(SIR, bs="cr", k = 3, by = fshape) +
                      s(OmegaCalcite, bs="cr", k = 3, by = fshape) +
                      s(OmegaAragonite, bs="cr", k = 3, by = fshape) +
                      s(shell_length, bs="cr", k = 3, by = fshape) +
                      s(flocation, bs="re"),
                   method = "ML", 
                   family=gaussian(link="identity"),
                   data = df_pca_long, 
                   control = gam.control(trace = F))

# Model selection accounting for VIF factor among predictors
options(na.action = na.fail)
d_tab <- dredge(m_gam, rank = "AIC", m.lim = c(2, 5), extra = c(max_corr), trace = 2) 
# Retrieve best model from models with VIF < 3
m_gam <- get.models(d_tab, subset = max_corr < 3)[[1]] 

# Check summary
summary(m_gam)

# Fit preferred model with REML and validate ------------------------------
mody <- gam(svalue ~ fshape +
              s(Temp, bs="cr", k = 3, by=fshape) + 
              s(ALK, bs="cr", k = 3, by=fshape) + 
              s(SIR, bs="cr", k = 3, by=fshape) +
              s(shell_length, bs="cr", k = 3, by=fshape) +
              s(flocation, bs="re"),
            method = "REML", 
            family= gaussian(link="identity"),
            data = df_pca_long, 
            control = gam.control(trace = T))

summary(mody)
appraise(mody)

# Checking model assumptions
R <- resid(mody)
E <- fitted(mody)

# Normality
n1 <- data.frame(E,R) %>% 
  ggplot(aes(sample=R)) +
  geom_qq(alpha=0.25) +
  geom_qq_line(colour="red")

n2 <- data.frame(E,R) %>% 
  ggplot(aes(R)) +
  geom_histogram(bins = 50, colour="grey20", fill="lightblue") 

n1 + n2

# Homoscedasticity
data.frame(E,R) %>% 
ggplot(aes(E,R)) +
  geom_point(alpha=0.25) +
  geom_hline(yintercept = 0, colour="darkblue") +
  geom_smooth(method="loess", colour="red", linetype=2, span=1, se=F)

# Spatial autocorrelation
kachy <-data.frame(df_pca_long$lon,df_pca_long$lat,R)
colnames(kachy) <- c("lon", "lat", "resids")
coordinates(kachy) <- c('lon','lat')

vario <- variogram(R~1, data=kachy, alpha=c(0,45,90,135))
plot(vario)
# No spatial autocorrelation found. 

# Residuals vs Temperature per level
lattice::xyplot(R ~ Temp | fshape, data = df_pca_long, abline = 0, col = 1, cex = 0.5)
lattice::xyplot(R ~ ALK | fshape, data = df_pca_long, abline = 0, col = 1, cex = 0.5)
lattice::xyplot(R ~ SIR | fshape, data = df_pca_long, abline = 0, col = 1, cex = 0.5)
lattice::xyplot(R ~ shell_length | fshape, data = df_pca_long, abline = 0, col = 1, cex = 0.5)

# Predict model Temperature -----------------------------------------------
pdat <- expand.grid(Temp = seq(min(df_pca_long$Temp), max(df_pca_long$Temp),
                           by = ((max(df_pca_long$Temp) - min(df_pca_long$Temp))/(3000 - 1))),
                    ALK = mean(df_pca_long$ALK),
                    SIR = mean(df_pca_long$SIR),
                    shell_length = mean(df_pca_long$shell_length))

pdat_temp <- data.frame(Temp = rep(pdat$Temp, 3), 
                        ALK = rep(pdat$ALK, 3), 
                        SIR = rep(pdat$SIR, 3), 
                        shell_length = rep(pdat$shell_length, 3),
                        flocation = "Voe",
                        fshape = rep(c("PC1", "PC2", "PC3"), each = 3000))

pred_temp <- predict(mody, newdata = pdat_temp, exclude = c("s(flocation)"), 
                     type = "response", se.fit = T)

predframe_temp <- data.frame(Temp = pdat_temp$Temp,
                             ALK = pdat_temp$ALK, 
                             SIR = pdat_temp$SIR, 
                             shell_length = pdat_temp$shell_length,
                             fshape = pdat_temp$fshape, 
                             preds = pred_temp$fit, se = pred_temp$se.fit)

# function for having 1 decimal values
scaleFUN <- function(x) sprintf("%.1f", x)
ggplot(predframe_temp) + geom_smooth(aes(Temp, preds), col = 1, lwd = 0.5, se = F) +
  geom_smooth(aes(Temp, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(Temp, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  facet_wrap(~fshape, nrow = 2, scale = "free") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + 
  scale_y_continuous(labels = scaleFUN) +
  xlab("Temperature") + ylab("shape-PCs")

# Predict model Alkalinity --------------------------------------------------
pdat <- expand.grid(ALK = seq(min(df_pca_long$ALK), max(df_pca_long$ALK), 
                          by = ((max(df_pca_long$ALK) - min(df_pca_long$ALK))/(3000 - 1))),
                    Temp = mean(df_pca_long$Temp),
                    SIR = mean(df_pca_long$SIR),
                    shell_length = mean(df_pca_long$shell_length))

pdat_alk <- data.frame(Temp = rep(pdat$Temp, 3), 
                        ALK = rep(pdat$ALK, 3), 
                        SIR = rep(pdat$SIR, 3), 
                        shell_length = rep(pdat$shell_length, 3), 
                        flocation = "Voe",
                        fshape = rep(c("PC1", "PC2", "PC3"), each = 3000))

pred_alk <- predict(mody, newdata = pdat_alk, exclude = c("s(flocation)"), 
                    type = "response", se.fit = T)

predframe_alk <- data.frame(Temp = pdat_alk$Temp,
                             ALK = pdat_alk$ALK, 
                             SIR = pdat_alk$SIR, 
                             shell_length = pdat_alk$shell_length,
                             fshape = pdat_alk$fshape, 
                             preds = pred_alk$fit, se = pred_alk$se.fit)

ggplot(predframe_alk) + geom_smooth(aes(ALK, preds), col = 1, lwd = 0.5, se = F) +
  geom_smooth(aes(ALK, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(ALK, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  facet_wrap(~fshape, nrow = 2, scale = "free") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + 
  scale_y_continuous(labels = scaleFUN) +
  xlab("Salinity") + ylab("shape-PCs")

# Predict model SIR --------------------------------------------------------
pdat <- expand.grid(SIR = seq(min(df_pca_long$SIR), max(df_pca_long$SIR), 
                          by = ((max(df_pca_long$SIR) - min(df_pca_long$SIR))/(3000 - 1))),
                    Temp = mean(df_pca_long$Temp),
                    ALK = mean(df_pca_long$ALK),
                    shell_length = mean(df_pca_long$shell_length))

pdat_SIR <- data.frame(Temp = rep(pdat$Temp, 3), 
                        ALK = rep(pdat$ALK, 3), 
                        SIR = rep(pdat$SIR, 3), 
                        shell_length = rep(pdat$shell_length, 3), 
                        flocation = "Voe",
                        fshape = rep(c("PC1", "PC2", "PC3"), each = 3000))

pred_SIR <- predict(mody, newdata = pdat_SIR, exclude = c("s(flocation)"), 
                    type = "response", se.fit = T)

predframe_SIR <- data.frame(Temp = pdat_SIR$Temp,
                             ALK = pdat_SIR$ALK, 
                             SIR = pdat_SIR$SIR, 
                             shell_length = pdat_SIR$shell_length,
                             fshape = pdat_SIR$fshape, 
                             preds = pred_SIR$fit, se = pred_SIR$se.fit)

ggplot(predframe_SIR) + geom_smooth(aes(SIR, preds), col = 1, lwd = 0.5, se = F) +
  geom_smooth(aes(SIR, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(SIR, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  facet_wrap(~fshape, nrow = 2, scale = "free") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + 
  scale_y_continuous(labels = scaleFUN) +
  xlab("Salinity") + ylab("shape-PCs")

# Predict model shell length ----------------------------------------------
pdat <- expand.grid(shell_length = seq(min(df_pca_long$shell_length), max(df_pca_long$shell_length), 
                                   by = ((max(df_pca_long$shell_length) - min(df_pca_long$shell_length))/(3000 - 1))),
                    Temp = mean(df_pca_long$Temp),
                    ALK = mean(df_pca_long$ALK),
                    SIR = mean(df_pca_long$SIR))

pdat_ash <- data.frame(Temp = rep(pdat$Temp, 3), 
                       ALK = rep(pdat$ALK, 3), 
                       SIR = rep(pdat$SIR, 3), 
                       shell_length = rep(pdat$shell_length, 3), 
                       flocation = "Voe",
                       fshape = rep(c("PC1", "PC2", "PC3"), each = 3000))

pred_ash <- predict(mody, newdata = pdat_ash, exclude = c("s(flocation)"), 
                    type = "response", se.fit = T)

predframe_ash <- data.frame(Temp = pdat_ash$Temp,
                            ALK = pdat_ash$ALK, 
                            SIR = pdat_ash$SIR, 
                            shell_length = pdat_ash$shell_length,
                            fshape = pdat_ash$fshape, 
                            preds = pred_ash$fit, se = pred_ash$se.fit)

ggplot(predframe_ash) + geom_smooth(aes(shell_length, preds), col = 1, lwd = 0.5, se = F) +
  geom_smooth(aes(shell_length, preds + 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  geom_smooth(aes(shell_length, preds - 1.96 * se), col = 1, lty = 2, lwd = 0.5, se = F) +
  facet_wrap(~fshape, nrow = 2, scale = "free") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + 
  scale_y_continuous(labels = scaleFUN) +
  xlab("Salinity") + ylab("shape-PCs")

# Estimate effect size --------------------------------------------------
# Estimating effect size as in Telesca et al., 2018

# Function from Telesca et al. 2018 to normalize predictors
my_norm <- function(x) (x-mean(x))/sd(x)

# Normalise variables
df_pca_long_norm <- df_pca_long %>% 
  mutate(Temp = my_norm(Temp),
         ALK = my_norm(ALK),
         SIR = my_norm(SIR),
         shell_length = my_norm(shell_length))

# Fit preferred model with normalised variables as gamm to extract lme component
mody_norm <- gamm(svalue ~ fshape +
              s(Temp, bs="cr", k = 3, by=fshape) + 
              s(ALK, bs="cr", k = 3, by=fshape) + 
              s(SIR, bs="cr", k = 3, by=fshape) +
              s(shell_length, bs="cr", k = 3, by=fshape),
            method = "REML", 
            random = list(flocation = ~1),
            family= gaussian(link="identity"),
            data = df_pca_long_norm, 
            control = gam.control(trace = T))

summary(mody_norm$gam)

# Calculate confidence intervals for effects
ci_mody_norm <- intervals(mody_norm$lme, which = "fixed")

# As data frame & remove intercepts
ci_mody_norm <- data.frame(ci_mody_norm$fixed) %>% 
  slice(4:n()) %>% 
  rownames_to_column(var="names") %>% 
  mutate(names = as.factor(names))

# Add PCs names
ci_mody_norm$PC <- rep(c("shape-PC1", "shape-PC2","shape-PC3"), rep=5)

# Add variable abbreviations
ci_mody_norm$vars <- c(rep("T",3),rep("Alk",3),rep("SIR",3),rep("SH",3))

# Order variables for ggplot
ci_mody_norm$vars<- factor(ci_mody_norm$vars, levels=c("T", "Alk", "SIR", "SH"))

# Plot
ggplot(ci_mody_norm, aes(vars, ymin = lower, ymax = upper, geom = "pointrange")) +
  geom_hline(yintercept = 0, alpha = I(5/12), lty = 2) + 
  geom_errorbar(position = position_dodge(width = 0.3), width = 0.1) + 
  geom_point(aes(vars, est.), position = position_dodge(width = 0.3), size = 2.6, shape = 21, fill = "white") + 
  facet_wrap(~PC, ncol = 3) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Coefficients") + 
  ylab("Standardised regression estimates")
     
# Export model output -------------------------------------------------
# Export environmental data
write_csv(predframe_temp, "data/predframe_temp.csv")
write_csv(predframe_alk, "data/predframe_alk.csv")
write_csv(predframe_SIR, "data/predframe_SIR.csv")
write_csv(predframe_ash, "data/predframe_ash.csv")
write_csv(ci_mody_norm, "data/predframe_eff_size.csv")
saveRDS(mody, "data/shell_shape_model.rds")
saveRDS(mody_norm, "data/norm_shell_shape_model.rds")

                                                                                                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                                                                                                               



