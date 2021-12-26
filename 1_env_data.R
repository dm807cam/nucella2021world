library("tidyverse")
library("rnaturalearth")
#devtools::install_github("ropensci/rnaturalearthhires")
library("cmocean")
library("sf")
library("raster")
library("seacarb")
library("latex2exp")

# Import helper file. This file contain proprietary code and
# will not be supplied with the rest of the code.
# All functions used from this file will be highlighted in the script.
source("imp_func.R", local = TRUE)

# Read netCDF and import data --------------------------------------------------------
# netCDF import function
get_netCDF <- function(path_to_file, variable) {
  require("ncdf4")
  require("raster")
  
  quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }
  
  br <-
    quiet(brick(path_to_file,
                varname = variable))
  
  nc <-
    nc_open(path_to_file)
  depth <- ncvar_get(nc, "Depth")
  nc_close(nc)
  
  br <- setZ(br, depth, name = 'water_depth')
  return(br)
}

# Import oxygen netCDF data
brenv_oxygen <-
  get_netCDF(
    "data/env/GLODAP_REP_GRIDDED_FIELDS/CLIMATOLOGY/GLODAPv2.2016b.oxygen_CMEMS.nc",
    "oxygen"
  )

# Import temperature netCDF data
brenv_temperature <-
  get_netCDF(
    "data/env/GLODAP_REP_GRIDDED_FIELDS/CLIMATOLOGY/GLODAPv2.2016b.temperature_CMEMS.nc",
    "temperature"
  )

# Import pH netCDF data
brenv_pH <-
  get_netCDF(
    "data/env/GLODAP_REP_GRIDDED_FIELDS/CLIMATOLOGY/GLODAPv2.2016b.pHtsinsitutp_CMEMS.nc",
    "pHtsinsitutp"
  )

# Import salinity netCDF data
brenv_salinity <-
  get_netCDF(
    "data/env/GLODAP_REP_GRIDDED_FIELDS/CLIMATOLOGY/GLODAPv2.2016b.salinity_CMEMS.nc",
    "salinity"
  )

# Import alkalinity netCDF data
brenv_alklinity <-
  get_netCDF(
    "data/env/GLODAP_REP_GRIDDED_FIELDS/CLIMATOLOGY/GLODAPv2.2016b.Talk_CMEMS.nc",
    "TAlk"
  )

# Import CO2 netCDF data
brenv_CO2 <-
  get_netCDF(
    "data/env/GLODAP_REP_GRIDDED_FIELDS/CLIMATOLOGY/GLODAPv2.2016b.TCO2_CMEMS.nc",
    "TCO2"
  )

# Import world map as sf
world <- ne_countries(scale = 10, returnclass = "sf")

# Convert raster brick to data frame ---------------------------------------------
# Get all raster bricks (all bricks start with "brenv_")
dat <- ls(pattern = "brenv_")

# Extend coordinate system and 
# export raster bricks to individual data frames (var, lon, lat)
for (ii in 1:length(dat)) {
  kachy <- get(dat[ii])$X1
  
  # Extend the coordinates
  e1 = extent(c(
    xmin = 20,
    xmax = 360,
    ymin = -90,
    ymax = 90
  ))
  e2 = extent(c(
    xmin = 360,
    xmax = 360 + 20,
    ymin = -90,
    ymax = 90
  ))
  
  d1 = crop(kachy, e1)
  d2 = crop(kachy, e2)
  
  xmin(d2) = xmin(d2) - 360
  xmax(d2) = xmax(d2) - 360
  
  kachy = rotate(merge(d1, d2))
  
  varname <- c(str_split_fixed(dat[ii], "_", 2)[2], "lon", "lat")
  
  kachy <- as(kachy, "SpatialPixelsDataFrame")
  assign(paste("dfenv", varname[1], sep = "_"),
         data.frame(kachy) %>%
           setNames(varname))
}

# Merge all data frames to single data frame
dat <- ls(pattern = "dfenv_")
dfenv_all <- get(dat[1]) %>%
  full_join(get(dat[2])) %>%
  full_join(get(dat[3])) %>%
  full_join(get(dat[4])) %>%
  full_join(get(dat[5])) %>%
  full_join(get(dat[6]))

# Calculate carbonate system variables using "seacarb" ------------------------------

# Check range
round(range(dfenv_all$salinity),2)
round(range(dfenv_all$temperature),2)

# Remove temperature values (<=0) to work with seacarb "preferred" equations.
dfenv_all <-dfenv_all %>% 
  filter(temperature > 0)

# Calculate the carbonate variables
df_carb <- carb(
  8,
  dfenv_all$pH,
  dfenv_all$alklinity / 10 ^ 6,
  S = dfenv_all$salinity,
  T = dfenv_all$temperature,
  Patm = 1,
  P = 0,
  Pt = 0,
  Sit = 0,
  pHscale = "T",
  eos = "teos10",
  long = dfenv_all$lon,
  lat = dfenv_all$lat,
  warn  = "y"
)

# Calculate SIR
df_carb$SIR <- df_carb$HCO3 / ((10^-df_carb$pH)*1000000)

# Seacarb does not contain coordinates in data frame
# Re-include coordinates
df_carb <- cbind(df_carb,
                 lon = dfenv_all$lon,
                 lat = dfenv_all$lat)

# To make the next step computationally more efficient 
# reduce data frame to sampling range
df_carb <-
  filter(df_carb, lon >= -35 & 
                  lon <= 24 & 
                  lat >= 36 & 
                  lat <= 72)

# Plot environmental data maps --------------------------------------------
# Select variables to plot
dat <- c("S", "T", "pH", "SIR", "OmegaAragonite", "OmegaCalcite")

# Define color scheme
col <- c("haline", "thermal", "speed", "turbid", "dense", "dense")

# Define variable names as LATEX expression
datname <-
  c("S (psu)",
    "T $(^oC)$",
    "pH",
    "SIR",
    "$\\Omega_{Ar}$",
    "$\\Omega_{Ca}$")

# Plot all maps into large panel
datalist <- list()
for (ii in 1:length(dat)) {
  ij <- unique(dat)[ii]
  
  tmp <- ggplot(data = world) +
    geom_raster(data = df_carb,
                aes_string("lon", "lat", fill = ij)) +
    geom_sf(fill = "white", size = 0.2) +
    scale_fill_cmocean(name = col[ii]) +
    coord_sf(xlim = c(-35, 24),
             ylim = c(36, 72),
             expand = F) +
    labs(fill = unname(TeX(datname[ii]))) +
    xlab("longitude (decimal degree)") +
    ylab("latitude (decimal degree)")
  
  datalist[[ij]] <-
    tmp + theme_set(theme_bw(base_size = 10, base_family = 'Times')) +
    theme(legend.position = "right",
          axis.title = element_blank(),
          legend.key.width = unit(0.4,"line"),
          legend.title = element_text(size=8),
          legend.text = element_text(size=7)) 
}
nCol <- floor(sqrt(length(datalist)))

pdf("plots/all_maps.pdf",
    onefile = F)

do.call(cowplot::plot_grid, c(datalist, ncol = nCol, align = "hv"))

dev.off()

# Select nearest env parameters to sampling sites ---------------------------
# Import caliper data
df_caliper <- read_csv("data/shell_dimension_caliper.csv")

# Select and number unique sampling locations
df_env <- df_caliper %>%
  group_by(location) %>%
  summarise(lon = lon[1],
            lat = lat[1]) %>% 
  arrange(desc(lat)) %>% 
  mutate(nlocation = 1:nrow(.))

# Find closest env points to sampling locations
for (ii in 1:nrow(df_env)) {
  kachy <- df_carb
  kachy$dist <- pointDistance(df_env[ii,c(-1,-4)],
                              df_carb[, c(21, 22)],
                              lonlat = T,
                              allpairs = FALSE)
  
  tmp <- kachy %>%
    slice(which.min(kachy$dist)) %>% 
    dplyr::select(-lon, -lat)
  
  datalist[[ii]] <- cbind(tmp, df_env[ii, ])
}
df_env <- do.call(rbind, datalist)

# Rename temperature and salinity variable
df_env <- df_env %>%
  rename(Temp = T,
         Psal = S)

# Check if points were selected successfully
ggplot(data = world) +
  geom_point(data = df_env,
              aes(lon, lat, colour = Psal),
              show.legend = F,
             size = 3) +
  geom_sf(fill = "white") +
  scale_colour_cmocean(name = "thermal") +
  geom_point(
    df_env,
    mapping = aes(lon, lat),
    colour = "black",
    size = 1
  ) +
  coord_sf(xlim = c(-35, 24),
           ylim = c(36, 72),
           expand = F) 
# All looking good!! :)

# Order by lat
df_env <- df_env[order(df_env$lat),]

# Export environmental data
write_csv(df_env, "data/environmental_comp.csv")
write_csv(df_carb, "data/environmental_all.csv")
