# -------------------------------------------------------------------- #
#                           Jessica FR Coelho                          #
#                               09/03/2022                             #
#                 LGM modelling with occurrences filtered              #
#                      Suitability maps & metrics                      #
# -------------------------------------------------------------------- #
getwd()
#DELL, 15-08-2022
#R version 4.0.5, 32 bits
#C:/Users/jessf/Documents/UFRN_PhD/LGM_modelling

#Get error messages in English
Sys.setenv(lang = "en_US")

#Load libraries
library(sp)
library(raster)
library(dismo)
library(rgeos)
library(maptools)
library(rJava)
library(rgdal)
library(dplyr)
library(sf)
library(sdmpredictors)
set.seed(123)

# -------------------------------------------------------------------- #
#                        Maps and visualization
# -------------------------------------------------------------------- #
#Visualization purposes and delimitation of model calibration area
#Import Spalding ecoregions, ocean and countries shapes
countries <- readOGR('maps_shapes/countries/10m_admin_0_countries.shp')
ecoregions <- readOGR('maps_shapes/spalding_ecoregions/meow_ecos_expl_clipped_expl.shp')
ocean <- st_read(dsn = "maps_shapes/ne_110m_ocean",
                 layer = "ne_110m_ocean",
                 quiet = TRUE)
#ocean_ogr <- readOGR(dsn = "maps_shapes/ne_110m_ocean", layer = "ne_110m_ocean")

#Visualize map
plot(st_geometry(ocean),
     col = "grey80",
     border = "grey",
     axes = T,
     xlim = c(-70, -25),
     ylim = c(-45, 30),
     lwd = 0.8)

# -------------------------------------------------------------------- #
#             Load abiotic data and explore marine datasets 
# -------------------------------------------------------------------- #
library(sdmpredictors) #see vignette for reference
marine <- list_datasets(terrestrial = FALSE, marine = TRUE)

#Explore layers 
layers <- list_layers(marine)
paleo <- list_layers_paleo(terrestrial = FALSE)
unique(paleo$epoch)
unique(paleo$model_name) 

#Example to get info of different layers
get_layers_info(c("MS_bathy_5m", #current
                  "BO_B1_2100_sstmax", #future
                  "MS_bathy_21kya"))$common #paleo

#(down)load specific layers 
#4 is the number of variables
present4 <- load_layers(c("MS_bathy_5m",
                          "MS_biogeo06_bathy_slope_5m",
                          "MS_biogeo08_sss_mean_5m",
                          "MS_biogeo13_sst_mean_5m"))
plot(present4)

#Check for correlation
test_cor <- getValues(present4)
correlation.matrix <- cor(test_cor, use = "complete.obs")
#write.csv(correlation.matrix, 'correlation_matrix_present.csv')
cor.table = as.data.frame(as.table(correlation.matrix))
high.cor <- subset(cor.table, abs(Freq) > 0.7 & Freq != 1)
high.cor #0

#Get the equivalent paleo layer code for a current climate layer 
get_paleo_layers(c("MS_bathy_5m",
                   "MS_biogeo06_bathy_slope_5m",
                   "MS_biogeo08_sss_mean_5m",
                   "MS_biogeo13_sst_mean_5m"), 
                 model_name = c("21kya_geophysical",
                                "21kya_ensemble_adjCCSM"), 
                 years_ago = 21000)$layer_code

#Then load these layers
paleo4 <- load_layers(c("MS_bathy_21kya",
                        "MS_biogeo06_bathy_slope_21kya",
                        "MS_biogeo08_sss_mean_21kya_adjCCSM",
                        "MS_biogeo13_sst_mean_21kya_adjCCSM"))
plot(paleo4)

#Check correlation
test_cor <- getValues(paleo4)
correlation.matrix <- cor(test_cor, use = "complete.obs")
cor.table = as.data.frame(as.table(correlation.matrix))
high.cor <- subset(cor.table, abs(Freq) > 0.8 & Freq != 1)
high.cor #0

#Get citation to include in paper
print(layer_citations("MS_biogeo13_sst_mean_21kya_adjCCSM"))

#Cutting layers to ecoregions' shape and crop to central/south Americas extent
e <- extent(-90, -25, -45, 40)
study_area <- crop(present4, e) #Crop to the extent
pres_layers <- mask(study_area, mask = ecoregions) #Delimitation to shape geometry
plot(pres_layers)

study_area_lgm <- crop(paleo4, e)
lgm_layers <- mask(study_area_lgm, mask = ecoregions)
plot(lgm_layers)

#Rename layers equally in both climates
names(pres_layers) <- c("bathy",
                        "slope",
                        "sss_mean",
                        "sst_mean")
names(lgm_layers) <- c("bathy",
                       "slope",
                       "sss_mean",
                       "sst_mean")
plot(pres_layers)
plot(lgm_layers)

# Export rasters to run Circuitscape
# Present/current layers
writeRaster(pres_layers[[1]], 
            filename = 'abiotic_data/current_bioregion/bathy.asc',
            overwrite = T)
writeRaster(pres_layers[[2]], 
            filename = 'abiotic_data/current_bioregion/slope.asc',
            overwrite=T)
writeRaster(pres_layers[[3]], 
            filename = 'abiotic_data/current_bioregion/sss_mean.asc',
            overwrite = T)
writeRaster(pres_layers[[4]], 
            filename = 'abiotic_data/current_bioregion/sst_mean.asc',
            overwrite = T)

# LGM layers
writeRaster(lgm_layers[[1]], 
            filename = 'abiotic_data/lgm_bioregion/bathy.asc',
            overwrite = T)
writeRaster(lgm_layers[[2]], 
            filename = 'abiotic_data/lgm_bioregion/slope.asc',
            overwrite=T)
writeRaster(lgm_layers[[3]], 
            filename = 'abiotic_data/lgm_bioregion/sss_mean.asc',
            overwrite = T)
writeRaster(lgm_layers[[4]], 
            filename = 'abiotic_data/lgm_bioregion/sst_mean.asc',
            overwrite = T)

# -------------------------------------------------------------------- #
#                        Load occurrence records 
# -------------------------------------------------------------------- #
occ <- read.csv("Harengula_sp_joint_30.csv", h=T)

#Visualize
plot(pres_layers[[4]])
points(occ[ ,2:3], 
       pch = 21, 
       cex = 0.8,
       col = "red",
       bg = "red")

# Spatial thinning
library(spThin)
occ_subset <- occ_suppl[ ,2:3]
occ_suppl_t <- thin.algorithm(rec.df.orig = occ_subset, 
                              thin.par = 100,
                              reps = 10)
View(occ_suppl_t)
thin_rep1 <- occ_suppl_t[[5]] #get any random set: 28 occ, 100km

plot(pres_layers[[1]]) #Visualize
points(thin_rep1, 
       pch = 21, 
       cex = 0.8,
       col = "red",
       bg = "red")

# -------------------------------------------------------------------- #
#        Simple Modeling: do this first, then alter parameters
# -------------------------------------------------------------------- #
library(dismo)
predictors <- stack(pres_layers)
sp.occ <- occ[ ,2:3]
sp.occ <- thin_rep1

me <- maxent(predictors,
             sp.occ,
             path = 'results_SM')
model = predict(me, predictors)

#Modeling with a few altered parameters
library(dismo)
me2 <- maxent(predictors, 
              sp.occ,
              args = c("-J",
                       "-P",
                       "randomtestpoints=25",
                       "replicates=20",
                       "replicatetype=bootstrap",
                       "randomseed",
                       "removeduplicates",
                       "doclamp=false", 
                       "extrapolate=true",
                       "convergencethreshold=0.00001",
                       "maximumiterations=10000",
                       "maximumbackground=10000"),
              path='results_PAR2')
model.rep <- predict(me2, predictors)
model.mean <- mean(model.rep)
plot(model.mean)

#Save raster to open on QGIS and Circuitscape later
writeRaster(model.mean,
            filename = 'results/Hsp_present.grd',
            overwrite = F)

#LGM projection
lgm <- predict(me2, lgm_layers)
lgm.model <- mean(lgm)
plot(lgm.model)

writeRaster(lgm.model,
            filename = "results_PAR2/Hsp_LGM42",
            format = "ascii",
            NAFlag = "-9999")
plot(lgm.model)

# -------------------------------------------------------------------- #
#                         Visualize models
# -------------------------------------------------------------------- #
par(mfrow = c(1,2))
plot(lgm.model, main = "LGM")        #Plot LGM model
plot(model.mean, main = "Present")   #Plot present model

# -------------------------------------------------------------------- #
#                         Suitability differences
# -------------------------------------------------------------------- #
# Suitability map change
pres_lgm <- model.mean - lgm.mean
plot(pres_lgm)

# -------------------------------------------------------------------- #
#               Extracting suitability metrics - Present
# -------------------------------------------------------------------- #
occ_dart <- read.delim(file = "occ_dartseq.txt",
                       header = TRUE,
                       sep = ',')
View(occ_dart)

# Get model from present and transform grids lower than 0 to NA
values(model.mean)[values(model.mean) < 0] = NA

# Extracting mean suitability values per polygon around sites
# Get EPSG here: http://projfinder.com/
library(sf)
dat_sf <- st_as_sf(occ_dart[ ,3:4],
                   coords = c("Lon", "Lat"),
                   crs = 4055)

# Buffer circles 100.000 m (100 km)
circles <- st_buffer(dat_sf, dist = 100000)
plot(circles)

circles_pos = st_point_on_surface(circles)

# Visualize
par(mfrow = c(1,1))
plot(model.mean,
     main = "Present",
     font = 2)
points(circles,
       pch = 21,
       cex = 1,
       col = "black",
       bg = "white")    
points(occ_dart[ ,3:4], 
       pch = 21, 
       cex = 0.6,
       col = "red",
       bg = "red")

# Extract suitability
suit.pres.mean <- extract(model.mean,
                          circles,
                          fun = mean,
                          na.rm = T,
                          df = T,
                          sp = T)

# -------------------------------------------------------------------- #
#               Extracting suitability metrics - LGM
# -------------------------------------------------------------------- #
# Extract suitability
# Circles buffered 100km
suit.lgm.mean <- extract(lgm.model,
                         circles,
                         fun = mean,
                         na.rm = T,
                         df = T,
                         sp = T)

# -------------------------------------------------------------------- #
#             Suitability differences present minus LGM                #
# -------------------------------------------------------------------- #
pres_lgm2 <- model.mean - lgm.model
plot(pres_lgm) #Present minus LGM model

# Extract suitability of this difference
# Circles buffered 100km
suit.diff.mean <- extract(pres_lgm,
                          circles,
                          fun = mean,
                          na.rm = T,
                          df = T,
                          sp = T)
occ.diff.mean <- cbind(occ_dart, suit.diff.mean)
View(occ.diff.mean)

# --------------------------------------------------------------------------- #
#                            jessicovsky@gmail.com                            #
# --------------------------------------------------------------------------- #