rm(list=ls());gc()
library(raster)
library(rgdal)
library(maptools)
data(wrld_simpl)

#------------------------------------------------------------------------------------
# Non-Bioscore predictors based on SRTM

# Elevation Heterogeneity: Defined as the variability in elevations coming from
#                          the SRTM-30m map <https://www2.jpl.nasa.gov/srtm/>.
#------------------------------------------------------------------------------------
# Load the SRTM 90m resolution map
DNKSRTM30m <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Denmark SRTM/DK-SRTM30m.tif")

# Estimate elevation variability
DK_ELEVRng_SRTM <- aggregate(DNKSRTM30m,
                             fact = c(30*(1/(60*60)))/res(DNKSRTM30m),
                             fun = sd,
                             na.rm=T)
#------------------------------------------------------------------------------------
## Save the Raster
writeRaster(x = DK_ELEVRng_SRTM, # raster to save
            filename = "~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/Elevantion Heterogenity/SRTM30_SD_1km.tif", # Name of the raster
            format="GTiff", overwrite=TRUE) # Arguments

plot(DK_ELEVRng_SRTM)
removeTmpFiles(h=1/60)
gc()


# Build DK SRTM30 map
# rm(list=ls());gc()
# library(raster)
# library(rgdal)
# library(maptools)
# data(wrld_simpl)
# 
# #Load SRTM raster
# SRTM30m <- raster("~/Downloads/output_SRTMGL1.tif")
# plot(SRTM30m)
# # Denmark's polygon
# ## Load high resolution map of Denmark
# DKMap <- shapefile("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/BaselineShapes/DNK_admin/GADML_DNK0.shp")
# # crop to DKMap
# SRTM30mDKreg <- crop(SRTM30m, DKMap)
# #mask to DKMap
# SRTM30mDKMask <- mask(SRTM30mDKreg, DKMap)
# # reportject to UTM
# SRTM30mDKMask2 <- projectRaster(SRTM30mDKMask, crs= CRS("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs"), method="ngb")
# # Save
# writeRaster(SRTM30mDKMask2,"/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Denmark SRTM/DK-SRTM30m.tif", format="GTiff", overwrite=TRUE)
