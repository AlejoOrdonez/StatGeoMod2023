rm(list=ls());gc()
library(raster)
library(rgdal)
library(maptools)
data(wrld_simpl)

#------------------------------------------------------------------------------------
# Non-Bioscore predictors based on SRTM

# Soil heterogeneity: 

#------------------------------------------------------------------------------------
# Load Data
SOILGRID <- raster("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2020/Courses/Statistical & Geospatial modeling/Lectures-Practical/16. Project/Data/Rasters/Soils/SoilsGrid/wrb_0-5cm_mean.tif")
## re-project
SOILGRID2 <- projectRaster(SOILGRID,
                           crs = CRS(EU.LAEA),
                           method= "ngb")
## estimate the number of Soil Types
DK_SOIL_SOILGRD <- aggregate(SOILGRID2,
                             fact = 1000/res(SOILGRID2),
                             fun = function(x , ...){a <- length(unique(na.omit(x)))
                             a <- ifelse(a==0,NA,a)
                             return(a)})
