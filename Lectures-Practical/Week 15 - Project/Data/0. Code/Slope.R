rm(list=ls());gc()
library(raster)
library(rgdal)
library(maptools)
data(wrld_simpl)

#------------------------------------------------------------------------------------
# Non-Bioscore predictors based on SRTM

# Slopes: rate of change of elevation coming from the SRTM-30m map <https://www2.jpl.nasa.gov/srtm/>.
#         in Bio----- are defined as areas that, after modeling the digital elevation model
#         (http://download.kortforsyningen.dk/content/dhmterr%C3%A6n-16-m-grid)
#         on a 9.6 x 9.6 m scale, has an average slope of more than 15 degrees.
#         The layer is set to 0 for the places where intensive fields overlap
#         with the slopes. Sloping terrain is included as a potential indicator
#         of biodiversity because sloping areas are more difficult to cultivate
#         than flat terrain, and thus have an increased probability of containing
#         important habitats. --> in Bioscore 1
#------------------------------------------------------------------------------------
# Load the SRTM 30m resolution map
DNKSRTM30m <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Denmark SRTM/DK-SRTM30m.tif")

# Estimate Mean elevation
DNKSRTM1km <- aggregate(DNKSRTM30m,
                        fact = c(30*(1/(60*60)))/res(DNKSRTM30m),
                        fun = mean)

# Estimate Slope
DNKSRTM1kmSlope <- terrain(DNKSRTM1km,opt="slope", unit="degrees", neighbors=8)

## Save the Raster
writeRaster(x = DNKSRTM1kmSlope, # raster to save
            filename = "~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/Slope/Slope.tif", # Name of the raster
            format="GTiff", overwrite=TRUE) # Arguments


# Estimate Slope 30m
DNKSRTM30mSlope <- terrain(DNKSRTM30m,opt="slope", unit="degrees", neighbors=8)

# Aggregate the Slopes
DNKSRTMSlope30mAgg <- aggregate(DNKSRTM30mSlope,
                                fact = c(30*(1/(60*60)))/res(DNKSRTM30m),
                                fun = mean)

## Save the Raster
writeRaster(x = DNKSRTMSlope30mAgg, # raster to save
            filename = "~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/Slope/Slope_30mAgg.tif", # Name of the raster
            format="GTiff", overwrite=TRUE) # Arguments