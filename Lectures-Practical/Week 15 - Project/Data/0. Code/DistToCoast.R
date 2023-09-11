rm(list=ls());gc()
library(raster)
library(rgdal)
library(maptools)
data(wrld_simpl)

#------------------------------------------------------------------------------------
# Bioscore predictors based on Base cartography

# Distance to the coast: Areas that are less than 1 km from the coast. As for the
#                        other proxies, no points are given for intensively
#                       cultivated fields.

#------------------------------------------------------------------------------------
# Denmark's polygon
## Load high resolution map of Denmark
DKMap <- shapefile("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/BaselineShapes/DNK_admin/GADML_DNK0.shp")

## Re-project the spatial polygon of Denmark to UTM
DK.Map.UTM <- spTransform(x = DKMap,
                          CRSobj = CRS("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs"))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#### Create a grided map of Vascular plant richness using the raster package
DMK.raster <- raster(DMI.DK.Shp.Proj) # Build a raster to summarise observations and observations
# Make a 1km Resolution raster
res(DMK.raster) <- 10000
DMK.raster[] <- 1:ncell(DMK.raster)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Transform Denmark from a polygon to a line
bound <- as(DK.Map.UTM, 'SpatialLinesDataFrame')
# Get the points for each truning point of the shore line
bound2 <- do.call("rbind",coordinates(bound)[[1]])

# Get the UTM lat long for all 1km cells in DK
DK.In <- coordinates(DMK.raster)[!is.na(DMK.raster[]),]

# Estimate the min distance for each cell to the shoreline
DK.dist <- apply(DK.In, 1, function(x){
  min(spDistsN1(bound2,x))
})
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create a Summary raster
DK.DistRast <- raster(DMK.raster)
# attach the Distances to the summary raster
DK.DistRast[which(!is.na(DMK.raster[]))] <- DK.dist

#### Save the Raster
writeRaster(DK.DistRast,
            "~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/Distance to Coast/DK_DISTtoCOAST.tif",
            format="GTiff", overwrite=TRUE)

