rm(list=ls());gc()
library(raster)
library(rgdal)
library(maptools)
data(wrld_simpl)

#------------------------------------------------------------------------------------
# Bioscore predictors based on Corine

## Coniferous Density:  Coniferous Density is calculated as the proportion 
##                      of a specific forests type [Broad-leaved forest, Coniferous
##                      forest, Mixed forest] in a nationwide network of cells
##                      of 1 x 1 km.
#------------------------------------------------------------------------------------
# Denmark's polygon
## Load high resolution map of Denmark
DKMap <- shapefile("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/BaselineShapes/DNK_admin/GADML_DNK0.shp")

## Re-project the spatial polygon of Denmark to laea
DK.Map.laea <- spTransform(x = DKMap,
                           CRSobj = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Land Use Land Cover Data - Corine 2012

if("CorineDK.2012.tif"%in%dir("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Corine/")){
  ## Load the Corine data for Denmark
  CorineDK.2012 <- raster("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Corine/CorineDK.2012.tif")
}

if(!"CorineDK.2012.tif"%in%dir("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Corine/")){
  ## Load the Full Corine data
  CorineDK.2012 <- raster("/Volumes/Seagate Expansion/Data/EUROPE/LAND COVER/PRESENT/Corine Land Cover/v2020/2012/DATA/U2018_CLC2012_V2020_20u1.tif")
  
  ## Crop the data to only DK
  CorineDK.Crop <-crop(x = CorineEU.2012, # Corine Map
                       y = extent(DK.Map.laea) # The extent of a shapefile to define the crop 
  )
  
  ## Mask the data to only DK
  CorineDK.2012 <-mask(x = CorineDK.Crop, # Corine Map
                       mask = DK.Map.laea # a shapefile to define what to mask out 
  )
  rm(CorineDK.Crop)
  ## A quick visualization
  plot(CorineDK.2012)
  writeRaster(x = CorineDK.2012, # raster to save
              filename = "/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Corine/CorineDK.2012.tif", # Name of the raster
              format="GTiff", overwrite=TRUE) # Arguments
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Land Use Land Cover Data - Corine 2018
if("CorineDK.2018.tif"%in%dir("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Corine/")){
  ## Load the Corine data for Denmark
  CorineDK.2018 <- raster("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Corine/CorineDK.2018.tif")
}

if(!"CorineDK.2018.tif"%in%dir("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Corine/")){
  ## Load the Full Corine data
  CorineEU.2018 <- raster("/Volumes/Seagate Expansion/Data/EUROPE/LAND COVER/PRESENT/Corine Land Cover/v2020/2018/DATA/U2018_CLC2018_V2020_20u1.tif")
  
  ## Crop the data to only DK
  CorineDK.Crop <-crop(x = CorineEU.2018, # Corine Map
                       y = extent(DK.Map.laea) # The extent of a shapefile to define the crop 
  )
  
  ## Mask the data to only DK
  CorineDK.2018 <-mask(x = CorineDK.Crop, # Corine Map
                       mask = DK.Map.laea # a shapefile to define what to mask out 
  )
  rm(CorineDK.Crop)
  ## A quick visualization
  plot(CorineDK.2018)
  writeRaster(x = CorineDK.2018, # raster to save
              filename = "/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Corine/CorineDK.2018.tif", # Name of the raster
              format="GTiff", overwrite=TRUE) # Arguments
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Load Corine legend file
Legend <- read.csv("/Volumes/Seagate Expansion/Data/EUROPE/LAND COVER/PRESENT/Corine Land Cover/v2020/2012/Legend/CLC2018_CLC2012_V2018_20_QGIS.txt",header=F)
head(Legend) 
### Categories - the First column is a tree number code defining a hierarchy of the levels [see CLC use guide appendix] <https://land.copernicus.eu/user-corner/technical-library/clc-product-user-manual>
### the first aggregation
### 1. ARTIFICIAL SURFACES
### 2. AGRICULTURAL AREAS
### 3. FOREST AND SEMI- NATURAL AREAS
### 4. WETLANDS
### 5. WATER BODIES
### 9. NA

### the Second aggregation for FOREST AND SEMI- NATURAL AREAS is:
### 3.1. Forests
### 3.2. Scrub and/or herbaceous associations
### 3.3. Scrub and/or herbaceous associations

### the Third aggregation for FOREST is:
### 3.1.1. Broad-leaved forest
### 3.1.2. Coniferous forest
### 3.1.3. Mixed forest

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Coniferous forest Density based on Corine

## Define the value for Coniferous forest areas 
NatAreaVal <- c(312) 

## Build a Reclassification matrix
ReclassMrtx <- cbind(1:dim(Legend)[1],  # These are the "is" values [those in the raster]
                     Legend[,1]%in%NatAreaVal # These are the "becomes" values [the NEW values]
)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Estimations Based on Corine 2012 
### Reclassify the raster so that it ONLY shows the areas considered FOREST AND SEMI- NATURAL AREAS 
CorineDK.2012.CF.Forest <- reclassify(x = CorineDK.2012, # raster to reclasify
                                      rcl = ReclassMrtx # the Reclassification matrix
)
### Now lets look at it
plot(CorineDK.2012.CF.Forest)

### Now we are ready to aggregate the 100x100m base Corine raster into a 1x1km.
### Here you will determine the density of Nature areas per 1km [km of Nature per km]

### Aggregate the 100x100m dataset
Summ.Func <- function(x, na.rm=TRUE){sum(x, na.rm=TRUE)/sum(!is.na(x))} # Define the density of points of a class [sum(x)] in all land sites [sum(!is.na(x))]

CorineDK.2012.CF.Forest1k <- aggregate(x = CorineDK.2012.CF.Forest, # raster to aggregate
                                       fact = 10, # Aggregation factor (number of cells in each direction to merge together)
                                       fun = Summ.Func) # how to aggregate the cells
### Now lets look at this reclasified map
plot(CorineDK.2012.CF.Forest1k)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Estimations Based on Corine 2018 
### Reclassify the raster so that it ONLY shows the areas considered FOREST AND SEMI- NATURAL AREAS 
CorineDK.2018.CF.Forest <- reclassify(x = CorineDK.2018, # raster to reclasify
                                      rcl = ReclassMrtx # the Reclassification matrix
)
### Now lets look at it
plot(CorineDK.2018.CF.Forest)

### Now we are ready to aggregate the 100x100m base Corine raster into a 1x1km.
### Here you will determine the density of Nature areas per 1km [km of Nature per km]

### Aggregate the 100x100m dataset
Summ.Func <- function(x, na.rm=TRUE){sum(x, na.rm=TRUE)/sum(!is.na(x))} # Define the density of points of a class [sum(x)] in all land sites [sum(!is.na(x))]

CorineDK.2018.CF.Forest1k <- aggregate(x = CorineDK.2018.CF.Forest, # raster to aggregate
                                       fact = 10, # Aggregation factor (number of cells in each direction to merge together)
                                       fun = Summ.Func) # how to aggregate the cells
### Now lets look at this reclasified map
plot(CorineDK.2018.CF.Forest1k)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Average of Open Natural density in 2012 and 2018
CorineDK.CF.Forest1k <- mean(stack(CorineDK.2018.CF.Forest1k,CorineDK.2012.CF.Forest1k))

plot(CorineDK.CF.Forest1k)

## Save the Raster
writeRaster(x = CorineDK.CF.Forest1k, # raster to save
            filename = "/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/Forest Type Density/CFForestDens_Cor.tif", # Name of the raster
            format="GTiff", overwrite=TRUE) # Arguments
removeTmpFiles(h=1/60)
gc()
