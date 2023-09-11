rm(list=ls());gc()
library(raster)
library(rgdal)
library(maptools)
data(wrld_simpl)

#------------------------------------------------------------------------------------
# Bioscore predictors based on Basemap

## Nature Density:  Nature density is calculated as the proportion of forests and 
##                  protected light-open nature types in a nationwide network of
##                  cells of 1 x 1 km. 
#------------------------------------------------------------------------------------
# Denmark's polygon
## Load high resolution map of Denmark
DKMap <- shapefile("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/BaselineShapes/DNK_admin/GADML_DNK0.shp")

## Re-project the spatial polygon of Denmark to UTM
DK.Map.UTM <- spTransform(x = DKMap,
                          CRSobj = CRS("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs"))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Land Use Land Cover Data - basemap 2011

if("Basemap.2011.tif"%in%dir("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Basemap/")){
  ## Load the Corine data for Denmark
  Basemap.2011 <- raster("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Corine/Basemap.2011.tif")
}

if(!"Basemap.2011.tif"%in%dir("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Basemap/")){
  ## Load the Full Basemap data
  BasemapRaw.2011 <- raster("/Volumes/Seagate Expansion/Data/DENMARK/LAND COVER/BASEMAP/Data/basemap03_2011_2016_2018/lu_agg_2011.tif")
  ## Mask the data to only DK
  Basemap.2011 <- mask(x = BasemapRaw.2011, # basemap raster
                       mask = DK.Map.UTM,# a shapefile to define what to mask out 
                       filename = "/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Basemap/Basemap.2011.tif", # Name of the raster
                       format="GTiff", overwrite=TRUE) # Arguments
  rm(BasemapRaw.2011);gc()
  ## A quick visualization
  #plot(Basemap.2011)
  
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Land Use Land Cover Data - basemap 2016

if("Basemap.2016.tif"%in%dir("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Basemap/")){
  ## Load the Corine data for Denmark
  Basemap.2016 <- raster("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Corine/Basemap.2016.tif")
}


if(!"Basemap.2016.tif"%in%dir("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Basemap/")){
  ## Load the Full Basemap data
  BasemapRaw.2016 <- raster("/Volumes/Seagate Expansion/Data/DENMARK/LAND COVER/BASEMAP/Data/basemap03_2011_2016_2018/lu_agg_2016.tif")
  ## Mask the data to only DK
  Basemap.2016 <-mask(x = BasemapRaw.2016, # Corine Map
                      mask = DK.Map.UTM, # a shapefile to define what to mask out 
                      filename = "/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Basemap/Basemap.2016.tif", # Name of the raster
                      format="GTiff", overwrite=TRUE) # Arguments
  rm(BasemapRaw.2016);gc()
  ## A quick visualization
  #plot(Basemap.2016)
  
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Land Use Land Cover Data - basemap 2018

if("Basemap.2018.tif"%in%dir("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Basemap/")){
  ## Load the Corine data for Denmark
  Basemap.2018 <- raster("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Corine/Basemap.2018.tif")
}

if(!"Basemap.2018.tif"%in%dir("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Basemap/")){
  ## Load the Full Basemap data
  BasemapRaw.2018 <- raster("/Volumes/Seagate Expansion/Data/DENMARK/LAND COVER/BASEMAP/Data/basemap03_2011_2016_2018/lu_agg_2018.tif")
  ## Mask the data to only DK
  Basemap.2018 <-mask(x = BasemapRaw.2018, # Corine Map
                      mask = DK.Map.UTM, # a shapefile to define what to mask out 
                      filename = "/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Basemap/Basemap.2018.tif", # Name of the raster
                      format="GTiff", overwrite=TRUE) # Arguments
  rm(BasemapRaw.2018);gc()
  ## A quick visualization
  #plot(Basemap.2018)
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++