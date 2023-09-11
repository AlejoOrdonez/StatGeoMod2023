rm(list=ls());gc()
library(raster)
library(rgdal)
library(maptools)
data(wrld_simpl)

#------------------------------------------------------------------------------------
# Bioscore predictors based on Basemap

# Wetlands: Points for areas that lie on low-lying land, but not intensively
#                 cultivated fields.
#------------------------------------------------------------------------------------
# Denmark's polygon
## Load high resolution map of Denmark
DKMap <- shapefile("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/BaselineShapes/DNK_admin/GADML_DNK0.shp")

## Re-project the spatial polygon of Denmark to UTM
DK.Map.UTM <- spTransform(x = DKMap,
                          CRSobj = CRS("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs"))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Estimations Based on Basemap
## Load the Full Basemap data
BasemapRaw.2011 <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Basemap/lu_agg_2011.tif")
BasemapRaw.2016 <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Basemap/lu_agg_2016.tif")
BasemapRaw.2018 <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Basemap/lu_agg_2018.tif")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Load legend of basemap
BasemapLeg <- readxl::read_excel("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/Basemap/Aggregated Legend.xlsx")
BasemapLeg <- as.data.frame(BasemapLeg)

### Categories - the First column is a tree number code defining a hierarchy of the levels [see basemap report use guide appendix] <https://dce2.au.dk/pub/TR159.pdf>
# The base levels are 
### 1 Building-Urban
### 2 Agriculture
### 3 Wetlands
### 4 Water

### Wetlands include
### 31	Forest
### 32	Open Nature

### Open Nature include
#### 321 Nature, dry; Agriculture, extensive
#### 322 Nature, wet -->  WETLANDS
#### 323 Nature, wet; Agriculture, extensive -->  WETLANDS


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Nature Density based on Basemap

## Define the value for Wetlands [There have a value of 3 at level 1]
NatAreaVal <- c(322,323) 

## Build a Reclassification matrix
BaseMapUnique <- unique(getValues(BasemapRaw.2011))
ReclassMrtx <- cbind(as.numeric(BaseMapUnique), # These are the "is" values [those in the raster]
                     ifelse(substring(as.vector(BaseMapUnique),1,3)==NatAreaVal,1, # These are the "becomes" values [the NEW values]
                            ifelse(substring(as.vector(BaseMapUnique),1,1)%in%c(4,8,9),NA,0))) # Water/unmap and Non-Dk values are set tp NA
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Estimations Based on Basemap 2011 
a<-Sys.time()
### Reclassify the raster so that it ONLY shows the areas considered FOREST AND SEMI- Wetlands 
Basemap.2011.WetLnd <- reclassify(x = BasemapRaw.2011, # raster to reclasify
                               rcl = ReclassMrtx # the Reclassification matrix
)
### Now lets look at it
plot(Basemap.2011.WetLnd)

### Now we are ready to aggregate the 10x10m reclassified basemap raster into a 1x1km.
### Here you will determine the density of Nature areas per 1km [km of Nature per km

# Aggregate the 100x100m dataset
Summ.Func <- function(x, na.rm=TRUE){sum(x, na.rm=TRUE)/sum(!is.na(x))} # Define the density of points of a class [sum(x)] in all land sites [sum(!is.na(x))]
Basemap.WetLnd1k.2011 <- aggregate(x = Basemap.2011.WetLnd, # raster to aggregate
                                fact = 100, # Aggregation factor (number of cells in each direction to merge together)
                                fun = Summ.Func) # how to aggregate the cells
### Now lets look at this reclassified map
plot(Basemap.WetLnd1k.2011)
Sys.time()-a
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Estimations Based on Basemap 2016 
a<-Sys.time()
### Reclassify the raster so that it ONLY shows the areas considered FOREST AND SEMI- Wetlands 
Basemap.2016.WetLnd <- reclassify(x = BasemapRaw.2016, # raster to reclasify
                               rcl = ReclassMrtx # the Reclassification matrix
)
### Now lets look at it
plot(Basemap.2016.WetLnd)

### Now we are ready to aggregate the 10x10m reclassified basemap raster into a 1x1km.
### Here you will determine the density of Nature areas per 1km [km of Nature per km

# Aggregate the 100x100m dataset
Summ.Func <- function(x, na.rm=TRUE){sum(x, na.rm=TRUE)/sum(!is.na(x))} # Define the density of points of a class [sum(x)] in all land sites [sum(!is.na(x))]
Basemap.WetLnd1k.2016 <- aggregate(x = Basemap.2016.WetLnd, # raster to aggregate
                                fact = 100, # Aggregation factor (number of cells in each direction to merge together)
                                fun = Summ.Func) # how to aggregate the cells
### Now lets look at this reclassified map
plot(Basemap.WetLnd1k.2016)
Sys.time()-a
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Estimations Based on Basemap 2018 
a<-Sys.time()
### Reclassify the raster so that it ONLY shows the areas considered FOREST AND SEMI- Wetlands 
Basemap.2018.WetLnd <- reclassify(x = BasemapRaw.2018, # raster to reclasify
                               rcl = ReclassMrtx # the Reclassification matrix
)
### Now lets look at it
plot(Basemap.2018.WetLnd)

### Now we are ready to aggregate the 10x10m reclassified basemap raster into a 1x1km.
### Here you will determine the density of Nature areas per 1km [km of Nature per km

# Aggregate the 100x100m dataset
Summ.Func <- function(x, na.rm=TRUE){sum(x, na.rm=TRUE)/sum(!is.na(x))} # Define the density of points of a class [sum(x)] in all land sites [sum(!is.na(x))]
Basemap.WetLnd1k.2018 <- aggregate(x = Basemap.2018.WetLnd, # raster to aggregate
                                fact = 100, # Aggregation factor (number of cells in each direction to merge together)
                                fun = Summ.Func) # how to aggregate the cells
### Now lets look at this reclassified map
plot(Basemap.WetLnd1k.2018)
Sys.time()-a

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Average of Wetlands density in 2012 and 2018
Basemap.WetLnd1k <- mean(stack(Basemap.WetLnd1k.2018,Basemap.WetLnd1k.2016,Basemap.WetLnd1k.2011))

plot(Basemap.WetLnd1k)

## Save the Raster
writeRaster(x = Basemap.WetLnd1k, # raster to save
            filename = "~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/Wetland Density/WetLndDens_BaseMap.tif", # Name of the raster
            format="GTiff", overwrite=TRUE) # Arguments
removeTmpFiles(h=1/60)
gc()

