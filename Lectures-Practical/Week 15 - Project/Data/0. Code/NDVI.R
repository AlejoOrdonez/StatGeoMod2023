rm(list=ls());gc()
library(raster)
library(rgdal)
library(maptools)
data(wrld_simpl)

rm(list=ls());gc()
library(raster)
library(rgdal)
library(maptools)
data(wrld_simpl)

#------------------------------------------------------------------------------------
# Non-Bioscore predictors based on SRTM

# NDVI: The Normalized Difference Vegetation Index (NDVI) is an indicator of the greenness of the biomes.
#------------------------------------------------------------------------------------
# Denmark's polygon
## Load high resolution map of Denmark
DKMap <- shapefile("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/BaselineShapes/DNK_admin/GADML_DNK0.shp")

## Re-project the spatial polygon of Denmark to UTM
DK.Map.UTM <- spTransform(x = DKMap,
                          CRSobj = CRS("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs"))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load the NDVI Data
Dir.WData <- "/Volumes/Seagate Expansion/Data/DENMARK/NDVI/Data/SHORT TERM STATISTICS/"

NDVI.DirList <- dir(Dir.WData)

NDVI.List <- lapply(NDVI.DirList, 
                    function(x){
                      NDVITmp <- raster(dir(paste0(Dir.WData,x,"/20150101"),pattern="mean",full.names =T))
                      NDVITmp2 <-mask(NDVITmp, DKMap)
                      return(NDVITmp2)
                    })

NDVI.Stack <- do.call("stack",NDVI.List) 

NDVI.Mean <- calc(NDVI.Stack,mean,na.rm=T)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Save the Raster
writeRaster(x = NDVI.Mean, # raster to save
            filename = "~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/NDVI/NDVI.tif", # Name of the raster
            format="GTiff", overwrite=TRUE) # Arguments
removeTmpFiles(h=1/60)
