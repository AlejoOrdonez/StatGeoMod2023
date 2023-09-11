rm(list=ls());gc()
library(raster)
library(rgdal)
library(maptools)
data(wrld_simpl)
#------------------------------------------------------------------------------------
# Non-Bioscore predictors based on SRTM

# Bioclimatic variables
## Climate information - BIOCLIM
setwd("/Volumes/Crucial X6/Data/GLOBAL/CLIMATE/CURRENT/CHELSA V1-2-1/Current climatology/BIOCLIM")
for(Var in dir()){#(Var <- dir()[1])
  ##  Load the Data
  RastTmp <- raster(Var)
  ##  Crop and mask to Denmark
  RastTmp2 <- mask(crop(RastTmp, DMK), DMK)
  ## re-project
  RastTmp2 <- projectRaster(RastTmp2,
                            crs = CRS(EU.LAEA))
  ## Align with richness raster
  RastTmp3 <- resample(RastTmp2, DK.Rich, method = "ngb")
  #### Save the Raster
  writeRaster(RastTmp3,
              paste0("~/Dropbox/Courses_&_Conferences/2021/Courses/Statistical and Geospatial Modeling/Lectures-Practical/14-15. Project/Data/Rasters/2. Climate/DK_BIO", gsub(".tif","",strsplit(Var,"_")[[1]][3]), "-CHELSA.tif"),
              format="GTiff", overwrite=TRUE)
}

## Climate information - extend BIOCLIM
#### Potential Evapotranspiration
setwd("/Volumes/Crucial X6/Data/GLOBAL/CLIMATE/CURRENT/CHELSA V1-2-1/Current climatology/PET")
### Get the DATA
DK_PET_CHELSA <- stack(dir())
##  Crop and mask to Denmark
DK_PET_CHELSA2 <- mask(crop(DK_PET_CHELSA, DMK), DMK)
## re-project
DK_PET_CHELSA2 <- projectRaster(DK_PET_CHELSA2,
                                crs = CRS(EU.LAEA))
### summary
DK_PET_CHELSA3 <- sum(DK_PET_CHELSA2)
## Align with richness raster
DK_PET_CHELSA3 <- resample(DK_PET_CHELSA3, DK.Rich, method = "ngb")
#### Save the Raster
writeRaster(DK_PET_CHELSA3,
            "~/Dropbox/Courses_&_Conferences/2021/Courses/Statistical and Geospatial Modeling/Lectures-Practical/14-15. Project/Data/Rasters/2. Climate/DK_PET-CHELSA.tif",
            format="GTiff", overwrite=TRUE)

#### Relative Humidity
setwd("/Volumes/Crucial X6/Data/GLOBAL/CLIMATE/CURRENT/CHELSA V1-2-1/Current climatology/RH")
### Get the DATA
DK_RH_CHELSA <- stack(dir())
##  Crop and mask to Denmark
DK_RH_CHELSA2 <- mask(crop(DK_RH_CHELSA, DMK), DMK)
## re-project
DK_RH_CHELSA2 <- projectRaster(DK_RH_CHELSA2,
                               crs = CRS(EU.LAEA))
### summary
DK_RH_CHELSA3 <- sum(DK_RH_CHELSA2)
## Align with richness raster
DK_RH_CHELSA3 <- resample(DK_RH_CHELSA3, DK.Rich, method = "ngb")

#### Save the Raster
writeRaster(DK_RH_CHELSA3,
            "~/Dropbox/Courses_&_Conferences/2021/Courses/Statistical and Geospatial Modeling/Lectures-Practical/14-15. Project/Data/Rasters/2. Climate/DK_RH-CHELSA.tif",
            format="GTiff", overwrite=TRUE)
