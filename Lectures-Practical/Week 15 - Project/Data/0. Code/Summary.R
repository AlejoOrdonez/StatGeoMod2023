rm(list=ls());gc()
library(raster)
library(rgdal)
library(maptools)
data(wrld_simpl)

#------------------------------------------------------------------------------------
# Add And values to Biodiveristy Data
# #### NOVANA
# DK.Biodiv.In <- readRDS("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Biodiveristy/NOVANA/NOVANA_NA.Rda")
# sort(names(DK.Biodiv.In))
# DK.Biodiv <- data.frame(ID = 1:dim(DK.Biodiv.In)[1],
#                         UTMLongitude = DK.Biodiv.In$x,
#                         UTMLatitude = DK.Biodiv.In$y,
#                         S.AllGrp = DK.Biodiv.In$antalarter # All evaluated Groups
#                         )
# 
# DK.Biodiv <- na.omit(DK.Biodiv)
# ## Re-project the spatial polygon of Denmark to UTM
# DK.BiodivPntSHP.UTM <- SpatialPointsDataFrame(coords = DK.Biodiv[,c("UTMLongitude","UTMLatitude")],
#                                               proj4string = CRS("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_def"),
#                                               data = DK.Biodiv)
# ## Make a spatial points dataframe
# DK.BiodivPntSHP <- spTransform(x = DK.BiodivPntSHP.UTM,
#                                CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# ## Re-porject the spatial points data.frame to UTM zone of Denmark https://epsg.io/25832 https://sdfe.dk/media/2917583/001-etrs89-utm.pdf
# DK.BiodivPntSHP.laea <- spTransform(x = DK.BiodivPntSHP.UTM,
#                                     CRSobj = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))

#------------------------------------------------------------------------------------
#### BIOWIDE
# load the Raw data
DK.Biodiv.In <- readxl::read_excel("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Biodiveristy/Biowide/occurrence.xlsx")
names(DK.Biodiv.In)

# Which groups where sampled?
unique(DK.Biodiv.In$phylum)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Summarize the richness
DK.BiodivLatLongPaste <- c(paste0(DK.Biodiv.In$decimalLongitude,"-",DK.Biodiv.In$decimalLatitude))
funcCont <- function(x){length(unique(x))}
RichSumm <- tapply(DK.Biodiv.In$scientificName,
                   INDEX = list(DK.BiodivLatLongPaste,
                                DK.Biodiv.In$phylum),
                   FUN = funcCont)
RichSumm <- as.data.frame(RichSumm)
RichSumm$AllGrp <- as.numeric(tapply(DK.Biodiv.In$scientificName,
                                        INDEX = list(DK.BiodivLatLongPaste),
                                        FUN = funcCont))

head(RichSumm)
# Make the first summary table focusing on Arthropoda and Plants
DK.Biodiv <- data.frame(ID = 1:130,
                        decimalLongitude = do.call("c",lapply(strsplit(row.names(RichSumm),"-"),function(x){as.numeric(x[1])})),
                        decimalLatitude = do.call("c",lapply(strsplit(row.names(RichSumm),"-"),function(x){as.numeric(x[2])})),
                        S.AllGrp = RichSumm[,"AllGrp"], # All evaluated Groups
                        S.Arthropoda = RichSumm[,"Arthropoda"], # Species richness Arthropoda
                        S.Spermatophyte = c(apply(RichSumm[,c("Magnoliophyta","Pinophyta")],1,sum, na.rm=T)), # Species richness vascular plant
                        S.Magnoliophyta = RichSumm[,"Magnoliophyta"] # # Species richness flowering plants
)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Prep to extract env data
## Make a spatial points dataframe
DK.BiodivPntSHP <- SpatialPointsDataFrame(coords = DK.Biodiv[,c("decimalLongitude","decimalLatitude")],
                                          proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),
                                          data = DK.Biodiv)
## Re-porject the spatial points data.frame to UTM zone of Denmark https://epsg.io/25832 https://sdfe.dk/media/2917583/001-etrs89-utm.pdf
DK.BiodivPntSHP.laea <- spTransform(x = DK.BiodivPntSHP,
                                    CRSobj = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))

## Re-project the spatial polygon of Denmark to UTM
DK.BiodivPntSHP.UTM <- spTransform(x = DK.BiodivPntSHP,
                                    CRSobj = CRS("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_def"))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Extract the Region
landsdeleSHP <- shapefile("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/BaselineShapes/DNK_admin/GADML_DNK1.shp")
landsdele <- over(DK.BiodivPntSHP, # The SpatialPointDataFrame
                            landsdeleSHP["NAME_1"] # The attribute to extract 
)
DK.Biodiv$landsdele <- as.vector(landsdele[,1])

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Extract the parental material
Geology <- shapefile("/Volumes/Seagate Expansion/Data/DENMARK/Geology/GEUS Surface Geology Map of Denmark/Data/Jordart_200000_Shape/jordart_200000.shp")
Geology.UTM <- spTransform(x = Geology,
                                   CRSobj = CRS("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_def"))
GeologyExt <- over(DK.BiodivPntSHP.UTM, # The SpatialPointDataFrame
                          Geology.UTM["TSYM"] # The attribute to extract
                          )
Legend <- rbind(c("ES","Flyvesand","Aeolian sand"),
                c("F","Ferskvandsdannelser","Freshwater deposits"),
                c("HV","Marsk","Marsh"),
                c("HSL","Marint sand og ler","Marine sand and clay"),
                c("HG","Strandvolde","Beach ridges"),
                c("MSG","Morænesand og grus","Sandy and gravelly Till"),
                c("ML","Moræneler","Clayey Till"),
                c("DSG","Smeltevandssand og -grus","Glaciofluvial"),
                c("DL","Smeltevandsler","Glaciolacustrine"),
                c("T","Extramarginale aflejringer","Downwash sandy deposits"),
                c("Y","Senglaciale havaflejringer","Late-Glacial marine deposits"),
                c("AF","Prækvartær","Pre-Quaternary"),
                c("SØ","Søer","Lake"),
                c("FYLD","Fyld, havne, diger m.m.","Antropogenic constructions"),
                c("SK","Skrivekridt med flint","Chalk"),
                c("ZK","Danienkalk med flint","Lime"))
names(Legend)<-c("code","Navn_DK", "Name_EN")

DK.Biodiv$Geology <- factor(as.vector(GeologyExt[,1]),Legend[,1],Legend[,3])


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Nature Density:  Nature density is calculated as the proportion of forests and 
##                  protected light-open nature types in a nationwide network of
##                  cells of 1 x 1 km. 

# NatDenCor <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/1. Nature Density/NatDens_Cor.tif")
# DK.Biodiv$NatDensCor <- extract(NatDenCor,DK.BiodivPntSHP.laea)
# plot(DK.Biodiv$S.AllGrp~DK.Biodiv$NatDensCor)

NatDenBasemap <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/1. Nature Density/NatDens_BaseMap.tif")
DK.Biodiv$NatDensBasMp <- extract(NatDenBasemap,DK.BiodivPntSHP.UTM)
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$NatDensBasMp)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##  Open nature:   Open Nature density is calculated as the proportion of
##                  protected light-open nature types in a nationwide network of
##                  cells of 1 x 1 km.

# OpnNatCor <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/2. Open Nature Density/OpnNatDens_Cor.tif")
# DK.Biodiv$OpnNatCor <- extract(OpnNatCor,DK.BiodivPntSHP.laea)
# plot(DK.Biodiv$S.AllGrp~DK.Biodiv$OpnNatCor)

OpnNatBasemap <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/2. Open Nature Density/OpnNatDens_Basemap.tif")
DK.Biodiv$OpnNatBasMp <- extract(OpnNatBasemap,DK.BiodivPntSHP.UTM)
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$OpnNatBasMp)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Forest Density:  Forest density is calculated as the proportion of forests in a
##                  nationwide network of cells of 1 x 1 km. 

# ForDenCor <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/3. Forest Density/ForestDens_Cor.tif")
# DK.Biodiv$ForDenCor <- extract(ForDenCor,DK.BiodivPntSHP.laea)
# plot(DK.Biodiv$S.AllGrp~DK.Biodiv$ForDenCor)

ForDenBasemap <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/3. Forest Density/ForestDens_Basemap.tif")
DK.Biodiv$ForDenBasMp <- extract(ForDenBasemap,DK.BiodivPntSHP.UTM)
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$ForDenBasMp)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Agricultural Density:  Agricultural density is calculated as the proportion of 
##                       areas classified as Agricultural in a nationwide network
##                       of cells of 1 x 1 km.

# AgrCor <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/4.Agro Density/AgroDens_Cor.tif")
# DK.Biodiv$AgrCor <- extract(AgrCor,DK.BiodivPntSHP.laea)
# plot(DK.Biodiv$S.AllGrp~DK.Biodiv$AgrCor)

AgrBasemap <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/4.Agro Density/AgroDens_Basemap.tif")
DK.Biodiv$AgrBasMp <- extract(AgrBasemap,DK.BiodivPntSHP.UTM)
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$AgrBasMp)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Wetlands: Points for areas that lie on low-lying land, but not intensively
#                 cultivated fields.

# WetLndCor <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/5. Wetland Density/WetLndDens_Cor.tif")
# DK.Biodiv$WetLndCor <- extract(AgrCor,DK.BiodivPntSHP.laea)
# plot(DK.Biodiv$S.AllGrp~DK.Biodiv$WetLndCor)

WetLndBasemap <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/5. Wetland Density/WetLndDens_BaseMap.tif")
DK.Biodiv$WetLndBasMp <- extract(WetLndBasemap,DK.BiodivPntSHP.UTM)
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$WetLndBasMp)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Distance to the coast: Areas that are less than 1 km from the coast. As for the
#                        other proxies, no points are given for intensively
#                       cultivated fields.
Dis2Cost <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/6. Distance to Coast/DISTtoCOAST.tif")
DK.Biodiv$Dis2Cost <- extract(Dis2Cost,DK.BiodivPntSHP.laea)
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$Dis2Cost)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Human Footprint: The density of man-made lines in the landscape such as roads,
#               ditches and field boundaries. The proxy gives points for areas
#               with less than 8 km of lines per 500 x 500 m.
HII <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/7. HumanFootprint/HII.tif")
DK.Biodiv$HII <- extract(HII,DK.BiodivPntSHP.laea)
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$HII)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Slopes: rate of change of elevation coming from the SRTM-30m map <https://www2.jpl.nasa.gov/srtm/>.
#         in Bio----- are defined as areas that, after modelling the digital elevation model

Slope30mAgg <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/8. Slope/Slope_30mAgg.tif")
DK.Biodiv$Slope30mAgg <- extract(Slope30mAgg,DK.BiodivPntSHP)
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$Slope30mAgg)

# Slope <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/8. Slope/Slope.tif")
# DK.Biodiv$Slope <- extract(Slope,DK.BiodivPntSHP)
# plot(DK.Biodiv$S.AllGrp~DK.Biodiv$Slope)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Elevation Heterogeneity: Defined as the variability in elevations coming from
#                          the SRTM-30m map <https://www2.jpl.nasa.gov/srtm/>.

ElevHet <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/9. Elevantion Heterogenity/SRTM30_SD_1km.tif")
DK.Biodiv$ElevHet <- extract(ElevHet,DK.BiodivPntSHP)
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$ElevHet)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Geology: Variability on Parental Material.

SoilTypHet <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/10. Soil type heterogenity/SOILType-SOILGRD.tif")
DK.Biodiv$SoilTypHet <- factor(extract(SoilTypHet,DK.BiodivPntSHP))
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$SoilTypHet)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# pH: Soil pH
SoilpH <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/11. pH/SoilpH.tif")
DK.Biodiv$SoilpH <- extract(SoilpH,DK.BiodivPntSHP)
DK.Biodiv$SoilpH[DK.Biodiv$SoilpH==0] <- NA
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$SoilpH)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# NDVI: The Normalized Difference Vegetation Index (NDVI) is an indicator of the greenness of the biomes.
NDVI <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/12. NDVI/NDVI.tif")
DK.Biodiv$NDVI <- extract(NDVI,DK.BiodivPntSHP)
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$NDVI)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Climate
## Mean annual Temperature
MAT <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/13. Climate/DK_BIO01-CHELSA.tif")
DK.Biodiv$MAT <- extract(MAT,DK.BiodivPntSHP.laea)/10
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$MAT)

## Temperature Seasonality
TempSeas <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/13. Climate/DK_BIO04-CHELSA.tif")
DK.Biodiv$TempSeas <- extract(TempSeas,DK.BiodivPntSHP.laea)/10
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$TempSeas)

## Annual Precipitation
AnnPrec <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/13. Climate/DK_BIO12-CHELSA.tif")
DK.Biodiv$AnnPrec <- extract(AnnPrec,DK.BiodivPntSHP.laea)
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$AnnPrec)

## Precipitation Seasonality
PrecSea <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/13. Climate/DK_BIO15-CHELSA.tif")
DK.Biodiv$PrecSea <- extract(PrecSea,DK.BiodivPntSHP.laea)
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$PrecSea)

## Potential Transpiration
PET <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Data/13. Climate/DK_PET-CHELSA.tif")
DK.Biodiv$PET <- extract(PET,DK.BiodivPntSHP.laea)
plot(DK.Biodiv$S.AllGrp~DK.Biodiv$PET)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Save the Dataframe
row.names(DK.Biodiv)<-1:dim(DK.Biodiv)[1]

write.csv(DK.Biodiv,"~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Biodiveristy/Biowide_With_EnvDta.csv")
