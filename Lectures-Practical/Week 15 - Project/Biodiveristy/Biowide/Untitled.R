rm(list=ls());gc()
library(raster)
library(rgdal)
library(maptools)
data(wrld_simpl)

#------------------------------------------------------------------------
#### NOVANA
# mydata <- readRDS("/Users/au467796/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Biodiveristy/NOVANA/NOVANA_NA.Rda")
# head(mydata)
# names(mydata)
# Loc <- unique(mydata[,c("plot","x","y")])
# MeanRich <- aggregate(as.numeric(mydata$antalarter), list(factor(mydata$plot)),mean)
# MeanRich <- MeanRich[match(Loc$plot,MeanRich$Group.1),]
# Loc$Rich <- MeanRich[,2]
# plot(Loc[,c("x","y")],
#      pch=19,
#      cex=c((Loc$Rich/max(Loc$Rich,na.rm=T))*2)
#      )
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#### BIOWIDE
# load the Raw data
DK.Biodiv.In <- readxl::read_excel("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Biodiveristy/Biowide/occurrence.xlsx")
names(DK.Biodiv.In)

# Which groups where sampled?
unique(DK.Biodiv.In$phylum)


# Summarize the richness
DK.BiodivLatLongPaste <- c(paste0(DK.Biodiv.In$decimalLongitude,"-",DK.Biodiv.In$decimalLatitude))
funcCont <- function(x){length(unique(x))}
RichSumm <- tapply(DK.Biodiv.In$scientificName,
                   INDEX = list(DK.BiodivLatLongPaste,
                                DK.Biodiv.In$phylum),
                   FUN = funcCont)

# Make the first summary table focusing on Arthropoda and Plants
DK.Biodiv <- data.frame(ID = 1:130,
                        decimalLongitude = do.call("c",lapply(strsplit(row.names(RichSumm),"-"),function(x){as.numeric(x[1])})),
                        decimalLatitude = do.call("c",lapply(strsplit(row.names(RichSumm),"-"),function(x){as.numeric(x[2])})),
                        S.Arthropoda = RichSumm[,"Arthropoda"], # Species richness Arthropoda
                        S.Spermatophyte = c(apply(RichSumm[,c("Magnoliophyta","Pinophyta")],1,sum, na.rm=T)), # Species richness vascular plant
                        S.Magnoliophyta = RichSumm[,"Magnoliophyta"] # # Species richness flowering plants
                        )

# The first visualization
## plot the richness of Arthropoda in sampled sites
plot(DK.Biodiv[,c("decimalLongitude","decimalLatitude")],
     pch=19,
     cex=c((DK.Biodiv$S.Arthropoda/max(DK.Biodiv$S.Arthropoda))*3),
     main = "Richness Arthropoda")
plot(wrld_simpl,add=T)
legend("topright",
       legend= round(seq(1,max(DK.Biodiv$S.Arthropoda),length.out=5)),
       pch=19,
       pt.cex=c((seq(1,max(DK.Biodiv$S.Arthropoda),length.out=5)/max(DK.Biodiv$S.Arthropoda))*3)
)


## plot the richness of Spermatophyte in sampled sites
plot(DK.Biodiv[,c("decimalLongitude","decimalLatitude")],
     pch=19,
     cex=c((DK.Biodiv$S.Spermatophyte/max(DK.Biodiv$S.Spermatophyte))*3),
     main = "Richness Spermatophyte")
plot(wrld_simpl,add=T)
legend("topright",
       legend= round(seq(1,max(DK.Biodiv$S.Spermatophyte),length.out=5)),
       pch=19,
       pt.cex=c((seq(1,max(DK.Biodiv$S.Spermatophyte),length.out=5)/max(DK.Biodiv$S.Spermatophyte))*3)
       )


#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# Prep to extract env data
## Make a spatial points dataframe
DK.BiodivPntSHP <- SpatialPointsDataFrame(coords = DK.Biodiv[,c("decimalLongitude","decimalLatitude")],
                                          proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),
                                          data = DK.Biodiv)

plot(DK.BiodivPntSHP,
     pch=19,
     cex=c(3*(DK.Biodiv$S.Magnoliophyta/max(DK.Biodiv$S.Magnoliophyta))),
     main = "Richness Spermatophyte")
plot(wrld_simpl,add=T)

# for contrast now add a more detailed DK map
DKMap <- shapefile("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/BaselineShapes/DNK_admin/GADML_DNK0.shp")
plot(DKMap,add=T,border="red")

## Now you are ready to extract some basic positinal data 
## Lest start with extract the Region
landsdeleSHP <- shapefile("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/BaselineShapes/DNK_admin/GADML_DNK1.shp")
DK.Biodiv$landsdele <- over(DK.BiodivPntSHP, # The SpatialPointDataFrame
                            landsdeleSHP["NAME_1"] # The attribute to extract 
                            )
# now lets visualize it

plot(DK.BiodivPntSHP,
     pch=c(15:19)[as.numeric(factor(DK.Biodiv[,"landsdele"][,1]))],
     cex=c(3*(DK.Biodiv$S.Magnoliophyta/max(DK.Biodiv$S.Magnoliophyta))))

plot(landsdeleSHP,
     col=hcl.colors(length(landsdeleSHP[,"NAME_1"][,1]), # number of colors is the number of regions
                    palette = "Roma",alpha = 0.5), # the could palette to use
     
     add=T)

#------------------------------------------------------------------------------------
## Bioscore predictors

# Nature Density:  Nature density is calculated by interpolating the landscape's
#                  nature density calculated as the proportion of forests and 
#                  protected light-open nature types in a nationwide network of
#                  cells of 1 x 1 km.  DONE

#### How to: Summarize the % (or density) of Natural forest areas in Corine at 1km
#Load the data
CorineEU.2012 <- raster("/Volumes/Crucial X6/Data/EUROPE/LAND COVER/PRESENT/Corine Land Cover/100m/1. Land Cover/u2018_clc2012_v2020_20u1_raster100m 2/DATA/U2018_CLC2012_V2020_20u1.tif")

## Re-project the spatial polygon of Denmark to laea
DK.Map.laea <- spTransform(x = DKMap,
                           CRSobj = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))

# crop the data to only DK
CorineDK.Crop <-crop(x = CorineEU.2012, # Corine Map
                     y = extent(DK.Map.laea) # The extent of a shapefile to define the crop 
)

# Mask the data to only DK
CorineDK.2012 <-mask(x = CorineDK.Crop, # Corine Map
                     mask = DK.Map.laea # a shapefile to define what to mask out 
                     )
plot(CorineDK.2012)

# lets check what the values mean
Legend <- read.csv("/Volumes/Crucial X6/Data/EUROPE/LAND COVER/PRESENT/Corine Land Cover/100m/1. Land Cover/u2018_clc2012_v2020_20u1_raster100m 2/Legend/CLC2018_CLC2012_V2018_20_QGIS.txt",header=F)
head(Legend)
## Categories - the First column is a tree number code defining a hierarchy of the levels [see CLC use guide appendix] <https://land.copernicus.eu/user-corner/technical-library/clc-product-user-manual>
# Here I will use the first aggregation
# 1. ARTIFICIAL SURFACES
# 2. AGRICULTURAL AREAS
# 3. FOREST AND SEMI- NATURAL AREAS
# 4. WETLANDS
# 5. WATER BODIES
# 9. NA

NatAreaVal <- 3 # define that natural areas [There have a value of 3 at level 1]

ReclassMrtx <- matrix(c(1:length(as.numeric(substring(Legend[,1],1,1))), # These are the "is" values
                       (as.numeric(substring(Legend[,1],1,1))==NatAreaVal)*1), # These are the "becomes" values
                      ncol = 2, # You define that the matrix has two columns
                      nrow = length(as.numeric(substring(Legend[,1],1,1))) # and that it was the same number of rows as classes
                      )

# Reclassify the raster so that it ONLY shows the areas considered FOREST AND SEMI- NATURAL AREAS 
CorineDK.2012.Nat <- reclassify(x = CorineDK.2012, # raster to reclasify
                                rcl = ReclassMrtx # the Reclassification matrix
                                )
# Now lets look at it
plot(CorineDK.2012.Nat)

# Now we are ready to aggregate the 100x100m base Corine raster into a 1x1km.
# Here you will determine the density of Nature areas per 100km [km of Nature per km]

Summ.Func <- function(x, na.rm=TRUE){sum(x, na.rm=TRUE)/sum(!is.na(x))} # Define the density of points of a class [sum(x)] in all land sites [sum(!is.na(x))]

# Aggregate the 100x100m dataset
CorineDK.2012.Nat1k <- aggregate(x = CorineDK.2012.Nat, # raster to aggregate
                                 fact = 10, # Aggregation factor (number of cells in each direction to merge together)
                                 fun = Summ.Func) # how to aggregate the cells

plot(CorineDK.2012.Nat1k)

# Extract Nature areas per 100km
## Re-porject the spatial points data.frame to UTM zone of Denmark https://epsg.io/25832 https://sdfe.dk/media/2917583/001-etrs89-utm.pdf
DK.BiodivPntSHP.laea <- spTransform(x = DK.BiodivPntSHP,
                                    CRSobj = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))


DK.Biodiv$NatureDens <- extract(CorineDK.2012.Nat1k,DK.BiodivPntSHP.laea)

plot(x = DK.Biodiv$NatureDens,
     y = DK.Biodiv$S.Magnoliophyta,
     pch = 19, 
     col = factor(DK.Biodiv$landsdele[,1]))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#### How to: Summarize the % (or density) of Cells classified as Forest (wet),
#            Nature-dry, and Nature-wet in basemap.
Basemap <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/3. Land Cover/Basemap/Basemap03_public_geotiff/basemap03_2018/lu_agg_2018.tif")

## Re-project the spatial polygon of Denmark to UTM
DK.Map.UTM <- spTransform(x = DKMap,
                          CRSobj = CRS("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_def"))

# The legend of basemap
# Those with level 1 value of 3 are natural areas
BasemapLeg <- readxl::read_excel("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/3. Land Cover/Basemap/Basemap03_public_geotiff/Aggregated Legend.xlsx")
BasemapLeg <- as.data.frame(BasemapLeg)


# Reclassify the raster so that it ONLY shows the areas considered FOREST AND SEMI- NATURAL AREAS 
# the Reclassify matrix for Basemap
a <- unique(getValues(Basemap))
ReclassMrtx <- cbind(as.numeric(a),
                      ifelse(substring(as.vector(a),1,1)==3,1,
                             ifelse(substring(as.vector(a),1,1)%in%c(4,8,9),NA,0)))
# Reclassify Basemap
Basemap.Nat <- reclassify(x = Basemap, # raster to reclasify
                          rcl = ReclassMrtx # the Reclassification matrix
                          )
plot(Basemap.Nat)


# Now we are ready to aggregate the 10x10m basemap raster into a 1x1km.
# Here you will determine the density of Nature areas per 100km [km of Nature per km]

Summ.Func <- function(x, na.rm=TRUE){sum(x, na.rm=TRUE)/sum(!is.na(x))} # Define the density of points of a class [sum(x)] in all land sites [sum(!is.na(x))]

# Aggregate the 100x100m dataset
Basemap.Nat1k <- aggregate(x = Basemap.Nat, # raster to aggregate
                           fact = 100, # Aggregation factor (number of cells in each direction to merge together)
                           fun = Summ.Func) # how to aggregate the cells

plot(Basemap.Nat)

# Extract Nature areas per 100km
## Re-project the spatial points data.frame to UTM zone of Denmark https://epsg.io/25832 https://sdfe.dk/media/2917583/001-etrs89-utm.pdf
DK.BiodivPntSHP.UTM <- spTransform(x = DK.BiodivPntSHP,
                                    CRSobj = CRS("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs"))

DK.Biodiv$NatureDensBaseMap <- extract(Basemap.Nat1k,DK.BiodivPntSHP.UTM)
p <- DK.Biodiv$NatureDensBaseMap
plot(x = log(p/(1-p)),
     y = DK.Biodiv$S.Magnoliophyta,
     pch = 19, 
     col = factor(DK.Biodiv$landsdele[,1]))
#....................................................................................
# Mapped nature:  Density of protected light-open nature types with the natural
#                 forest strategy's mapped layout of biodiversity forest. We have
#                 excluded the selective logging category, but in return included
#                 the state's mapping of §25 forests carried out in 2015-2016.
#                 Furthermore, areas mapped as oak thickets are included,
#                 cf. Section 26 of the Forest Act. The proxy is important for 
#                 all species groups and land types, but especially for plants
#                 and for open land. DONE

#### How to: Summarize the % of Natural open areas in Corine at 1km

#Load the data
CorineEU.2012 <- raster("/Volumes/Crucial X6/Data/EUROPE/LAND COVER/PRESENT/Corine Land Cover/100m/1. Land Cover/u2018_clc2012_v2020_20u1_raster100m 2/DATA/U2018_CLC2012_V2020_20u1.tif")

## Re-project the spatial polygon of Denmark to laea
DK.Map.laea <- spTransform(x = DKMap,
                           CRSobj = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))

# crop the data to only DK
CorineDK.Crop <-crop(x = CorineEU.2012, # Corine Map
                     y = extent(DK.Map.laea) # The extent of a shapefile to define the crop 
)

# Mask the data to only DK
CorineDK.2012 <-mask(x = CorineDK.Crop, # Corine Map
                     mask = DK.Map.laea # a shapefile to define what to mask out 
)
plot(CorineDK.2012)

# lets check what the values mean
Legend <- read.csv("/Volumes/Crucial X6/Data/EUROPE/LAND COVER/PRESENT/Corine Land Cover/100m/1. Land Cover/u2018_clc2012_v2020_20u1_raster100m 2/Legend/CLC2018_CLC2012_V2018_20_QGIS.txt",header=F)
head(Legend)
## Categories - the First column is a tree number code defining a hierarchy of the levels [see CLC use guide appendix] <https://land.copernicus.eu/user-corner/technical-library/clc-product-user-manual>
NatAreaVal <- c(32,33) # define that OPEN natural areas are 1 in the level 2

ReclassMrtx <- cbind(c(1:dim(Legend)[1]),
                     c((substring(Legend[,1],1,2)%in%NatAreaVal)*1))

# Reclassify the raster so that it ONLY shows the areas considered FOREST AND SEMI- NATURAL AREAS 
CorineDK.2012.OpnNat <- reclassify(x = CorineDK.2012, # raster to reclasify
                                   rcl = ReclassMrtx # the Reclassification matrix
)
# Now lets look at it
plot(CorineDK.2012.OpnNat)

# Now we are ready to aggregate the 100x100m base Corine raster into a 1x1km.
# Here you will determine the density of Nature areas per 100km [km of Nature per km]

Summ.Func <- function(x, na.rm=TRUE){sum(x, na.rm=TRUE)/sum(!is.na(x))} # Define the density of points of a class [sum(x)] in all land sites [sum(!is.na(x))]

# Aggregate the 100x100m dataset
CorineDK.2012.OpnNatNat1k <- aggregate(x = CorineDK.2012.OpnNat, # raster to aggregate
                                       fact = 10, # Aggregation factor (number of cells in each direction to merge together)
                                       fun = Summ.Func) # how to aggregate the cells

plot(CorineDK.2012.OpnNatNat1k)

# Extract open Nature areas per 1km
## Re-porject the spatial points data.frame to UTM zone of Denmark https://epsg.io/25832 https://sdfe.dk/media/2917583/001-etrs89-utm.pdf
DK.BiodivPntSHP.laea <- spTransform(x = DK.BiodivPntSHP,
                                    CRSobj = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))


DK.Biodiv$OpnNatureDens <- extract(CorineDK.2012.OpnNatNat1k,DK.BiodivPntSHP.laea)

plot(x = log(p/(1-p)),
     y = DK.Biodiv$S.Magnoliophyta,
     pch = 19, 
     col = factor(DK.Biodiv$landsdele[,1]))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#### How to: Summarize the % (or density) of Cells classified as
#            Nature-dry, and Nature-wet in basemap.
Basemap <- raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/3. Land Cover/Basemap/Basemap03_public_geotiff/basemap03_2018/lu_agg_2018.tif")
## Re-project the spatial polygon of Denmark to UTM
DK.Map.UTM <- spTransform(x = DKMap,
                          CRSobj = CRS("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_def"))
# The legend of basemap
# Those with level 1 value of 3 are natural areas
BasemapLeg <- readxl::read_excel("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/3. Land Cover/Basemap/Basemap03_public_geotiff/Aggregated Legend.xlsx")
BasemapLeg <- as.data.frame(BasemapLeg)


# Reclassify the raster so that it ONLY shows the areas considered open SEMI- NATURAL AREAS 
# the Reclassify matrix for Basemap
a <- unique(getValues(Basemap))

ReclassMrtx <- cbind(as.numeric(a),
                     ifelse(substring(as.vector(a),1,2)==32,1,
                            ifelse(substring(as.vector(a),1,1)%in%c(4,8,9),NA,0)))
# Reclassify Basemap
Basemap.OpenNat <- reclassify(x = Basemap, # raster to reclasify
                              rcl = ReclassMrtx # the Reclassification matrix
)
plot(Basemap.OpenNat)


# Now we are ready to aggregate the 10x10m basemap raster into a 1x1km.
# Here you will determine the density of Nature areas per 100km [km of Nature per km]

Summ.Func <- function(x, na.rm=TRUE){sum(x, na.rm=TRUE)/sum(!is.na(x))} # Define the density of points of a class [sum(x)] in all land sites [sum(!is.na(x))]

# Aggregate the 100x100m dataset
Basemap.OpenNatNat1k <- aggregate(x = Basemap.OpenNat, # raster to aggregate
                                  fact = 100, # Aggregation factor (number of cells in each direction to merge together)
                                  fun = Summ.Func) # how to aggregate the cells

plot(Basemap.OpenNatNat1k)

# Extract Nature areas per 100km
## Re-porject the spatial points data.frame to UTM zone of Denmark https://epsg.io/25832 https://sdfe.dk/media/2917583/001-etrs89-utm.pdf
DK.BiodivPntSHP.UTM <- spTransform(x = DK.BiodivPntSHP,
                                   CRSobj = CRS("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs"))

DK.Biodiv$OpnNatureDensBaseMap <- extract(Basemap.OpenNatNat1k,DK.BiodivPntSHP.UTM)
p <- DK.Biodiv$OpnNatureDensBaseMap
plot(x = log(p/(1-p)),
     y = DK.Biodiv$S.Magnoliophyta,
     pch = 19, 
     col = factor(DK.Biodiv$landsdele[,1]))

#....................................................................................
#Forest structure: New proxy that has been developed on the basis of nationwide
#                   LiDAR data for the height of the crown roof. The proxy has
#                   been developed to reflect the variation in crown height
#                   within a 50 m radius (Groom et al. 2018).

#### How to: Summarize variation in canopy heigh at 1km - use https://glad.umd.edu/dataset/gedi





#....................................................................................
# Low-lying land: Points for areas that lie on low-lying land, but not intensively
#                 cultivated fields.

#### How to: Summarize the % of wetlands open areas in Corine at 1km.
#            alternative use Using 1km means of the SRTM 90m database estimate the areas in the
#            first quantile or bewteen -7 and 10m. --> will be categorical

#....................................................................................
# Distance to the coast: Areas that are less than 1 km from the coast. As for the
#                        other proxies, no points are given for intensively
#                       cultivated fields.

#....................................................................................
# Nitrogen deposition:  Based on model calculations for atmospheric deposition of
#                       nitrogen, carried out as part of the National Program for
#                       Monitoring the Water Environment and Nature, NOVANA (Ellermann et al. 2019).
#                       The deposit is calculated in 400 x 400 meter grid cells.
#                       Average deposition for the grid cells has been calculated
#                       based on the vegetation types found in each grid cell,
#                       based on the area categories from Basemap03 (Levin 2019),
#                       compared with the deposition model's deposition levels for
#                       different vegetation types. Proxy points are given to cells
#                       with nitrogen deposition of <10 kg/ha

### How to: Use the nitrogen deposition map of Tan et al (2018)
# https://thredds.met.no/thredds/catalog/data/EMEP/Articles_data/Schwede_etal_Ndep_2018/catalog.html

#....................................................................................
# Line density: The density of man-made lines in the landscape such as roads,
#               ditches and field boundaries. The proxy gives points for areas
#               with less than 8 km of lines per 500 x 500 m.

### How to: Use the Human Footprint index <https://sedac.ciesin.columbia.edu/data/set/wildareas-v3-2009-human-footprint>

#....................................................................................
# Habitat nature: Areas that are mapped as one of the habitat directive's protected
#                 nature types according to Annex I.

### How to: use the Natura2000 areas shapefile (https://www.eea.europa.eu/data-and-maps/data/natura-13/natura-2000-spatial-data/natura-2000-shapefile-1)
#   or https://www.arcgis.com/home/item.html?id=fa8b41fdc9e74841a054a4946e944717 also look at
# https://www.eea.europa.eu/data-and-maps/data/article-17-database-habitats-directive-92-43-eec-2


#------------------------------------------------------------------------
## Not in Bioscore 2020

# Slopes: are defined as areas that, after modeling the digital elevation model
#         (http://download.kortforsyningen.dk/content/dhmterr%C3%A6n-16-m-grid)
#         on a 9.6 x 9.6 m scale, has an average slope of more than 15 degrees.
#         The layer is set to 0 for the places where intensive fields overlap
#         with the slopes. Sloping terrain is included as a potential indicator
#         of biodiversity because sloping areas are more difficult to cultivate
#         than flat terrain, and thus have an increased probability of containing
#         important habitats.

#### How to: Using 1km means of the SRTM 90m database estimate the slope

#....................................................................................
# Forest Type: Type of forest vegetation.
#### How to: Check Copernicus

#....................................................................................
# Cultivated area: Density of cultivated areas.

#### How to: Summarize the % of cultivated areas in Corine at 1km.

#....................................................................................
# Geology: Parental Material.
#### How to: Extract the parental material from Denmark soil map.

#....................................................................................
# Soils: pH.
#### How to: Extract the pH from SoilGrids.

#....................................................................................
# Productivity: NDVI 
#### How to: Extract NDVI values from <https://land.copernicus.eu/global/products/ndvi>

#....................................................................................
# Productivity: Vegetation productivity 
#### How to: Extract productivity values from <https://www.eea.europa.eu/data-and-maps/data/annual-above-ground-vegetation-productivity>










## extract Surface_Geology
Surface_GeologySHP <- shapefile("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/7. soils/Jordartskort 1_200000/Jordart_200000_Shape/Jordart_200000.shp")
Surface_GeologySHPProj <- spTransform(Surface_GeologySHP, CRSobj = CRS("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs"))
DK.Biodiv$Surface_Geology <- over(DK.BiodivPntSHPProj,Surface_GeologySHPProj["TSYM"])

# Land cover - Corine
Corine <- raster("/Volumes/Seagate Expansion/Data/EUROPE/LAND COVER/PRESENT/Corine Land Cover/v2020/2018/DATA/U2018_CLC2018_V2020_20u1.tif")
DK.BiodivPntSHPProj2 <- spTransform(DK.BiodivPntSHP, CRSobj = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))

DK.Biodiv$LansUse <- factor(read.table("/Volumes/Seagate Expansion/Data/EUROPE/LAND COVER/PRESENT/Corine Land Cover/v2020/2018/Legend/CLC2018_CLC2018_V2018_20_QGIS.txt",sep=",")[,6][DK.Biodiv$LansUse])

Corine2 <- crop(Corine,extent(DK.BiodivPntSHPProj2))
plot(Corine2)
plot(DK.BiodivPntSHPProj2,pch=18,add=T)

# Land cover - basemap
basemap <- raster("/Users/alejandroordonez/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/3. Land Cover/Basemap/Basemap03_public_geotiff/basemap03_2018/lu_agg_2018.tif")
DK.Biodiv$LansUse <- 
unique(as.factor(extract(basemap,DK.BiodivPntSHPProj)))
142000 # Road, not paved
211000 # Agriculture, intensive, temporary crops
311000 # Forest
312000 # Forest wet
321000 # Nature dry
321220 # Nature dry & Agriculture extensive
322000 # Nature wet
322220 # Nature wet & Agriculture extensive
411000 # Lake


## extract Climate
Climate <- stack(list.files("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/2. Climate",full.names = T))
DK.Biodiv <- data.frame(DK.Biodiv,extract(Climate,DK.BiodivPntSHPProj))

## extract Human Influence
DK.Biodiv <- data.frame(DK.Biodiv,
                        HII = extract(raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/5. Anthropogenic conditions/DK_HII.tif"),DK.BiodivPntSHPProj))

## extract Distance to Coast
DK.Biodiv <- data.frame(DK.Biodiv,
                        DISTtoCOAST = extract(raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/6. Topographic features/DK_DISTtoCOAST.tif"),DK.BiodivPntSHPProj))

## extract Elevation Range
DK.Biodiv <- data.frame(DK.Biodiv,
                        ELEVRng = extract(raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/4. Landscape heterogeneity/DK_ELEVRng-SRTM.tif"),DK.BiodivPntSHPProj))

## extract Nat Veg Heterogeneity
DK.Biodiv <- data.frame(DK.Biodiv,
                        VegHet = extract(raster("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/4. Landscape heterogeneity/DK_VegHet-CLC.tif"),DK.BiodivPntSHPProj))

## extract Land Cover
LandCov <- stack(list.files("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/Rasters/3. Land Cover",full.names = T))
DK.Biodiv <- data.frame(DK.Biodiv,extract(LandCov,DK.BiodivPntSHPProj))

names(DK.Biodiv) <- gsub("DK_","",names(DK.Biodiv))7 geomorpho

### BASE CARTOGRAPHY
# Get the "landsdele" data - from DWA
# Here I define the url to get data from DAWA/DAGI for landsdele. I define that the format should be geojson
#url = "https://api.dataforsyningen.dk/landsdele?format=geojson"
# I will download the geojson to a temporary file
#download.file(url, "~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/BaselineShapes/landsdele")

# Then, I use the rgdal package to read in the downloaded file
#geodata <- rgdal::readOGR("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/BaselineShapes/landsdele")
#saveRDS(geodata, file = "~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/BaselineShapes/landsdele.rds")
#setwd("~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/BaselineShapes")
#landsdele.shp <- getData('GADM', country='DNK', level=1)
#shapefile(DKGADML3,"~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/BaselineShapes/landsdele.shp")

write.csv(DK.Biodiv,"~/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Project/DK.Biodiv.csv")