rm(list=ls());gc()
require(terra)

# Load the predictors
setwd("/Users/alejandroordonez/Library/CloudStorage/Dropbox/Courses_&_Conferences/2022/Courses/Statistical and Geospatial Modeling/Lectures-Practical/Week 14 - Project/Data/")
NatArea <- rast("1. Nature Density/NatDens_BaseMap.tif")
plot(NatArea)

DistCoast <- rast("6. Distance to Coast/DISTtoCOAST.tif")
DistCoast <- project(DistCoast,NatArea)
plot(DistCoast)


HII <- rast("7. HumanFootprint/HII.tif")
HII <- project(HII,NatArea)
plot(HII)

Slope <- rast("8. Slope/Slope_30mAgg.tif")
Slope <- project(Slope,NatArea)
plot(Slope)

ElvHet <- rast("9. Elevantion Heterogenity/SRTM30_SD_1km.tif")
ElvHet <- project(ElvHet,NatArea)
plot(ElvHet)


SoilHet <- rast("10. Soil type heterogenity/SOILType-SOILGRD.tif")
SoilHet <- project(SoilHet,NatArea)
plot(SoilHet)

MAT <- rast("13. Climate/DK_BIO01-CHELSA.tif")
MAT <- project(MAT,NatArea)
plot(MAT)

TAP <- rast("13. Climate/DK_BIO12-CHELSA.tif")
TAP <- project(TAP,NatArea)
plot(TAP)


PET <- rast("13. Climate/DK_PET-CHELSA.tif")
PET <- project(PET,NatArea)
plot(PET)


AllVars <- c(NatArea,DistCoast,HII,Slope,ElvHet,SoilHet,MAT,TAP,PET)
#AllVars <- c(NatArea,scale(c(DistCoast,HII,Slope,SoilHet,MAT,TAP)))
names(AllVars) <- c("NatDens_BaseMap","DISTtoCOAST","HII","Slope_30mAgg","SRTM30_SD_1km","SOILType.SOILGRD","DK_BIO01.CHELSA","DK_BIO12.CHELSA","DK_PET.CHELSA")
Temp <- na.omit(AllVars[])

# LM model per parameter
Temp2 <- Temp[sample(dim(Temp)[1],1000),c(6,1)]
plot(lowess(Temp2),
     col = "red",
     lwd = 2,
     typ="b",
     pch=19,
     cex=0.5)
plot(Temp2,pch=19,cex=0.6)
lines(lowess(Temp2), lwd=2,col="red")

Temp2 <- Temp[sample(dim(Temp)[1],10000),]
Temp2 <- data.frame(scale(Temp2))
summary(lm(NatDens_BaseMap~. ,data=Temp2))

# GLM model ALL parameters
a <- lapply(1:10, function(x){
  Temp2 <- Temp[sample(dim(Temp)[1],10000),]
  Temp2 <- data.frame(NatDens_BaseMap=Temp2[,1],scale(Temp2[,-1]))
  #Temp2 <- data.frame(Temp2)
  Mod1 <- glm(NatDens_BaseMap~. ,data=Temp2, family="binomial")
  coef(Mod1)
  Predict <- predict(scale(AllVars),
                     Mod1,
                     type = "response")
  Predict})

plot((mean(do.call("c",a))-NatArea)<0)




# Boxplot model per parameter
AllVars <- c(NatArea>0.9,DistCoast,HII,Slope,ElvHet,SoilHet,MAT,TAP,PET)
#AllVars <- c(NatArea,scale(c(DistCoast,HII,Slope,SoilHet,MAT,TAP)))
names(AllVars) <- c("NatDens_BaseMap","DISTtoCOAST","HII","Slope_30mAgg","SRTM30_SD_1km","SOILType.SOILGRD","DK_BIO01.CHELSA","DK_BIO12.CHELSA","DK_PET.CHELSA")
Temp <- na.omit(AllVars[])
Temp2 <- Temp[,c(2,1)]
anova(lm(Temp2[,1]~Temp2[,2]))
coef(lm(Temp2[,1]~Temp2[,2]))
boxplot(Temp2[,1]~Temp2[,2])

# GLM model ALL parameters
Temp2 <- data.frame(NatDens_BaseMap=Temp[,1],scale(Temp[,-1]))
  Mod1 <- glm(NatDens_BaseMap~. ,data=Temp2, family="binomial")
  coef(Mod1)
  
  Predict <- predict(scale(AllVars),
                     Mod1,
                     type = "response")
  plot(Predict)
hist(Predict)


