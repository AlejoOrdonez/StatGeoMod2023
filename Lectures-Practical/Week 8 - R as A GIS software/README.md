# **StaGeoMod2022 - week 2 -[Using `Spatial*` and `Raster` Data in R]**

## Objectives of this week.

1. Load, create, manipulate, and save `Spatial*` object using `R`.
2. Load, create, manipulate/extract, and save `Raster*` object using `R`.
3. Use map algebra to create and manipulate `Raster*` object using `R`.
4. Extract information from `Spatial*` and `Raster*`object using R.

## Using `Spatial*` and `Raster` Data in R in a nutshell

*Text from  <https://edrub.in/ARE212/section12.html>*

Spatial data are increasingly important within research—economics, environmental studies, development, agriculture, political science all are increasingly utilizing spatial variation. Spatial variation is important in many data-generating processes (e.g., the effects of air pollution), and it also can provide some compelling natural experiments (e.g., the rollout of a new policy throughout a country). Furthermore, some topics simply cannot be analyzed well without analyzing the topic in space, e.g., air pollution, gerrymandering, segregation, or well depletion.

Spatial data are fairly distinct from other types of data in their structure. These structural differences tend to require unique techniques for compiling, describing, and analyzing the data. Why the difference? Spatial data have more dimensions than most other datasets that we use. Rather than a value connected to a time, we usually have a value connected to a latitude and a longitude–perhaps coupled with time and/or elevation. This difference is not huge; it just requires a little more machinery.

Of course, there are other tools. ArcGIS is probably the most well-known name in GIS (geographic information systems), but it is proprietary, expensive, and only works with Windows. If you really need it, there are several labs on campus that provide access—see the D-Lab. QGIS provides a free, open-source, cross-platform alternative to ArcGIS. We’ll stick with R today; it’s free, it’s open source, and you are fairly familiar with it. Plus there’s a nice online community of people using R for GIS, so you can (usually/eventually) find solutions when things go wrong. Finally, because you will likely use R for the econometrics in your research, using are for your GIS allows you to minimize the number of programs/scripts required during the course of your analysis

