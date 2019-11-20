# Mapping2Sampe_R
R-script to map points for installing moisture probes 

by Tomás Roquette Tenreiro

Institute for Sustainable Agriculture (IAS-CSIC)
Córdoba, 2019

## 1 Introduction

### 1.1 R Markdown

This is a R Markdown (V3.6) / LaTeX type presentation. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. The generated document includes both
content and output of any embedded R code chunks within the document. For more details see
http://rmarkdown.rstudio.com.

### 1.2 General introduction

The following R-script aims to describe an analytical procedure that combined multiple spatial data
in order to classify management zones for sampling and soil moisture probe installation. The main
goal is to distinguish management zones based on physical (i.e. orientation, elevation, texture,
ECa), chemical (pH), and biological (plant vigor) properties within a crop field. The selected field is
about 9.5 ha, located in the arable region of Cordoba. This document aims also to deliver a script
that can be adjusted to similar analyses or function as a guide to conduct geospatial analysis with
R. The reader can use this document both as a dissemination and a decision supporting tool.

### 1.3 Necessary material - R libraries

The first step consists on updating all necessary libraries for this analysis.

```
install.packages(”rmarkdown”)
install.packages(”dplyr”)
install.packages(”plyr”)
install.packages(”reshape2”)
install.packages(”agricolae”)
install.packages(”quantreg”)
install.packages(”ploty”)
install.packages(”sf”)
install.packages(”raster”)
install.packages(”spData”)
install.packages(’spDataLarge’, repos=’https://nowosad.github.io/drat/’, type=’source’)
install.packages(”gridExtra”)
install.packages(”RColorBrewer”)
install.packages(”root.dir”)
install.packages(”tidyverse”)
install.packages(”ggplot2”)
install.packages(”wesanderson”)
install.packages(”ggpmisc”)
install.packages(”knitr”)
install.packages(”installr”)
install.packages(”lmtest”, repos = ”http://cran.us.r-project.org”)
install.packages(”tinytex”)
```

The second step consists on calling all libraries to the script.
```
library(knitr)
library(sf)
library(dplyr)
library(plyr)
library(reshape2)
library(agricolae)
library(quantreg)
library(plotly)
library(sp)
library(spDataLarge)
library(spData)
library(tmap)
library(raster)
library(gridExtra)
library(quantreg)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(gridExtra)
library(wesanderson)
library(ggpmisc)
library(markdown)
library(tinytex)
library(lmtest)
```
### 1.4 Initial details - working directory

In this section we set initial details to specify the working directory; in this particular case the analysis
was linked to the internal folder "Experimental Catchment Cordoba 19 20" where input and output
data is saved. To run this code please specify the working directory where your input files are saved.

```
rm(list=ls())
getwd()
knitr::opts chunk$set(echo = TRUE)
knitr::opts knit$set(root.dir =
”C:/Users/Tomas R. Tenreiro/Desktop/Experimental Catchment Cordoba 19 20” )
```

### 1.5 Input material

This section uploads all input material. In this particular case, we will work with satellite data
(Sentinel-2), atmospherically corrected (and cloud cover < 4%), that was downloaded from Sentinel-
HUB browser https://apps.sentinel-hub.com/. The script considers imagery from two growing seasons
(i.e. 2017/18 and 2018/19). In order to explore spatial correlation between plant vigor and soil
properties under rain-fed conditions, we focused on late phenological stages before crop senescence.
First year imagery corresponds to a winter wheat crop and second year imagery corresponds to a
rapeseed crop.

```
# Date 1
Sentinel Red 19.04.2018 <- raster(”sentinel/R analysis/Sentinel Red 19.04.2018.tiff”)
Sentinel NIR 19.04.2018 <- raster(”sentinel/R analysis/Sentinel NIR 19.04.2018.tiff”)
# Date 2
Sentinel Red 07.05.2018 <- raster(”sentinel/R analysis/Sentinel Red 07.05.2018.tiff”)
Sentinel NIR 07.05.2018 <- raster(”sentinel/R analysis/Sentinel NIR 07.05.2018.tiff”)
# Date 3
Sentinel Red 14.05.2018 <- raster(”sentinel/R analysis/Sentinel Red 14.05.2018.tiff”)
Sentinel NIR 14.05.2018 <- raster(”sentinel/R analysis/Sentinel NIR 14.05.2018.tiff”)
# Date 4
Sentinel Red 16.06.2018 <- raster(”sentinel/R analysis/Sentinel Red 16.06.2018.tiff”)
Sentinel NIR 16.06.2018 <- raster(”sentinel/R analysis/Sentinel NIR 16.06.2018.tiff”)
# Date 5
Sentinel Red 04.04.2019 <- raster(”sentinel/R analysis/Sentinel Red 04.04.2019.tiff”)
Sentinel NIR 04.04.2019 <- raster(”sentinel/R analysis/Sentinel NIR 04.04.2019.tiff”)
# Date 6
Sentinel Red 14.04.2019 <- raster(”sentinel/R analysis/Sentinel Red 14.04.2019.tiff”)
Sentinel NIR 14.04.2019 <- raster(”sentinel/R analysis/Sentinel NIR 14.04.2019.tiff”)
# Date 7
Sentinel Red 27.04.2019 <- raster(”sentinel/R analysis/Sentinel Red 27.04.2019.tiff”)
Sentinel NIR 27.04.2019 <- raster(”sentinel/R analysis/Sentinel NIR 27.04.2019.tiff”)
# Date 8
Sentinel Red 14.05.2019 <- raster(”sentinel/R analysis/Sentinel Red 14.05.2019.tiff”)
Sentinel NIR 14.05.2019 <- raster(”sentinel/R analysis/Sentinel NIR 14.05.2019.tiff”)
```

## 2 Geo-spatial analysis of satellite NDVI

### 2.1 Estimate NDVI

This section estimates NDVI considering that Satellite B4 corresponds to Red wave length and B8
to NIR wave length. Data is in raster format with a spatial resolution of 10x10 m.

```
# NDVI 2018
Sentinel NDVI 19.04.2018 <- (Sentinel NIR 19.04.2018 - Sentinel Red 19.04.2018) / (Sentinel NIR 19.04.2018 +
Sentinel Red 19.04.2018)
Sentinel NDVI 07.05.2018 <- (Sentinel NIR 07.05.2018 - Sentinel Red 07.05.2018) / (Sentinel NIR 07.05.2018 +
Sentinel Red 07.05.2018)
Sentinel NDVI 14.05.2018 <- (Sentinel NIR 14.05.2018 - Sentinel Red 14.05.2018) / (Sentinel NIR 14.05.2018 +
Sentinel Red 14.05.2018)
Sentinel NDVI 16.06.2018 <- (Sentinel NIR 16.06.2018 - Sentinel Red 16.06.2018) / (Sentinel NIR 16.06.2018 +
Sentinel Red 16.06.2018)
# NDVI 2019
Sentinel NDVI 04.04.2019 <- (Sentinel NIR 04.04.2019 - Sentinel Red 04.04.2019) / (Sentinel NIR 04.04.2019 +
Sentinel Red 04.04.2019)
Sentinel NDVI 14.04.2019 <- (Sentinel NIR 14.04.2019 - Sentinel Red 14.04.2019) / (Sentinel NIR 14.04.2019 +
Sentinel Red 14.04.2019)
Sentinel NDVI 27.04.2019 <- (Sentinel NIR 27.04.2019 - Sentinel Red 27.04.2019) / (Sentinel NIR 27.04.2019 +
Sentinel Red 27.04.2019)
Sentinel NDVI 14.05.2019 <- (Sentinel NIR 14.05.2019 - Sentinel Red 14.05.2019) / (Sentinel NIR 14.05.2019 +
Sentinel Red 14.05.2019)
```

### 2.2 Create maps of NDVI

Here we use ’tm shape’ function applied to raster data (tm raster). It is important to confirm we
have uploaded the library ’tmap’.

```
# Map NDVI 2018
NDVI 19.04.2018 <- tm shape(Sentinel NDVI 19.04.2018) + tm raster(palette=”YlGn”,n=5) + tm legend(outside
= TRUE, text.size = 1.2)
NDVI 07.05.2018 <- tm shape(Sentinel NDVI 07.05.2018) + tm raster(palette=”YlGn”,n=5) + tm legend(outside
= TRUE, text.size = 1.2)
NDVI 14.05.2018 <- tm shape(Sentinel NDVI 14.05.2018) + tm raster(palette=”YlGn”,n=5) + tm legend(outside
= TRUE, text.size = 1.2)
NDVI 16.06.2018 <- tm shape(Sentinel NDVI 16.06.2018) + tm raster(palette=”YlGn”,n=5) + tm legend(outside
= TRUE, text.size = 1.2)
# Map NDVI 2019
NDVI 04.04.2019 <- tm shape(Sentinel NDVI 04.04.2019) + tm raster(palette=”YlGn”,n=5) + tm legend(outside
= TRUE, text.size = 1.2)
NDVI 14.04.2019 <- tm shape(Sentinel NDVI 14.04.2019) + tm raster(palette=”YlGn”,n=5) + tm legend(outside
= TRUE, text.size = 1.2)
NDVI 27.04.2019 <- tm shape(Sentinel NDVI 27.04.2019) + tm raster(palette=”YlGn”,n=5) + tm legend(outside
= TRUE, text.size = 1.2)
NDVI 14.05.2019 <- tm shape(Sentinel NDVI 14.05.2019) + tm raster(palette=”YlGn”,n=5) + tm legend(outside
= TRUE, text.size = 1.2)
```

### 2.3 Visualize NDVI raster maps

First for 2018:

Note: ”ncol” means the number of columns up to a maximum limit of 4 features.

```
tmap arrange(NDVI 19.04.2018, NDVI 07.05.2018, NDVI 14.05.2018, NDVI 16.06.2018, nncol=2)
```
Second for 2019:

```
tmap arrange(NDVI 04.04.2019, NDVI 14.04.2019, NDVI 27.04.2019, NDVI 14.05.2019, ncol=2)
```

### 2.4 Upload the field vector for clip operations

The selected field vector should be the shapefile polygon of the field. In this particular case,
the shapefile was obtained with QGIS considering a topographic characterization of the selected
catchment. The catchment was spatially defined according to a map of flow directions (raster)
obtained with the SAGA - Wang Liu algorithm. The algorithm defines water flow and hillshades
orientation in a particular area of interest from a DEM. We used a DEM with 5m spatial resolution
obtained with LiDAR from CNIG ( http://centrodedescargas.cnig.es/CentroDescargas /index.jsp ).
A catchment was isolated with an area of approximately 9.5 ha.

```
# Upload field vector
field vector <- st read(”sentinel/R analysis/Field vector.shp”)
plot (field vector$geometry)
```

In this step it is also important to check whether the coordinate system is the same as the satellite
imagery; in this particular case we will work with the system WGS84 EPSG 4326.

```
# Check coordinate system
crs(field vector)
crs(Sentinel NDVI 27.04.2019) # as an example to check coordinate system of rasters
```




