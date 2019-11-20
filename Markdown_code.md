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

### 2.5 Cropping and Masking tools

In the following section the script crops and masks the NDVI maps by the field vector shapefile.
For more info about these functions please check the outcome and have a look at the following link
https://rpubs.com/ricardo ochoa/416711

```
# Crop NDVI data 2018
cropped feature 1 = crop(Sentinel NDVI 19.04.2018, field vector)
plot (cropped feature 1)
cropped feature 2 = crop(Sentinel NDVI 07.05.2018, field vector)
plot (cropped feature 2)
cropped feature 3 = crop(Sentinel NDVI 14.05.2018, field vector)
plot (cropped feature 3)
cropped feature 4 = crop(Sentinel NDVI 16.06.2018, field vector)
plot (cropped feature 4)

# Mask NDVI data 2018
masked feature 1 = mask(Sentinel NDVI 19.04.2018, field vector)
plot(masked feature 1)
masked feature 2 = mask(Sentinel NDVI 07.05.2018, field vector)
plot(masked feature 2)
masked feature 3 = mask(Sentinel NDVI 14.05.2018, field vector)
plot(masked feature 3)
masked feature 4 = mask(Sentinel NDVI 16.06.2018, field vector)
plot(masked feature 4)

# Crop NDVI data 2019
cropped feature 5 = crop(Sentinel NDVI 04.04.2019, field vector)
plot (cropped feature 5)
cropped feature 6 = crop(Sentinel NDVI 14.04.2019, field vector)
plot (cropped feature 6)
cropped feature 7 = crop(Sentinel NDVI 27.04.2019, field vector)
plot (cropped feature 7)
cropped feature 8 = crop(Sentinel NDVI 14.05.2019, field vector)
plot (cropped feature 8)

# Mask NDVI data 2019
masked feature 5 = mask(Sentinel NDVI 04.04.2019, field vector)
plot(masked feature 5)
masked feature 6 = mask(Sentinel NDVI 14.04.2019, field vector)
plot(masked feature 6)
masked feature 7 = mask(Sentinel NDVI 27.04.2019, field vector)
plot(masked feature 7)
masked feature 8 = mask(Sentinel NDVI 14.05.2019, field vector)
plot(masked feature 8)
```

### 2.6 Visualization mode: Interactive (view) vs. Static (plot) viewing

Here we switch to an interactive viewing mode of the outcomes (please run this script in R-studio
to check zooming options). If you want to get back to the ”static” mode please switch ’view’ by
’plot’ in the function term ’tmap mode()’.

```
tmap mode(”view”)
# or tmap mode(”plot”) instead
```

### 2.7 Interactive mapping of masked NDVI rasters

```
# For 2018
masked 1 <- tm shape(masked feature 1) + tm raster(palette=”YlGn”,n=10, title=”NDVI 19.04.2018”) +
tm legend(outside = TRUE, text.size = 1.2)
masked 2 <- tm shape(masked feature 2) + tm raster(palette=”YlGn”,n=10, title=”NDVI 07.05.2018”) +
tm legend(outside = TRUE, text.size = 1.2)
masked 3 <- tm shape(masked feature 3) + tm raster(palette=”YlGn”,n=10, title=”NDVI 14.05.2018”) +
tm legend(outside = TRUE, text.size = 1.2)
masked 4 <- tm shape(masked feature 4) + tm raster(palette=”YlGn”,n=10, title=”NDVI 16.06.2018”) +
tm legend(outside = TRUE, text.size = 1.2)

# For 2019
masked 5 <- tm shape(masked feature 5) + tm raster(palette=”YlGn”,n=10, title=”NDVI 04.04.2019”) +
tm legend(outside = TRUE, text.size = 1.2)
masked 6 <- tm shape(masked feature 6) + tm raster(palette=”YlGn”,n=10, title=”NDVI 14.04.2019”) +
tm legend(outside = TRUE, text.size = 1.2)
masked 7 <- tm shape(masked feature 7) + tm raster(palette=”YlGn”,n=10, title=”NDVI 27.04.2019”) +
tm legend(outside = TRUE, text.size = 1.2)
masked 8 <- tm shape(masked feature 8) + tm raster(palette=”YlGn”,n=10, title=”NDVI 14.05.2019”) +
tm legend(outside = TRUE, text.size = 1.2)
```

### 2.8 Print NDVI masked map for 2018

It is possible to zoom each facet separately in R-studio and move the corresponding maps.

```
tmap arrange(masked 1, masked 2, masked 3, masked 4, ncol=2)
```
### 2.9 Print NDVI masked map for 2019

```
tmap arrange(masked 5, masked 6, masked 7, masked 8, ncol=2)
```

### 2.10 Display facets for 2018

1) We convert 2018 raster data into vectorial point-based data while keeping the same spatial
resolution (10x10m); 

2) Use of ’tmap arrange’ to display facets for NDVI vectorial mapping (2018).

```
# From raster to point
NDVI vector 19.04.2018 <- rasterToPoints(masked feature 1, spatial = TRUE) %>% st as sf()
NDVI vector 07.05.2018 <- rasterToPoints(masked feature 2, spatial = TRUE) %>% st as sf()
NDVI vector 14.05.2018 <- rasterToPoints(masked feature 3, spatial = TRUE) %>% st as sf()
NDVI vector 16.06.2018 <- rasterToPoints(masked feature 4, spatial = TRUE) %>% st as sf()

# Correct feature name
names(NDVI vector 19.04.2018)[names(NDVI vector 19.04.2018) == ”layer”] <- ”NDVI 19.04.2018”
names(NDVI vector 07.05.2018)[names(NDVI vector 07.05.2018) == ”layer”] <- ”NDVI 07.05.2018”
names(NDVI vector 14.05.2018)[names(NDVI vector 14.05.2018) == ”layer”] <- ”NDVI 14.05.2018”
names(NDVI vector 16.06.2018)[names(NDVI vector 16.06.2018) == ”layer”] <- ”NDVI 16.06.2018”

# Create maps for 2018
NDVI points 1 <- tm shape(NDVI vector 19.04.2018) + tm dots(col=”NDVI 19.04.2018”, palette=”YlGn”, n=10)
+ tm style(”cobalt”) + tm legend(outside = TRUE, text.size = 1.2)
NDVI points 2 <- tm shape(NDVI vector 07.05.2018) + tm dots(col=”NDVI 07.05.2018”, palette=”YlGn”, n=10)
+ tm style(”cobalt”) + tm legend(outside = TRUE, text.size = 1.2)
NDVI points 3 <- tm shape(NDVI vector 14.05.2018) + tm dots(col=”NDVI 14.05.2018”, palette=”YlGn”, n=10)
+ tm style(”cobalt”) + tm legend(outside = TRUE, text.size = 1.2)
NDVI points 4 <- tm shape(NDVI vector 16.06.2018) + tm dots(col=”NDVI 16.06.2018”, palette=”YlGn”, n=10)
+ tm style(”cobalt”) + tm legend(outside = TRUE, text.size = 1.2)

# Display facets
tmap arrange(NDVI points 1, NDVI points 2, NDVI points 3, NDVI points 4, ncol=4)
```

The previous script is a ’tmap-element’ that specifies facets (small multiples). Small multiples
can be created in two ways: 1) by specifying the by argument with one or two variable names, by which the data is grouped, 2) by specifying multiple variable names in any of the aesthetic argument
of the layer functions (for instance, the argument col in ’tm fill’). For more information please
check https://www.rdocumentation.org/packages/tmap/versions/2.3-1/topics/tm facets.

### 2.11 Display facets for 2019

1) We convert 2019 raster data into vectorial point-based data while keeping the same spatial
resolution (10x10m); 

2) Use of ’tmap arrange’ to display facets for NDVI vectorial mapping (2019).

```
# From raster to point
NDVI vector 04.04.2019 <- rasterToPoints(masked feature 5, spatial = TRUE) %>% st as sf()
NDVI vector 14.04.2019 <- rasterToPoints(masked feature 6, spatial = TRUE) %>% st as sf()
NDVI vector 27.04.2019 <- rasterToPoints(masked feature 7, spatial = TRUE) %>% st as sf()
NDVI vector 14.05.2019 <- rasterToPoints(masked feature 8, spatial = TRUE) %>% st as sf()

# Correct feature name
names(NDVI vector 04.04.2019)[names(NDVI vector 04.04.2019) == ”layer”] <- ”NDVI 04.04.2019”
names(NDVI vector 14.04.2019)[names(NDVI vector 14.04.2019) == ”layer”] <- ”NDVI 14.04.2019”
names(NDVI vector 27.04.2019)[names(NDVI vector 27.04.2019) == ”layer”] <- ”NDVI 27.04.2019”
names(NDVI vector 14.05.2019)[names(NDVI vector 14.05.2019) == ”layer”] <- ”NDVI 14.05.2019”

# Create maps for 2019
NDVI points 5 <- tm shape(NDVI vector 04.04.2019) + tm dots(col=”NDVI 04.04.2019”, palette=”YlGn”, n=10)
+ tm style(”cobalt”) + tm legend(outside = TRUE, text.size = 1.2)
NDVI points 6 <- tm shape(NDVI vector 14.04.2019) + tm dots(col=”NDVI 14.04.2019”, palette=”YlGn”, n=10)
+ tm style(”cobalt”) + tm legend(outside = TRUE, text.size = 1.2)
NDVI points 7 <- tm shape(NDVI vector 27.04.2019) + tm dots(col=”NDVI 27.04.2019”, palette=”YlGn”, n=10)
+ tm style(”cobalt”) + tm legend(outside = TRUE, text.size = 1.2)
NDVI points 8 <- tm shape(NDVI vector 14.05.2019) + tm dots(col=”NDVI 14.05.2019”, palette=”YlGn”, n=10)
+ tm style(”cobalt”) + tm legend(outside = TRUE, text.size = 1.2)

# Display facets
tmap arrange(NDVI points 5, NDVI points 6, NDVI points 7, NDVI points 8, ncol=4)
```

## 3 Geo-spatial analysis of field data

### 3.1 Display sampling and ECa measuring photos

Some photos of the measures taken with a Dualem sensor (i.e. electromagnetic induction sensor),
spatial resolution of 1x15m. Soil ECa was measured at 35 and 85 cm depth before sowing and a
few days after a rainfall event of aprox. 10 mm. Soil samples were also collected to estimate pH
and clay content.

For more info in regard to ECa sensing/mapping please read:

- Johnson, C. K., Doran, J. W., Duke, H. R., Wienhold, B. J., Eskridge, K. M., & Shanahan,
J. F. (2001). Field-scale electrical conductivity mapping for delineating soil condition. Soil Science
Society of America Journal, 65(6), 1829-1837.

- McCutcheon, M. C., Farahani, H. J., Stednick, J. D., Buchleiter, G. W., & Green, T. R.
(2006). Effect of soil water on apparent soil electrical conductivity and texture relationships in a
dryland field. Biosystems Engineering, 94(1), 19-32.

```
# Upload libraries ’imager’ and ’png’
library(knitr) # For knitting document and include graphics function
library(ggplot2) # For plotting
library(png) # For grabbing the dimensions of png files
```
