This repository contains scripts and functions for analyzing the occurrence of *Amaranthus tuberculatus* (common waterhemp) and its surrounding land cover characteristics using geospatial and ecological data. The project integrates data from sources like iNaturalist, [GBIF](https://www.gbif.org/), and the [National Land Cover Database](https://www.usgs.gov/centers/eros/science/national-land-cover-database) to assess land-use changes and species distribution. 

The notebooks/ directory contains notebooks (specifically, zonal_stats.ipynb) that detail the process of overlaying *Amaranthus tuberculatus* observation data with National Land Cover Database land cover raster files which you can download from the [MRLC](https://www.mrlc.gov/data?f%5B0%5D=category%3ALand%20Cover&f%5B1%5D=project_tax_term_term_parents_tax_term_name%3AAnnual%20NLCD&f%5B2%5D=region%3Aconus) website or [NLCD viewer](https://www.mrlc.gov/viewer/), which allows for more specific geography selection.

The scripts/ directory has modules to fetch and format the locations of open-location iNaturalist observations of *Amaranthus tuberculatus* and preserved specimen records from GBIF within a given bounding box, in additon to modules to calculate Simpson and Shannon diversity indices of a list of population proportions.


#### Data Sources:
- iNaturalist: Citizen science for recording species occurrence with geolocation and other data.

- GBIF: Global Biodiversity Information Facility dataset.

- NLCD: Land cover classification rasters from NASA/USGS.
