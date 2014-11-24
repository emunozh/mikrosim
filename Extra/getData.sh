#!/bin/sh
## A simple scrip to document the data extraction and preparation to run a
## microsimulation in Germany
## Esteban Munoz H. emunozh@gmail.com

## Make folder to store the data
mkdir ../Data
## Go into that folder
cd ../Data
## Download the shapefiles for visualization
wget https://www.zensus2011.de/SharedDocs/Downloads/DE/Shapefile/VG250_1Jan2011_UTM32.zip?__blob=publicationFile&v=25
## Download some ready to go data from the census 2011
wget https://www.zensus2011.de/SharedDocs/Downloads/DE/Pressemitteilung/DemografischeGrunddaten/2G_BevoelkerungStaatsangehoerigkeitGeschlecht.zip?__blob=publicationFile&v=4
wget https://www.zensus2011.de/SharedDocs/Downloads/DE/Pressemitteilung/DemografischeGrunddaten/2F_BevoelkerungAlterGeschlecht.zip?__blob=publicationFile&v=6

## Make a folder to decompress the shapefiles
mkdir Shapefiles
## Decompress the file
unzip VG250_1Jan2011_UTM32.zip -d Shapefiles
## Remove the zip file
rm VG250_1Jan2011_UTM32.zip

## Make a folder to store the data
mkdir Gemeinden
## Decompress the folders
unzip 2G_BevoelkerungStaatsangehoerigkeitGeschlecht.zip -d Gemeinden
unzip 2F_BevoelkerungAlterGeschlecht.zip -d Gemeinden
## Remove the zip files
rm 2G_BevoelkerungStaatsangehoerigkeitGeschlecht.zip
rm 2F_BevoelkerungAlterGeschlecht.zip

## Download some more data that requires some more attention, we use python to
## Do this. It used python3 with libraries: pandas and urllib. This scrip need
## Modification to be used by a lower version of python.
python cleanAGS.py

## Make a folder to download the data
mkdir Survey
## Download the micro census as csv file.
wget http://www.forschungsdatenzentrum.de/bestand/mikrozensus/cf/2002/fdz_mikrozensus_cf_2002_ascii-csv.zip
## Decompress the data
unzip fdz_mikrozensus_cf_2002_ascii-csv.zip -d Survey
## Remove the zip file
rm fdz_mikrozensus_cf_2002_ascii-csv.zip
