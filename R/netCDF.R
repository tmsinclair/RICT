#Analyse NetCDF files

#Load the ncdf4 package
library(ncdf4)
library(raster)
library(sf)

#read in an example file
ncin <- nc_open(file.choose(), verbose = TRUE, write = FALSE)

print(nc_file)

lon <- ncvar_get(nc_file, "grid_longitude")
nlon <- dim(lon)
head(lon)

time <- ncvar_get(nc_file, "time")
ntime <- dim(time)
head(time)

tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
nt

tunits

nc_close()
ncvar_get(nc_file, "sfcwind")
