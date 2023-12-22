

# Make a mosaic raster of raster tiles

# Data: forest cover over Germany, deciduous vs coniferous forests

# Process:
# read tiles over DE
# clip with Bavaria
# export forest cover deciduous vs coniferous as a raster


library(terra)


r <- rast(ncols=100, nrows=100)
values(r) <- 1:ncell(r)

x <- rast(ncols=2, nrows=2)
b <- rast(ncols=4, nrows=4)

filename <- paste0(tempfile(), "_.tif")
ff <- makeTiles(r, b, filename, overwrite = TRUE)
ff

vrt(ff)


# Have completed the tast in QGI, now I can use the extracted data for Bavaria: forest cover


# Need to rasmple raster to 30 m resolution

library(terra)



