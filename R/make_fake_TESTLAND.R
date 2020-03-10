
regional_landcover_path <-  './input_data/CDL_2018_42001.tif'

#reassign 'nodata' or NA values. By default CDL uses zeros around outside of map instead of NA
#try 255 as'nodata' value. That is value used by sample LU raster
state_CDL <- raster::raster(regional_landcover_path)
state_CDL_NA <- raster::reclassify(state_CDL, cbind(0, NA)) 
raster::colortable(state_CDL_NA) <- raster::colortable(state_CDL)

#make fake raster with all zeros and one cell with a values of 1.
rastval <- raster::values(state_CDL_NA)
rastval[!is.na(rastval)] <- 1
rastval[length(rastval)/2] <- 2

raster::values(state_CDL_NA) <- rastval
raster::plot(state_CDL_NA)

raster::writeRaster(state_CDL_NA, './input_data/CDL_2018_42001_TESTLAND.tif', NAflag=255, overwrite=T)
