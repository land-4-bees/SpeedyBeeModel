#' Calculate landscape nesting quality index
#'
#' Variation of landscape nesting quality index (Lonsdorf et al 2009). Uses FFT convolution to improve runtime. 
#' @param output_dir Path to directory for output files.
#' @param landcover_path Path to land cover raster, including base file name
#' @param habitattable_path Path to habitat quality by land use table (CDL class integers should be in column called 'LULC')
#' @param habitat_table Habitat quality by land use table (instead of habitattable_path)
#' @param nest_locations Nesting types to include. Must match names of habitat table.
#' @param forage_range Foraging range (in m) to use for distance weighting scores surrounding focal cell.
#' @param guild_table Bee community to use to model foraging activity. Includes foraging range and relative abundance of each species.
#' 
#' #optional parameters
#' @param agg_factor Aggregation factor for large rasters (use 4 to convert 50m to 120m resolution CDL)
#' @param normalize Normalize values by the number of cells within each moving window?
#' @param rastertag Text string to include in name of output raster
#' @param verbose Include more log messages from model run?
#'
#' @keywords nesting index
#' @export
#' @details 
#' It is necessary to specify 'forage_range' OR 'guild_table,' not both.
#'
#' 
nesting_index <- function(output_dir, landcover_path, habitattable_path = NA, 
                          habitat_table, nest_locations, forage_range = NA, guild_table = NA, 
                          agg_factor=NA, normalize=T, useW=F, 
                          rastertag=NA, verbose=T) {
  library(logger)
  #set up logging
  logger::log_threshold(DEBUG)
  if (verbose == T) {logger::log_info('Starting setup.')}
  
  #save landscape name, use to label results rasters
  land_name <- gsub(basename(landcover_path), pattern=".tif", replacement="")
  
  #read habitat table
  if (!is.na(habitattable_path)) {
    habitat_table <- read.csv(habitattable_path)
    if(!'LULC' %in% names(habitat_table)) {
      stop('Habitat table must have column called LULC with integer land use code.')
    }
  }
  
  #check to see if output directory already exists. Make if it doesn't exist.
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  if (is.na(forage_range)) {
    #is the foraging activity for all species in all nest_locations the same?
    guild <- read.csv(guild_table)
    #scale relative abundance so it sums to one
    guild$relative_abundance <- guild$relative_abundance/sum(guild$relative_abundance)
    
    #calculate abundance weighted average foraging range of all species in guild table
    mean_FR <- sum(guild$alpha*guild$relative_abundance)
    
    #if multiple species in the guild table, summarize as one species
    #it is not necessary to summarize by nesting guild, because nesting is not included
    #this is a variation from Lonsdorf et al (2009), which estimated bee abundance as product of nesting and forage
    if (length(guild[,1]) > 1) {
      guild <- dplyr::select(guild, -dplyr::starts_with('foraging')) %>%
        dplyr::mutate(dplyr::across(dplyr::starts_with('nesting')|alpha, 
                                    function(x){x*relative_abundance})) %>%
        dplyr::summarise(dplyr::across(dplyr::starts_with('nesting')|alpha, sum))
    }
    
  } else {
    mean_FR <- forage_range
  }

  #read land use raster
  hab.r <- raster::raster(landcover_path) 
  
  #does habtiat table contain the same classes as land cover raster?
  same <- unique(raster::values(hab.r)) %in% habitat_table$LULC
  
  #raster land covers that are NOT in landcover table
  missing <- unique(raster::values(hab.r))[!same]
  missing <- missing[!is.na(missing)]
  
  #warn if land cover raster has extra classes not in habitat table
  if (length(missing) > 0) {
    stop("Land cover raster has classes that are not in the habitat table.")
  }
  
  #####set up moving window
  
  #set a maximum foraging distance to be twice the foraging range
  #will be used later to determine size of moving window
  maxforage.dist <- mean_FR*2
  
  ###set up moving window specifications to be used for distance weighting later
  #store cell size of raster
  if (!is.na(agg_factor)) {
    c.size <- raster::res(hab.r)[1] * agg_factor
  } else {
    c.size <- raster::res(hab.r)[1] }
  
  #matrix boundary for moving window (in number of pixels)
  radius <- round(maxforage.dist/c.size)
  
  # size of moving window is set to twice the radius plus one
  # create matrices for moving window distance matrix
  weight.m <- dist.m <- matrix(data=1, nrow= radius*2+1, ncol= radius*2+1)
  focal.c <- median(1:nrow(dist.m)) # focal.c is row number of focal cell in the distance matrix
  
  # calculating distance all cells from the focal cell in the nested matrix
  # loops through unique combination of all row and column numbers (cells), uses pythag theorum to calculate distance from cell center to focal cell center
  for(i in 1:nrow(dist.m)) {
    for (j in 1:ncol(dist.m)) {
      
      dist.m[i,j] <- sqrt( ((i-0.5)*c.size - (focal.c - 0.5)*c.size)^2 +
                             ((j-0.5)*c.size - (focal.c - 0.5)*c.size)^2 )
    }
  }
  
  effdist.v <- exp(-dist.m / mean_FR) # calculate effective distance (moving window weights)
  
  # set 0, when effdist > (2xforaging distance)
  effdist.v[which(dist.m > maxforage.dist)] <- 0
  
  if (verbose == T) {logger::log_info('Setup complete, ready to begin distance weighting.')}
  
  
  #set up window sum raster (used in normalization later)
  if (normalize == T) {
    
    #convert land use raster > 0 into all ones (will use to clip output raster and calc moving window normalization values)
    mask_land <- raster::reclassify(hab.r, cbind(1, max(raster::values(hab.r), na.rm=T), 1))
    
    windowsum_path <- paste0(output_dir,"/intermediate/", land_name,"_windowsum.tif")
    
    if (file.exists(windowsum_path)) {
      window_sum <- raster::raster(windowsum_path)
    } else {
      
      #calculate moving window weight sums across raster (will be the same value for most of the raster except along edges)
      mask <- raster::as.matrix(mask_land)
      
      if (useW == T) {
        mask_dw <- smoothie::kernel2dsmooth(x=mask, W=FFT_matrix )
      } else {
        mask_dw <- smoothie::kernel2dsmooth(x=mask, K=effdist.v)
      }
      
      #translate moving window sums into raster
      window_sum <- raster::raster(mask_dw, template=mask_land)
      
      if (!dir.exists(paste0(output_dir,"/intermediate"))) {
        dir.create(paste0(output_dir,"/intermediate"))
      }
    
      raster::writeRaster(window_sum, windowsum_path, overwrite=T)
      if (verbose == T) {logger::log_info('Generated window-sum raster. First FFT complete.')}
    }
  }
  
  if (verbose == T) {logger::log_info('Starting nest index calculations for each nesting type. Reclass land use to appropriate nesting value.')}
  
  for (nest_location in nest_locations) {
    fcolumn <- names(dplyr::select(habitat_table, dplyr::contains(nest_location, ignore.case = T)))
    
    #reclassify land use to nesting index
    nest.r <- raster::reclassify(hab.r, habitat_table[,c("LULC", fcolumn)])
    
    #if specified, aggregate nesting raster to larger cell size
    if (!is.na(agg_factor)) {
      nest.r <- raster::aggregate(nest.r, fact = agg_factor, fun = mean)
      
      #strange behavior with applying moving window analysis to raster necessitates writing and reading aggregated rasters again
      raster::writeRaster(nest.r, paste0(output_dir, "/nest_reclass_agg_", land_name, ".tif"),
                          overwrite=F, NAflag=255)
      nest.r <- raster::raster(paste0(output_dir, "/nest_reclass_agg_", land_name, ".tif"))
  
      #if aggregation factor is supplied, reduce land use raster resolution
      hab.r <- raster::aggregate(hab.r, fact = agg_factor,fun = raster::modal)
      raster::writeRaster(hab.r, paste0(output_dir, "/hab_agg_", land_name, ".tif"),
                          overwrite=F, NAflag=255)
      hab.r <- raster::raster(paste0(output_dir, "/hab_agg_", land_name, ".tif"))
    }
  
    #####calculate distance weighted forage index
    
    if (verbose == T){logger::log_info('Translate raster to matrix and begin distance-weighting of nest index.')}
    
    #distance weighting with FFT convolution
    nesting <- raster::as.matrix(nest.r)
    
    if (useW == T) {
      FFT_matrix <- smoothie::kernel2dsmooth(x=nesting, K=effdist.v, setup=T)
      nesting_dw <- smoothie::kernel2dsmooth(x=nesting, W=FFT_matrix)
    } else {
      nesting_dw <- smoothie::kernel2dsmooth(x=nesting, K=effdist.v)
    }
    
    #translate matrix into raster
    simp.for <- raster::raster(nesting_dw, template=nest.r)
    
    #divide nesting raster by window sum to adjust cells on edge of landscape 
    #accounts for less available nesting due to lack of resources 'outside' county raster
    if (normalize==T){
      #divide the distance weighted nesting raster by the moving window sum
      simp.for <- simp.for/window_sum
      
      #clip distance weighted raster to boundary of land use raster
      simp.for <- simp.for * mask_land
    }
    if (verbose== T){logger::log_info('One nesting map is complete!')}
    
    if (!is.na(rastertag)) {
      #write output raster
      raster::writeRaster(simp.for, paste0(output_dir,"/", land_name, "_", nest_location, 
                         "_", rastertag, ".tif"), overwrite=T)
    } else {
      #write output raster
      raster::writeRaster(simp.for, paste0(output_dir,"/", land_name, "_", nest_location, 
                                           ".tif"), overwrite=T)
    }
  }
  if(normalize == T) {
    rm(hab.r, nest.r, nesting, nesting_dw, simp.for, mask_land, window_sum)
  } else {
    rm(hab.r, nest.r, nesting, nesting_dw, simp.for, weight.m, effdist.v)
  }
  gc()
  if (verbose ==T){logger::log_info('All nesting maps are complete!')}
  
}
