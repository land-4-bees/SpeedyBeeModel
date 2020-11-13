#' Calculate landscape forage quality index
#'
#' New implemenatation of landsacpe forage quality index (Lonsdorf et al 2009). Uses FFT convolution to improve runtime.
#' @param output_dir Path to directory for output files.
#' @param landcover_path Path to land cover raster, including base file name
#' @param foragetable_path Path to forage quality by land use table (CDL class integers should be in column called 'LULC')
#' @param forage_table Forage quality by land use table (insead of foragetable_path)
#' @param seasons seasons to include. Must match names of forage table.
#' @param forage_range Foraging range (in m) to use for distance weighting forage scores surrounding focal cell.
#' @param guild_table Bee community to use to model foraging activity. Includes foraging range and relative abundance of each species.
#' 
#' #optional parameters
#' @param agg_factor Aggregation factor for large rasters (use 4 to convert 50m to 120m resolution CDL)
#' @param normalize Normalize values by the number of cells within each moving window?
#' @param rastertag Text string to include in name of output raster
#'
#' @keywords forage index
#' @export
#' @examples
#' 
#' 
forage_index <- function(output_dir, landcover_path, foragetable_path = NA, 
                          forage_table, seasons, forage_range = NA, guild_table = NA, 
                          agg_factor=NA, normalize=T, useW=F, 
                          rastertag=NA, compress_rasters=T) {
    
  #set up logging
  logger::log_threshold(DEBUG)
  logger::log_info('Starting setup.')
  
  #save landscape name, use to label results rasters
  land_name <- gsub(basename(landcover_path), pattern=".tif", replacement="")
  
  #read forage quality table
  if (!is.na(foragetable_path)) {
    forage_table <- read.csv(foragetable_path)
    if(!'LULC' %in% names(forage_table)) {
      stop('Forage table must have column called LULC with integer land use code.')
    }
  }
  
  #check to see if output directory already exists. Make if it doesn't exist.
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  if (is.na(forage_range)) {
    #is the foraging activity for all species in all seasons the same?
    guild <- read.csv(guild_table)
    #scale relative abundance so it sums to one
    guild$relative_abundance <- guild$relative_abundance/sum(guild$relative_abundance)
    
    #calculate abundance weighted average foraging range of all species in guild table
    mean_FR <- sum(guild$alpha*guild$relative_abundance)
    
    #if multiple species in the guild table, summarize as one species
    #it is not necessary to summarize by nesting guild, because nesting is not included
    #this is a variation from Lonsdorf et al (2009), which estimated bee abundance as product of nesting and forage
    if (length(guild[,1]) > 1) {
      guild <- dplyr::select(guild, -dplyr::starts_with('nesting')) %>%
        dplyr::mutate(dplyr::across(dplyr::starts_with('foraging')|alpha, 
                                    function(x){x*relative_abundance})) %>%
        dplyr::summarise(dplyr::across(dplyr::starts_with('foraging')|alpha, sum))
    }
    
  } else {
    mean_FR <- forage_range
  }

  #read land use raster
  hab.r <- raster::raster(landcover_path) 
  
  #does forage table contain the same classes as land cover raster?
  same <- unique(raster::values(hab.r)) %in% forage_table$LULC
  
  #raster land covers that are NOT in landcover table
  missing <- unique(raster::values(hab.r))[!same]
  missing <- missing[!is.na(missing)]
  
  #stop model execution if land cover raster has extra classes not in forage table
  if (length(missing) > 0) {
    stop("Land cover raster has classes that are not in the forage quality table.")
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
  
  logger::log_info('Setup complete, ready to begin distance weighting.')
  
  
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
    
      raster::writeRaster(window_sum, windowsum_path, overwrite=T, options=tifoptions)
    }
  }
  
  for (season in seasons) {
    fcolumn <- names(dplyr::select(forage_table, dplyr::contains(season, ignore.case = T)))
    #reclassify land use to seasonal forage index
    for.r <- raster::reclassify(hab.r, forage_table[,c("LULC", fcolumn)])
    
    #if specified, aggregate forage raster to larger cell size
    if (!is.na(agg_factor)) {
      for.r <- raster::aggregate(for.r, fact = agg_factor, fun = mean)
      
      #strange behavior with applying moving window analysis to raster necessitates writing and reading aggregated rasters again
      raster::writeRaster(for.r, paste0(output_dir, "/for_reclass_agg_", land_name, ".tif"),
                          overwrite=F, NAflag=255, options=tifoptions)
      for.r <- raster::raster(paste0(output_dir, "/for_reclass_agg_", land_name, ".tif"))
  
      #if aggregation factor is supplied, reduce land use raster resolution
      hab.r <- raster::aggregate(hab.r, fact = agg_factor,fun = raster::modal)
      raster::writeRaster(hab.r, paste0(output_dir, "/hab_agg_", land_name, ".tif"),
                          overwrite=F, NAflag=255, options=tifoptions)
      hab.r <- raster::raster(paste0(output_dir, "/hab_agg_", land_name, ".tif"))
    }
  
    #####calculate distance weighted forage index
    
    #distance weighting with FFT convolution
    forage <- raster::as.matrix(for.r)
    
    if (useW == T) {
      FFT_matrix <- smoothie::kernel2dsmooth(x=forage, K=effdist.v, setup=T)
      forage_dw <- smoothie::kernel2dsmooth(x=forage, W=FFT_matrix)
    } else {
      forage_dw <- smoothie::kernel2dsmooth(x=forage, K=effdist.v)
    }
    
    #translate matrix into raster
    simp.for <- raster::raster(forage_dw, template=for.r)
    
    #divide forage raster by window sum to adjust cells on edge of landscape 
    #accounts for less available forage due to lack of resources 'outside' county raster
    if (normalize==T){
      #divide the distance weighted forage raster by the moving window sum
      simp.for <- simp.for/window_sum
      
      #clip distance weighted raster to boundary of land use raster
      simp.for <- simp.for * mask_land
    }
  
    if (!is.na(rastertag)) {
      #write output raster
      raster::writeRaster(simp.for, paste0(output_dir,"/", land_name, "_", season, 
                         "_", rastertag, ".tif"), overwrite=T, options=tifoptions)
    } else {
      #write output raster
      raster::writeRaster(simp.for, paste0(output_dir,"/", land_name, "_", season, 
                                           ".tif"), overwrite=T, options=tifoptions)
    }
  }
  if(normalize == T) {
    rm(hab.r, for.r, forage, forage_dw, simp.for, mask_land, mask, window_sum)
  } else {
    rm(hab.r, for.r, forage, forage_dw, simp.for, weight.m, effdist.v)
  }
  gc()
}
