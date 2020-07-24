#' Calculate bee landscape indices
#'
#' New implemenatation of landsacpe insecticide toxic load index (Douglas et al.). Uses FFT convolution to improve runtime.
#' @param output_dir Path to directory for output files.
#' @param landcover_path Path to land cover raster, including base file name
#' @param pesticide_path Path to pesticide load by land use table
#' @param forage_range Foraging range (in m) to use for distance weighting insecticide scores surrounding focal cell.
#' @param guild_table Bee community to use to model foraging activity. Includes foraging range and relative abundnace of each species.
#' @param ins_method Insecticide toxic load measure to use. Specify 'oral', 'contact' or 'mean' (mean of 'oral' and 'contact').  
#' @param agg_factor Aggregation factor for large rasters (use 4 to convert 50m to 120m resolution CDL)
#' @param normalize Normalize values by the number of cells within each moving window?
#'
#' @keywords insecticide index
#' @export
#' @examples
#' 
#' 
insecticide_index <- function(output_dir, landcover_path, pesticide_path, 
                              forage_range = NA, guild_table = NA, ins_method='mean',
                              agg_factor=NA, normalize=F, useW=F, check_pesttable=T) {
    
  #save landscape name, use to label results rasters
  land_name <- gsub(basename(landcover_path), pattern=".tif", replacement="")
  
  #read pesticide table
  pesticide_table <- read.csv(pesticide_path)
  
  #check to see if output directory already exists (if so, InVEST model will not overwrite)
  if (dir.exists(output_dir)){
    warning('Output directory already exists. New insecticide raster will not overwrite existing file with the same name.')
  } else {
    dir.create(output_dir)
  }
  
  if (is.na(forage_range)) {
    #is the foraging activity for all species in all seasons the same?
    guild <- read.csv(guild_table)
    #scale relative abundance so it sums to one
    guild$relative_abundance <- guild$relative_abundance/sum(guild$relative_abundance)
  
    # forage_values <- dplyr::select(guild, dplyr::starts_with('foraging_activity'))
    # fvalues <- unique(c(t(forage_values))) #unique values in 'foraging activity' columns
    # if (length(fvalues) == 1) {
    #   stop('The foraging activity constants in the guild table must not all be the same value.')
    # }
    
    #calculate abundance weighted average foraging range of all species in guild table
    mean_FR <- sum(guild$alpha*guild$relative_abundance) 
  } else {
    mean_FR <- forage_range
  }
  
  #read land use raster
  hab.r <- raster::raster(landcover_path) 
  
  #read pesticide table
  pestable <- dplyr::rename(pesticide_table, Value = value)
  
  if (check_pesttable == T) {
    #does pesticide table contain the same classes as land cover raster?
    same <- unique(raster::values(hab.r)) %in% pestable$Value
    
    #raster land covers that are NOT in pesticide table
    missing <- unique(raster::values(hab.r))[!same]
    missing <- missing[!is.na(missing)]
    
    
    #stop model execution if land cover raster has extra classes not in pesticide table
    if (length(missing) > 0) {
      stop("Land cover raster has classes that are not in the insecticide table.")
    }
  }
  
  #choose which column of pesticide table to use to reclassify land use
  if (ins_method == 'oral') {
    pcolumn <- 'ld50_or_ha_bil'
  } else if (ins_method == 'contact') {
    pcolumn <- 'ld50_ct_ha_bil'
  } else if (ins_method == 'mean') {
    pcolumn <- 'ins_int'
    
    #calculate intermediate pesticide index (mean of oral and contact toxicity)
    pestable <- dplyr::mutate(pestable, ins_int = (ld50_ct_ha_bil + ld50_or_ha_bil) /2) %>%
      dplyr::mutate(ins_int= tidyr::replace_na(ins_int, 0))
  }
  
  
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
  
  #reclassify land use to desired pesticide index
  ins.r <- raster::reclassify(hab.r, pestable[,c("Value", pcolumn)])
  
  #if specified, aggregate insecticide raster to larger cell size
  if (!is.na(agg_factor)) {
    ins.r <- raster::aggregate(ins.r, fact = agg_factor, fun = mean)
    
    #strange behavior with applying moving window analysis to raster necessitates writing and reading aggregated rasters again
    raster::writeRaster(ins.r, paste0(output_dir, "/ins_reclass_agg_", land_name, ".tif"), overwrite=F, NAflag=255)
    ins.r <- raster::raster(paste0(output_dir, "/ins_reclass_agg_", land_name, ".tif"))
  }
  
  #if aggregation factor is supplied, reduce land use raster resolution
  if (!is.na(agg_factor)) {
    hab.r <- raster::aggregate(hab.r, fact = agg_factor,fun = raster::modal)
    raster::writeRaster(hab.r, paste0(output_dir, "/hab_agg_", land_name, ".tif"), overwrite=F, NAflag=255)
    hab.r <- raster::raster(paste0(output_dir, "/hab_agg_", land_name, ".tif"))
  }
  
  
  #####calculate distance weighted pesticide index
  
  #try distance weighting with FFT convolution
  ins <- raster::as.matrix(ins.r)
  
  if (useW == T) {
    FFT_matrix <- smoothie::kernel2dsmooth(x=ins, K=effdist.v, setup=T)
    ins_dw <- smoothie::kernel2dsmooth(x=ins, W=FFT_matrix)
  } else {
    ins_dw <- smoothie::kernel2dsmooth(x=ins, K=effdist.v)
  }
  #translate into raster
  simp.ins <- raster::raster(ins_dw, template=ins.r)
  
  
  if (normalize == T) {
    #convert land use raster > 0 into all ones (will use to clip output raster and calc moving window normalization values)
    mask_land <- raster::reclassify(hab.r, cbind(1, max(raster::values(hab.r), na.rm=T), 1))
    
    #calculate moving window weight sums across raster (will be the same value for most of the raster except along edges)
    mask <- raster::as.matrix(mask_land)
    
    if (useW == T) {
      mask_dw <- smoothie::kernel2dsmooth(x=mask, W=FFT_matrix )
    } else {
      mask_dw <- smoothie::kernel2dsmooth(x=mask, K=effdist.v)
    }
    
    #translate moving window sums into raster
    window_sum <- raster::raster(mask_dw, template=mask_land)
    
    #divide the distance weighted insecticide raster by the moving window sum
    simp.ins <- simp.ins/window_sum
    
    #clip distance weighted raster to boundary of land use raster
    simp.ins <- simp.ins * mask_land
  }
  
  #write output raster
  raster::writeRaster(simp.ins, paste0(output_dir,"/", land_name, ".tif"), overwrite=T)

}