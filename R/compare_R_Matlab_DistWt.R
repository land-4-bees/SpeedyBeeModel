
state <- 'AdamsPA'
year <- 2018
insyear <- 2014
forage_range <- '3km'


#create raster that is all zero except for one pixel equal to 1

  #check that foraging range input is either 5km or 3km
  if (!forage_range %in% c('5km', '3km')) {
    stop('Foraging range not specified correctly. Must be either \'3km\' or \'5km\'.')
  }
  
  #arguments for pesticide and floral/nest landscape indices calculation (replace these with path names to downloaded data)
  output_dir <- './output/' #change output_dir in line 54!!!
  regional_landcover_path <-  './input_data/CDL_2018_42001_TESTLAND.tif'
  landcover_table <- './input_data/Koh_LU_table_APIS_TESTLAND.csv'
  if (forage_range == '5km') {
    guild_table <- 'D:/Documents/StateIndicesAzavea/landmodel_inputs/Apis2500_guild.csv'
  } else if (forage_range == '3km') {
    guild_table <- 'D:/Documents/StateIndicesAzavea/landmodel_inputs/Apis1500_guild.csv'
  }
  ins_method <- 'mean'
  raster_types <- c('floral_resources', 'nesting_distwt')
  pest_raster <- 'simple' #, 'floral_weighted'), 'normal_floral_weighted')
  
  out_nest <- T
  lg_raster <- T
  agg_factor <- NA
  
  #pesticide tables
  pesticide <- read.csv(paste0('./input_data/PA_CDL_PestReclass', insyear,'_TESTLAND.csv'), stringsAsFactors = F)
  
  part1end <- Sys.time()
  
  landcover_path <-regional_landcover_path 
  
  #store info on raster cell size
  lc <- raster::raster(landcover_path);  c.size <- raster::res(lc)[1] * agg_factor
  
  if (forage_range == '5km') {
    output_dir  <- paste0('./', state, "_TESTLAND_5km_res", c.size, "_", year, '/')
  } else if (forage_range == '3km') {
    output_dir  <- paste0('./', state, "_TESTLAND_3km_res", c.size, "_", year, '/')
  }
  
  if (!dir.exists(output_dir)) {
    #run honey bee InVEST model
    InVESTwrap::beeland_indices(landcover_path=landcover_path, pest_raster=pest_raster, 
                                output_dir=output_dir, landcover_table=landcover_table, 
                                pesticide_table=pesticide, guild_table=guild_table, ins_method=ins_method, 
                                lg_raster=lg_raster, agg_factor=agg_factor, out_nest=out_nest)
    
    part2end <- Sys.time()
    
    InVESTwrap::sum_output_lands(landcover_path = landcover_path, output_dir = output_dir,
                                 guild_table = guild_table, raster_types = raster_types,
                                 pest_raster=pest_raster, lg_raster=lg_raster, calc_sum = F,
                                 agg_factor=agg_factor, floral_byspecies=T)

    part3end <- Sys.time()
    
    ####### Fix raster file names
    
    #move pesticide maps into output directory
    pest_tomove <- list.files(output_dir, full.names=T)[grepl(list.files(output_dir), pattern='pesticide')]
    
    for (i in c(1:length(pest_tomove))) {
      file.rename(from=pest_tomove[i], to=paste0(output_dir, "summed_outputs/", basename(pest_tomove[i])))
    }
    
    # #move pesticide maps into output directory
    #read guild table to save foraging distance
    guild <- read.csv(guild_table)
    if (unique(guild$alpha) == 1500) {
      dist <-  "_3km"
    } else if (unique(guild$alpha) == 2500) {
      dist <-  "_5km"
    }
    
    #clean up file names in output directory to simplify
    files <- list.files(paste0(output_dir, "summed_outputs"), full.names=T)
    for (i in 1:length(files)) {
      file.rename(from=files[i], to= gsub(gsub(gsub(gsub(gsub(files[i], pattern="_NA", replacement=""), 
                                                    pattern="apismell", replacement=""), 
                                                    pattern="_distwt", replacement=""), 
                                                    pattern="total_", replacement=""), 
                                                    pattern="_resources", replacement=""))
    }
    
    #change file names to match Azavea conventions and add metadata tag
    files <- list.files(paste0(output_dir, "summed_outputs"), full.names=T)
    
    
    for (i in 1:length(files)) {
      path <- dirname(files[i])
      base <- gsub(gsub(gsub(gsub(gsub(basename(files[i]), pattern=paste0(state, "_"), replacement=""), 
                                  pattern="simple_pesticide", replacement="insecticide"), 
                             pattern=".tif", replacement=paste0(dist, ".tif")), 
                        pattern="mean_", replacement=""),
                   pattern=paste0("_CDL", year), replacement="")
      file.rename(from=files[i], to=paste0(path,"/",state,"_", base))
    }
  }
  
