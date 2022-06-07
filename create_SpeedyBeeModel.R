#translate spatial functions into correct format for R package
#this allows it to be loaded directly from Github

#adapted from this blog post: https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/

#install.packages("devtools")
library(devtools)
#devtools::install_github("klutometis/roxygen")
library(roxygen2)

#create and name R package (do not run again)
#create("beecoSp")

#add package dependencies
usethis::use_package("raster")
#updated 'description file'


#create documentation
devtools::load_all()
devtools::document()

# #add sample data
# devtools::use_data(ny_landuse, ny_samplesites)


#install package locally
setwd("..") # set wd to one above the current wd to be able to install beecoSp locally
devtools::install("beecoSp")


#try adding from Github
devtools::install_github("land-4-bees/beecoSp")



