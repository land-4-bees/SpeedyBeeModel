% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/insecticide_index.R
\name{insecticide_index}
\alias{insecticide_index}
\title{Calculate bee landscape indices}
\usage{
insecticide_index(
  output_dir,
  landcover_path,
  pesticide_path,
  forage_range = NA,
  guild_table = NA,
  ins_method = "mean",
  agg_factor = NA,
  normalize = F,
  useW = F,
  check_pesttable = T,
  rastertag = "insecticide",
  verbose = T
)
}
\arguments{
\item{output_dir}{Path to directory for output files.}

\item{landcover_path}{Path to land cover raster, including base file name}

\item{pesticide_path}{Path to pesticide load by land use table}

\item{forage_range}{Foraging range (in m) to use for distance weighting insecticide scores surrounding focal cell.}

\item{guild_table}{Bee community to use to model foraging activity. Includes foraging range and relative abundnace of each species.}

\item{ins_method}{Insecticide toxic load measure to use. Specify 'oral', 'contact' or 'mean' (mean of 'oral' and 'contact').}

\item{agg_factor}{Aggregation factor for large rasters (use 4 to convert 50m to 120m resolution CDL)}

\item{normalize}{Normalize values by the number of cells within each moving window?}

\item{rastertag}{Text string to include in name of output raster, optional}

\item{verbose}{Include more log messages from model run?}
}
\description{
New implemenatation of landsacpe insecticide toxic load index (Douglas et al.). Uses FFT convolution to improve runtime.
}
\keyword{index}
\keyword{insecticide}
