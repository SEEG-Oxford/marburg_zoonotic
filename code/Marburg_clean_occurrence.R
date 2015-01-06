# clean marburg occurrence data for modelling

# clear workspace
rm(list = ls())

# load packages
library(seegSDM)
library(snowfall)

# load functions from function file
source('code/R/functions.R')

# load data

# raw occurrence data
dat_marv <- read.csv('Marburg/marburg_raw_141001_v3.csv',
                     stringsAsFactors = FALSE)
# convert all column names to lowercase
colnames(dat_marv) <- tolower(colnames(dat_marv))

#convert human cases to FALSE and animal cases to TRUE
dat_marv[dat_marv == 'Human'] <- FALSE
dat_marv[dat_marv == 'Animal'] <- TRUE

dat_marv$poly<-dat_marv$shape == 'polygon'

# matching shapefile for polygons (all data)
shp_marv <- shapefile('Marburg/GIS_marv/Marburg_V2_polygon.shp')

# ~~~~~~~~~~~~~~
# tidy up the layers
# convert start and end dates into individual columns
dat_marv$start_date <- firstDay(dat_marv$year.start,
                                dat_marv$month.start)

# set up a date vector and populate only the non-missing elements
dat_marv$end_date <- rep(NA, nrow(dat_marv))
class(dat_marv$end_date) <- 'Date'

missing <- is.na(dat_marv$month.end) | is.na(dat_marv$year.end)

dat_marv$end_date[!missing] <- lastDay(dat_marv$year.end[!missing],
                                       dat_marv$month.end[!missing])

dat_marv$poly <- dat_marv$shape == 'polygon'

# keep only the ID, virus, lat, long and start/end dates
dat_marv <- dat_marv[, c('outbreak_id',
                         'virus',
                         'lat',
                         'long',
                         'start_date',
                         'poly',
                         'organism')]
#~~~~~~~~~~~~~~~~
#load covariates
# list the covariates of interest
covars <- c('5km_spatial_summary_rasters/EVI/mean/annual.synoptic.mean.tif',
            '5km_spatial_summary_rasters/EVI/range/annual.synoptic.mean.tif',
            '5km_spatial_summary_rasters/LST_day/mean/annual.synoptic.mean.tif',
            '5km_spatial_summary_rasters/LST_day/range/annual.synoptic.mean.tif',
            '5km_spatial_summary_rasters/LST_night/mean/annual.synoptic.mean.tif',
            '5km_spatial_summary_rasters/LST_night/range/annual.synoptic.mean.tif',
            '5km_spatial_summary_rasters/WorldClim_fourier_precip/prec57a0.tif',
            '5km_spatial_summary_rasters/WorldClim_fourier_precip/prec57a1.tif',
            '5km_spatial_summary_rasters/poverty/accessibility_50k_5km.mean.tif',
            '5km_spatial_summary_rasters/poverty/population_density_5km.mean.tif',
            '5km_spatial_summary_rasters/poverty/PET_1950-2000_5km.mean.tif',
            '5km_spatial_summary_rasters/poverty/srtm_1km_5km.mean.tif',
            '5km_spatial_summary_rasters/TCW/mean/annual.synoptic.mean.tif',
            '5km_spatial_summary_rasters/TCB/mean/annual.synoptic.mean.tif',
            '5km_spatial_summary_rasters/IGBP_Landcover/Majority/2012.majority.class.tif',
            '5km_spatial_summary_rasters/AI_Bare.tif',
            '5km_spatial_summary_rasters/KARST/karst_dist_mask.tif')

# give them readable names
names <- c('EVI_mean',
           'EVI_range',
           'LSTday_mean',
           'LSTday_range',
           'LSTnight_mean',
           'LSTnight_range',
           'precip_mean',
           'precip_annamp',
           'access_50k',
           'afripop',
           'PET_mean',
           'DEM',
           'TCW_mean',
           'TCB_mean',
           'LC_class',
           'arid_areas',
           'Karst')

# loop through opening links to the rasters
covs <- lapply(covars, raster)

# crop the final raster (a global map) to the same extent as the others
covs[[which(names == 'arid_areas')]] <- crop(covs[[which(names == 'arid_areas')]],
                                             covs[[1]])

# invert the arid areas so that 1 is arid
covs[[which(names == 'arid_areas')]] <- 1 - covs[[which(names == 'arid_areas')]]

# stack the rasters
covs <- stack(covs)

# give them nicer names
names(covs) <- names

# set the NA value
NAvalue(covs) <- -9999

# set the minimum value of population to a very small *positive* number
covs$afripop[getValues(covs$afripop) < 0.001] <- 0.001

template<-covs

# set up a cluster
sfInit(cpus = 8, parallel = TRUE)
sfLibrary(seegSDM)

# in parallel, loop through and get coordinates of all the cells covered
pt_list <- sfLapply(1:nrow(shp_marv),
                    function (i, shp_marv, template) {
                      
                      # get the polygon
                      poly <- shp_marv[i, ]
                      
                      # buffer it to make sure it covers the centre of at least one pixel
                      
                      # distance from the corner to the centre of a 5km pixel,
                      # in metres, and then some
                      d <- sqrt(2 * 2500 ^ 2) + 1
                      
                      # convert to decimal degrees (at equator; hacky)
                      d <- d * 10 ^ -5
                      
                      poly <- gBuffer(poly, width = d)
                      
                      # rasterize the layer
                      tmp <- rasterize(poly,
                                       template)
                      
                      # get coordinates of the cells there
                      pts <- xyFromCell(tmp,
                                        which(!is.na(getValues(tmp))))
                      
                      return (pts)
                    },
                    shp_marv,
                    template)

# stop the cluster
sfStop()

# split out the polygons
dat_poly <- dat_marv[dat_marv$poly, ]

# and the points
dat_pt <- dat_marv[!dat_marv$poly, ]

dat_new <- dat_pt

# add a weights column
dat_new$wt <- rep(1, nrow(dat_new))

# loop through each polygon (in the shapefile) and add it to dat
for(i in 1:nrow(shp_marv)) {
  
    
  # get the ID
  ID <- shp_marv@data$outbreakid[i]
  
  # the points
  pts <- pt_list[[i]]
  
  # the number of points
  n <- nrow(pts)
#   
  # if it's larger than 100, pick 100 of them at random
  if (n > 100) {
    pts <- pts[sample(1:n, 100, replace = FALSE), ]
    n <- 100
  }
  
  
  
  if (n == 0) {
    break
  }
  
  # the column info
  info <- dat_poly[dat_poly$outbreak_id == ID, ]
  
  # repeat it n times
  info_mat <- info[rep(1:nrow(info), each = n),] 
  
  # add the weights
  info_mat$wt <- 1 / n
  
  # stick the coordinates in
  info_mat[, c('lat', 'long')] <- pts[, 2:1]
  
  # append it to dat
  dat_new <- rbind(dat_new, info_mat)
}

# if the resulting dataset looks fairly clean, write it to disk
if (!(any(is.na(dat_new)))) {
  
  # output the resulting table
  write.csv(dat_new,
            file = 'Marburg/occurrence.csv',
            row.names = FALSE)
  
} else {
  
  # otherwise throw an error
  stop ('missing values in the dataset!')
  
}
