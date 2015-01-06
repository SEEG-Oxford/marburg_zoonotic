# Threshold the MARV risk map ensembles at different levels

# clear workspace
rm(list = ls())

# load packages
library(raster)
library(snowfall)

# load functions from function files
source('code/R/functions.R')

# load data

# get mean EBOV prediction raster
EBOV_risk <- brick('output/marburg_v3_human/marburg_v3_human.tif')[[1]]


# load the covariate rasters
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

# pull out the population layer
pop <- covs[[which(names(covs) == 'afripop')]]

# remove the other covariates
rm(covs)

# load an urban raster
urb <- raster('5km_spatial_summary_rasters/grumpALPHA/upr_u_5km.flt')

# and a peri-urban raster
peri <- raster('5km_spatial_summary_rasters/grumpALPHA/upr_p_5km.flt')

# crop these to match pop
urb <- crop(urb, pop)
peri <- crop(peri, pop)

# combine them to get a rural layer
tmp <- urb + peri
rural <- tmp == 0

# remove these other layers
rm(list = c('urb', 'peri', 'tmp'))

# create rural and urban population rasters
pop_urban <- pop * !rural
pop_rural <- pop * rural

# combine these
pop_all <- brick(pop,
                 pop_urban,
                 pop_rural)

names(pop_all) <- c('all', 'urban_periurban', 'rural')

# if an admin0 file doesn't already exist, make one
if (!(file.exists('tmp/5km/ad0.grd') & file.exists('tmp/5km/GAUL_lookup.csv'))) {
  
  # load an admin 0 layer shapefile
  ad0_shp <- shapefile('5km_spatial_summary_rasters/admin_units/admin2013_0.shp')
  
  # rasterize it to the same resolution and extent as the population
  # and prediction layers
  ad0_raster <- rasterize(ad0_shp,
                          pop,
                          field = 'GAUL_CODE')
  
  # save the lookup table for country names and GAUL codes
  GAUL <- ad0_shp@data[, c('NAME', 'GAUL_CODE')]
  
  # save these
  writeRaster(ad0_raster,
              file = 'tmp/5km/ad0.grd',
              overwrite = TRUE)
  
  write.csv(GAUL,
            file = 'tmp/5km/GAUL_lookup.csv',
            row.names = FALSE)
  
} else {
  
  # otherwise load them
  ad0_raster <- raster('tmp/5km/ad0.grd')
  GAUL <- read.csv('tmp/5km/GAUL_lookup.csv')
  
}

# occurrence data
occ <- read.csv('Marburg/occurrence.csv')

###remove animals
drop_animal=TRUE

if (drop_animal) {
  occ <- occ[!occ$organism, ]
}

# 
# # ###remove gabon and ROC survey data
# 
# occ<- occ[occ$outbreak_id != 26,]
# occ<- occ[occ$outbreak_id != 27,]
# occ <- occ[occ$outbreak_id != 28, ]
# occ<- occ[occ$outbreak_id != 29,]

## retrieve min value for polygon
shp_marv <- shapefile('Marburg/GIS_marv/Marburg_V2_polygon.shp')
mean_poly<-data.frame(extract(EBOV_risk, shp_marv, mean))
names(mean_poly)<-c('prob')

#new file for probabilities
prob_occ<-data.frame(unique(occ$outbreak_id), rep(0,length(unique(occ$outbreak_id))))
names(prob_occ)<-c('outbreak_id', 'probability')

# define outbreak_id of shp file
ID <- shp_marv@data$outbreakid


for (i in 1:nrow(prob_occ)){
  extract_ID<-prob_occ$outbreak_id[i]
  if(extract_ID %in% ID){
    shape_id<-which(ID==extract_ID)
    prob_occ$probability[i]<-mean_poly[shape_id,]
    
  }
  else {
   point_id<-which(occ$outbreak_id==extract_ID)
    prob_occ$probability[i]<-extract(EBOV_risk, occ[point_id, c('long','lat')]) 
  }
  
}

prob_occ<-na.omit(prob_occ)

proportion<-1

thresh<-quantile(prob_occ,
                 1-proportion,
                 na.rm=TRUE)

EBOV_pres<-EBOV_risk > thresh


# define the core set of countries
core_set <- c("Democratic Republic of the Congo",
              "Kenya",
              "Zimbabwe",
              "Angola",
              "Uganda")
              
# multiply this by population to get populations at risk in urban/peri_urban and rural areas
EBOV_PAR <- EBOV_pres * pop_all

# give the layers their names
names(EBOV_PAR) <- c('all', 'urban_periurban', 'rural')

# get the number of pixels predicted to be at risk in each country
risk_sums <- zonal(EBOV_pres,
                   ad0_raster,
                   fun = 'sum',
                   na.rm = TRUE)

# and the corresponding population at risk
PAR_sums <- zonal(EBOV_PAR,
                  ad0_raster,
                  fun = 'sum',
                  na.rm = TRUE)

# combine the two columns dataframes, if the first columns match
if (all.equal(risk_sums[, 1], PAR_sums[, 1])) {
  dat <- data.frame(country = risk_sums[, 1],
                    risk_pixels = risk_sums[, -1],
                    round(PAR_sums[, -1]))
} else {
  
  stop ('GAUL codes do not match!')
  
}

# lookup the country names
idx <- match(dat$country, GAUL[, 2])

# replace the GAUL zones with the country names
dat$country <- GAUL[idx, 1]

# subset this to only positives
dat <- dat[dat$risk_pixels > 0, ]

# add a column denoting whether they're in the core block
# start with false
dat$has_cases <- FALSE

# find those that are in the core block
in_core <- match(core_set, dat$country)

# and set these to true
dat$has_cases[in_core] <- TRUE

# sort by PAR, then whether cases are present
dat <- dat[order(dat$all, decreasing = TRUE), ]
dat <- dat[order(dat$has_cases, decreasing = TRUE), ]

# then save
write.csv(dat,
          file = 'output/model_v3_summary/country_risk_all.csv',
          row.names = FALSE)
