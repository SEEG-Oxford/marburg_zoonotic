
# ~~~~~~~~~~NICHE MODEL
# remove access, land cover class, tasseled-cap indices and arid areas
# also remove precipitation as it so spatially coarse and population
bad_covs <- c('access_50k',
              'LC_class',
              'TCW_mean',
              'TCB_mean',
              'arid_areas',
              'precip_mean',
              'precip_annamp',
              'afripop')

# get the population raster
pop <- covs[[which(names(covs) == 'afripop')]]

# chuck out the bad covariates
covs <- covs[[which(!(names(covs) %in% bad_covs))]]

# occurrence data
occ <- read.csv('Marburg/occurrence.csv')

###remove animals
drop_animal=FALSE

if (drop_animal) {
  occ <- occ[!occ$organism, ]
}

# ###remove gabon and ROC survey data

occ<- occ[occ$outbreak_id != 26,]
occ<- occ[occ$outbreak_id != 27,]
occ <- occ[occ$outbreak_id != 28, ]
occ<- occ[occ$outbreak_id != 29,]


# generate pseudo-absence data according to the population surface
# (assume perfect detection of human cases, so reported risk only
# a function of population density)

bg <- bgSample(pop,
               n = 10000,
               prob = TRUE,
               replace = TRUE,
               spatial = FALSE)

colnames(bg) <- c('long', 'lat')
bg <- data.frame(bg)

# add an outbreak id to this
bg$outbreak_id <- 0

# combine the occurrence and background records
dat <- rbind(cbind(PA = rep(1, nrow(occ)),
                   occ[, c('long', 'lat', 'outbreak_id')]),
             cbind(PA = rep(0, nrow(bg)),
                   bg))

# get the covariate values
dat_covs <- extract(covs, dat[, 2:3])

# and add them
dat_all <- cbind(dat, dat_covs)

# remove NAs
dat_all <- na.omit(dat_all)

ncpu <- 50
nboot <- ncpu * 10

# create a list with random permutations of dat_all, sampling one occurrence
# from each polygon in each iteration, the bootstrapping as usual.
# This way there's more jitter and it avoids the horrible normalization in
# gbm which otherwise gives some points a weight of ~250 (!).
data_list <- replicate(nboot,
                       subsamplePolys(dat_all,
                                      minimum = c(15, 15)),
                       simplify = FALSE)

# initialize the cluster
sfInit(parallel = TRUE, cpus = ncpu)
sfLibrary(seegSDM)
sfLibrary(PresenceAbsence)
sfLibrary(gbm)
#sfLibrary(snowfall)



#require n.fold = 3 argument as sample size is so small
model_list <- sfLapply(data_list,
                       runBRT,
                       gbm.x = 4:ncol(data_list[[1]]),
                       gbm.y = 1,
                       pred.raster = covs,
                       gbm.coords = 2:3,
                       wt = function(PA) ifelse(PA == 1, 1, sum(PA) / sum(1 - PA)),
                       n.folds = 3)

# get cv statistics in parallel
stat_lis <- sfLapply(model_list, getStats)

# summarise all the ensembles
preds <- stack(lapply(model_list, '[[', 4))

# summarise the predictions in parallel
preds_sry <- combinePreds(preds)

# stop the cluster
sfStop()



######## SAVING OUTPUTS #############
outpath <- 'output/marburg_v3_no_gabon/'


# convert the stats list into a matrix using the do.call function
stats <- do.call("rbind", stat_lis)

# save them
write.csv(stats,
          paste0(outpath, 'stats.csv'))

names(preds_sry) <- c('mean',
                      'median',
                      'lowerCI',
                      'upperCI')

# save the prediction summary
writeRaster(preds_sry,
            file = paste0(outpath,
                          'marburg_v3_no_gabon'),
            format = 'GTiff',
            overwrite = TRUE)

# save the relative influence scores
relinf <- getRelInf(model_list)
write.csv(relinf,
          file = paste0(outpath,
                        'relative_influence.csv'))