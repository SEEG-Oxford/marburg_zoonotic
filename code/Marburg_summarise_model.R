# extract statistics for the MARV model 

# clear workspace
rm(list = ls())

# source functions files
source('code/R/functions.R')

paths <- c('output/marburg_v3_all/',
           'output/marburg_v3_human/',
           'output/marburg_v3_no_gabon/')

names <- c('All MARV',
           'Human MARV',
           'No Gabon MARV')

sry <- lapply(paths, summarizeStats)

sry <- do.call(cbind, sry)
sry <- data.frame(sry)
names(sry) <- names

print(sry)

write.csv(sry,
          file = 'output/model_v3_summary/sry.csv')
