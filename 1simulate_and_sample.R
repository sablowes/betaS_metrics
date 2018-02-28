#-----first step for simulation study: simulate and sample communities-----------------
# tidy up 
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)
# load packages
rm(list=ls())
library(mobr)
library(mobsim)
library(tidyverse)
library(cowplot)

S_pool = 100
N = 1000

# want even and less even SADs
cv_abund = list(list(cv_abund=1), 
                list(cv_abund=2),
                list(cv_abund=4))

# generate maps
maps = NULL
maps$agg = lapply(cv_abund, function(x) 
  sim_thomas_community(S_pool, N, 'lnorm', x, fix_s_sim = T))
maps$poi = lapply(cv_abund, function(x) 
  sim_poisson_community(S_pool, N, 'lnorm', x, fix_s_sim = T))


# output maps
setwd('~/Dropbox/1current/conceptual/code/betaS_metrics/')
png('community_maps.png', width = 200, height = 200, units = 'mm', res = 150)

# maps of the four treatments
par(mfrow=c(3,2))
for(i in seq_along(cv_abund)) {
  plot(maps$agg[[i]], axes=F, xlab='', ylab='',
       main=paste('Aggregated (CV = ', cv_abund[[i]]$cv_abund,')'))
  plot(maps$poi[[i]], axes=F, xlab='', ylab='', 
       main=paste('Random (CV = ', cv_abund[[i]]$cv_abund,')'))
}
dev.off()

# sample n_quadrats from the mapped communities: want a perfect sample
n_quadrats = 25

comms = lapply(maps, function(x) 
  lapply(x, function(y) 
    sample_quadrats(y, n_quadrats, plot = F, quadrat_area = 0.04,
                    method = 'grid', delta_x = 0.2, delta_y = 0.2, avoid_overlap = T)))

# illustrate complete sample
png('sample_eg.png', width = 200, height = 200, units = 'mm', res = 150)
sample_quadrats(maps$agg[[1]], n_quadrats, plot = T, quadrat_area = 0.04,
                 method = 'grid', delta_x = 0.2, delta_y = 0.2, avoid_overlap = T)
dev.off()

# aggregate comms data into a community and attributes dataframes
spdat = rbind(dplyr::bind_rows(lapply(comms$agg, function(x) x$spec_dat)),
              dplyr::bind_rows(lapply(comms$poi, function(x) x$spec_dat)))

coords = rbind(dplyr::bind_rows(lapply(comms$agg, function(x) x$xy_dat)),
               dplyr::bind_rows(lapply(comms$poi, function(x) x$xy_dat)))

plot_attr = data.frame(coords, 
                       spatial = rep(c('agg', 'poi'), 
                                     each = n_quadrats * length(cv_abund)),
                       SAD_CV = rep(unlist(cv_abund), 
                                    each= n_quadrats, times = 2),
                       Replicate = rep(c(1:n_quadrats), times = length(cv_abund)*2))

plot_attr$group = paste(plot_attr$spatial, plot_attr$SAD_CV, sep='_')