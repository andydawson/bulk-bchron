library(stringr)
library(ggplot2)
library(overlapping)
library(reshape2)
library(Bchron)

chron_control_types <- read.csv("chroncontrol_types-edited.csv")

wang_fc = read.csv('wang/SiteInfo_fullcore.csv', stringsAsFactors = FALSE)

wang_cores  = list.files('wang/Cores_full')
wang_ncores = length(wang_cores)

site_meta = read.csv('bchron_report_v9.0.csv', stringsAsFactors = FALSE)


# there are 444 values for datasetid in wang_fc
which(wang_fc$datasetid %in% site_meta$datasetid)

length(unique(wang_fc$datasetid))


site_meta = read.csv('bchron_report_v9.0.csv', stringsAsFactors = FALSE)
site_meta = site_meta[which(site_meta$success == 1),]

datasetids = unique(site_meta$datasetid)
N_datasetids = length(datasetids)

which(wang_fc$datasetid %in% site_meta$datasetid)


perc_overlap = function(x.start, x.end, y.start, y.end){
  #  if(x.start == y.start & x.end == y.end){
  #    return(100)
  #  }
  x.len = abs(x.end - x.start)
  # largest start
  max.start = max(c(x.start, y.start))
  min.end = min(c(x.end, y.end))
  overlap = min.end - max.start
  overlap = ifelse(overlap <= 0, 0, overlap)
  perc_overlap = overlap / x.len * 100
  return(perc_overlap)
}

# # bacon_age_posts(d=geo_wang$depth[j], b.depths=geo_wang$depth, out=output, sections=sections)
# bacon_age_posts <- function(d, b.depths, out, sections)
# { 
#   its=out[,1:(ncol(out)-1)]
#   
#   dlo = min(b.depths)
#   dhi = max(b.depths)
#   thick = round_any((dhi - dlo) / sections, 5)
#   
#   elbows <- cbind(its[,1])
#   accs <- its[,2:(ncol(its)-1)]
#   for(i in 2:ncol(accs))
#     elbows <- cbind(elbows, elbows[,ncol(elbows)] + (thick * accs[,i-1]))
#   
#   if (d %in% b.depths){
#     ages <- elbows[,which(b.depths == d)] 
#   } else {
#     maxd <- max(which(b.depths < d))
#     ages <- elbows[,maxd] + ((d-b.depths[maxd]) * accs[,maxd])
#   }
#   ages
# }



diffs = data.frame(dsid = numeric(0),
                   depths = numeric(0),
                   olap_bchron = numeric(0),
                   olap_bacon  = numeric(0))


for (i in 1:20){#N_datasetids){
  
  print(i)
  
  dsid = datasetids[i]
  
  idx_dsid = which(wang_fc$datasetid == dsid)
  
  geochron = read.csv(paste0('Cores/', dsid, '/', dsid, '_prepared.csv'))
  geochron = geochron[which(geochron$keep == 1),]
  
  ageSds = geochron$error
  calCurves = rep(NA, nrow(geochron))
  calCurves[which(geochron$cc == 1)] = 'intcal20'
  calCurves[which(geochron$cc == 0)] = 'normal'
  
  geochron_cal <- BchronCalibrate(ages = geochron$age,
                                  ageSds = ageSds,
                                  calCurves = calCurves,
                                  allowOutside = TRUE)
  # # #  we want the weighted means from "calibrated"
  # wmean.date <- function(x) sum(x$ageGrid*x$densities / sum(x$densities))
  # control_young = wmean.date(cal)
  
  geochron_samples = sampleAges(geochron_cal)
  geochron_quants = t(apply(geochron_samples, 2, quantile, prob=c(0.025, 0.5, 0.975)))
  colnames(geochron_quants) = c('ylo', 'ymid', 'yhi')
  geo_quants = data.frame(depth = geochron$depth, geochron_quants)
  
  if (length(idx_dsid) == 0){
    print(paste0('No Bacon age-depth model for dataset id: ', dsid))
    next
  }
  
  files = list.files(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid]))                
  idx_file = which(str_sub(files,-8,-1) == 'rout.csv')
  wang_posts = read.csv(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid], '/', files[idx_file]))
  
  
  # geo_wang =  read.csv(paste0('wang/Cores_full/', wang_fc$handle[idx_dsid], '/', wang_fc$handle[idx_dsid], '.csv'))
  
  bacon_posts = read.csv(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid], '/', wang_fc$handle[idx_dsid], '_geo_samples.csv'))
  
  bacon_depths = bacon_posts[,1]
  
  bacon_posts = bacon_posts[,2:ncol(bacon_posts)]
  bacon_posts = bacon_posts[,sample(ncol(bacon_posts), 1000)]
  # 
  # sections = as.numeric(str_split(files[idx_out], "[._]")[[1]][2])
  

  bacon_posts = data.frame(depths=bacon_depths, bacon_posts)
  bacon_posts_long = melt(bacon_posts, id.vars = "depths")
  colnames(bacon_posts_long) = c('depths', 'iter', 'age')
  
  
  # foo = apply(bacon_posts[, -1], 1, mean)
  # bar = apply(bacon_posts[, -1], 1, sd)
  # plot(foo, bar)
  
  bacon_quants_row = apply(bacon_posts[,2:ncol(bacon_posts)], 1, function(x) quantile(x, c(0.4, 0.5, 0.6), na.rm = TRUE))
  bacon_quants = data.frame(depths = bacon_posts[,1], t(bacon_quants_row))
  colnames(bacon_quants) = c('depths', 'ylo', 'ymid', 'yhi')
  
  # ### Now do the same with bchron ### #
  
  bchron_posts = read.csv(paste0('Cores/', dsid, '/', dsid, '_bchron_geo_samples.csv'))
  bchron_posts = bchron_posts[,2:ncol(bchron_posts)]
  
  if (nrow(bchron_posts)==0){
    print(paste0('No posterior samples for dataset id: ', dsid))
    next
  }
  bchron_posts_long = melt(bchron_posts, id.vars = "depths")
  colnames(bchron_posts_long) = c('depths', 'iter', 'age')
  
  bchron_quants_row = apply(bchron_posts[,2:ncol(bchron_posts)], 1, function(x) quantile(x, c(0.4, 0.5, 0.6), na.rm = TRUE))
  bchron_quants = data.frame(depths = bchron_posts[,1], t(bchron_quants_row))
  colnames(bchron_quants) = c('depths', 'ylo', 'ymid', 'yhi')
  
  
  ngeo = nrow(geo_quants)
  
  for (n in 1:ngeo){
    
    geo_depth = geo_quants[n,1]
    geo_lo = geo_quants[n,2]
    geo_hi = geo_quants[n,4]
    
    if (!is.na(match(geo_depth, bchron_quants$depths))){
      bchron_lo = bchron_quants[match(geo_depth, bchron_quants$depths), 2]
      bchron_hi = bchron_quants[match(geo_depth, bchron_quants$depths), 4]
      # olap_bchron = perc_overlap(bchron_lo, bchron_hi, geo_lo, geo_hi)
      olap_bchron = perc_overlap(geo_lo, geo_hi, bchron_lo, bchron_hi)
      
    } else {
      olap_bchron = NA
    }
    
    if (!is.na(match(geo_depth, bacon_quants$depths))){
      bacon_lo = bacon_quants[match(geo_depth, bacon_quants$depths), 2]
      bacon_hi = bacon_quants[match(geo_depth, bacon_quants$depths), 4]
      # olap_bacon  = perc_overlap(bacon_lo, bacon_hi, geo_lo, geo_hi)
      olap_bacon  = perc_overlap(geo_lo, geo_hi, bacon_lo, bacon_hi)
    } else {
      olap_bacon = NA
    }
    
    # fraction of individually calibrated date range that falls in age-depth model pred ages    
    
    diffs = rbind(diffs,
                  data.frame(dsid = dsid,
                             depth = geo_depth,
                             olap_bchron = olap_bchron,
                             olap_bacon = olap_bacon))
    
  }
  
}


diffs_both = diffs[which((!is.na(diffs$olap_bacon)) & (!is.na(diffs$olap_bchron))),]

plot(diffs_both$olap_bacon, diffs_both$olap_bchron)

ggplot(data=diffs_both, aes(x=olap_bacon, y=olap_bchron)) +
  geom_point() +
  geom_smooth(method='lm', formula = y~x)

diffs_melt = melt(diffs_both, id.vars = c('dsid', 'depth'))

ggplot(data=diffs_melt, aes(x=variable, y=value)) + 
  geom_violin()

ggplot(data=diffs_melt, aes(x=depth, y=value)) + 
  geom_point(aes(colour=variable))

ggplot(data=diffs_melt, aes(x=value)) + 
  geom_histogram() + facet_grid(variable~.)

# 
# ggplot(data = diffs) +
#   geom_point(aes(x = depths, y = age_b), color = 'blue', alpha = .2) +
#   geom_point(aes(x = depths, y = age_w), color = 'red', alpha = .2)


# diffs$age_diff = diffs$age_b - diffs$age_w
# 
# ggplot(data=diffs) +
#   geom_histogram(aes(x=age_diff, y=..density..), bins=100)
# 
# summary(diffs$age_diff)
# 
# 
# 
# dsid = datasetids[1]
# idx_dsid = which(wang_fc$datasetid == dsid)
# 
# files = list.files(paste0('wang/Cores_full/', wang_fc$handle[idx_dsid]))                
# idx_file = which(str_sub(files,-8,-1) == 'rout.csv')
# wang_posts = read.csv(paste0('wang/Cores_full/', wang_fc$handle[idx_dsid], '/', files[idx_file]))
# wang_depths =  scan(paste0('wang/Cores_full/', wang_fc$handle[idx_dsid], '/', wang_fc$handle[idx_dsid], '_depths.txt'))
# wang_posts = data.frame(depths = wang_depths, wang_posts)
# 
# bchron_posts = read.csv(paste0('Cores/', dsid, '/', dsid, '_bchron_samples.csv'))
# 
# n_depths = nrow(bchron_posts)
# 
# 
# olap = data.frame(datasetid = numeric(0),
#                   depths = numeric(0),
#                   overlap = numeric(0))
# 
# for (i in 1:n_depths){
#   
#   # print(i)
#   
#   depth = bchron_posts[i,1]
#   
#   print(depth)
#   
#   
#   bchron_sub = bchron_posts[i,]
#   wang_sub   = wang_posts[i,]
#   
#   if (any(is.na(wang_sub))) {
#     next
#   }
#   
#   if (any(is.na(bchron_sub))) {
#     next
#   }
#   
#   d_list = list(wang = as.numeric(wang_sub)[-1], bchron = as.numeric(bchron_sub)[-1])
#   olap_index = overlap(d_list, plot=TRUE)$OV
#   
#   olap_site = data.frame(datasetid = dsid, depths = depth, overlap = olap_index)
#   
#   olap = rbind(olap,
#                olap_site)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# # ###Try overlap function###
# # library(reshape2)
# 
# 
# olap = data.frame(datasetid = numeric(0),
#                   depths = numeric(0),
#                   overlap = numeric(0))
# 
# 
# 
# bchron_sub = bchron_posts[6,]
# wang_sub   = wang_posts[6,]
# 
# d_list = list(wang = as.numeric(wang_sub)[-1], bchron = as.numeric(bchron_sub)[-1])
# olap_index = overlap(d_list, plot=TRUE)$OV
# 
# # 
# # for (j in 1:N_datasetids){
# #   d_list = list(as.numeric(wang_mean[j,2:ncol(wang_mean)]), as.numeric(bchron_mean[j,2:ncol(bchron_mean)]))
# #   
# #   olap_index = overlap(d_list, plot=TRUE)$OV
# #   
# #   olap = rbind(olap,
# #                olap_site)
# # 
# # 
# #   }
