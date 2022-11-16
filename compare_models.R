library(stringr)
library(ggplot2)
library(overlapping)
library(reshape2)
library(Bchron)

<<<<<<< HEAD
chron_control_types <- read.csv("chroncontrol_types-edited.csv")
=======
chron_control_types <- read.csv("data/chroncontrol_types-edited.csv")
>>>>>>> effea0c267fbe408395932b29c9af8a8f51d57b0

wang_fc = read.csv('Cores_bacon/SiteInfo_fullcore.csv', stringsAsFactors = FALSE)

wang_cores  = list.files('Cores_bacon/Cores_full')
wang_ncores = length(wang_cores)

site_meta = read.csv('bchron_report_v9.0.csv', stringsAsFactors = FALSE)

neo_dat = read.csv("data/pollen_north_america_lc6k2_v8.0.csv")
#neo_table = data.frame(neo_dat$depth, neo_dat$age, neo_dat$dataset_id)

# there are 444 values for datasetid in wang_fc
which(wang_fc$datasetid %in% site_meta$datasetid)
#which(neo_dat$dataset_id %in% site_meta$datasetid)
length(unique(wang_fc$datasetid))


site_meta = read.csv('bchron_report_v9.0.csv', stringsAsFactors = FALSE)
site_meta = site_meta[which(site_meta$success == 1),]

datasetids = unique(site_meta$datasetid)
N_datasetids = length(datasetids)

which(wang_fc$datasetid %in% site_meta$datasetid)
<<<<<<< HEAD



diffs = data.frame(dsid = numeric(0),
                   depths = numeric(0),
                   age_b  = numeric(0),
                   age_w  = numeric(0))

i = 1

pdf('figures/age_depth_compare.pdf', width=10, height=6)
for (i in 1:N_datasetids){#N_datasetids){
  
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
  geochron_quants = data.frame(depth = geochron$depth, geochron_quants)
  
  if (length(idx_dsid) == 0){
    print(paste0('No Bacon age-depth model for dataset id: ', dsid))
    next
  }
  
  geochron_bacon = read.csv(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid], '/', wang_fc$handle[idx_dsid], '.csv'))

  ageSds_bacon = geochron_bacon$error
  calCurves_bacon = rep(NA, nrow(geochron_bacon))
  calCurves_bacon[which(geochron_bacon$cc == 1)] = 'intcal20'
  calCurves_bacon[which(geochron_bacon$cc == 0)] = 'normal'
  
  geochron_bacon_cal <- BchronCalibrate(ages = geochron_bacon$age,
                                  ageSds = ageSds_bacon,
                                  calCurves = calCurves_bacon,
                                  allowOutside = TRUE)
  # # #  we want the weighted means from "calibrated"
  # wmean.date <- function(x) sum(x$ageGrid*x$densities / sum(x$densities))
  # control_young = wmean.date(cal)
  
  geochron_bacon_samples = sampleAges(geochron_bacon_cal)
  geochron_bacon_quants = t(apply(geochron_bacon_samples, 2, quantile, prob=c(0.025, 0.5, 0.975)))
  colnames(geochron_bacon_quants) = c('ylo', 'ymid', 'yhi')
  geochron_bacon_quants = data.frame(depth = geochron_bacon$depth, geochron_bacon_quants)
  
  files = list.files(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid]))                
  idx_file = which(str_sub(files,-8,-1) == 'rout.csv')
  wang_posts = read.csv(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid], '/', files[idx_file]))
  
  wang_depths =  scan(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid], '/', wang_fc$handle[idx_dsid], '_depths.txt'))
  
  # wang_posts = data.frame(depths = wang_depths, wang_posts)
  wang_posts = data.frame(depths = wang_depths, wang_posts)#[,1:100])
  
  wang_posts_long = melt(wang_posts, id.vars = "depths")
  colnames(wang_posts_long) = c('depths', 'iter', 'age')
  
  # wang_posts_long$iter 
  
  # ggplot() +
  #   geom_point(data = wang_posts, aes(x = depths, y = V1)) +
  #   geom_point(data = wang_posts, aes(x = depths, y = V2))
  # 
  # ggplot() +
  #   geom_point(data = wang_posts_long, aes(x = depths, y = age))
  # 
  # ggplot() +
  #   geom_line(data = wang_posts_long, aes(x = depths, y = age))
  # 
  # ggplot() +
  #   geom_line(data = wang_posts_long, aes(x = depths, y = age, group = iter))
  
  
  quantile(wang_posts$V1, c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  wang_quants_row = apply(wang_posts[,2:ncol(wang_posts)], 1, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  wang_quants = data.frame(depths = wang_posts[,1], t(wang_quants_row))
  colnames(wang_quants) = c('depths', 'ylo', 'ymid', 'yhi')
  
  adjustcolor( "red", alpha.f = 0.2)
  
  ggplot() +
    geom_ribbon(data = wang_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#FF000033") +
    geom_line(data = wang_quants, aes(x = depths, y = ymid))
  
  # ### Now do the same with bchron ### #
  
  bchron_posts = read.csv(paste0('Cores/', dsid, '/', dsid, '_bchron_samples.csv'))
  
  if (nrow(bchron_posts)==0){
    print(paste0('No posterior samples for dataset id: ', dsid))
    next
  }
  bchron_posts_long = melt(bchron_posts, id.vars = "depths")
  colnames(bchron_posts_long) = c('depths', 'iter', 'age')
  
  quantile(wang_posts$V1, c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  bchron_quants_row = apply(bchron_posts[,2:ncol(bchron_posts)], 1, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  bchron_quants = data.frame(depths = bchron_posts[,1], t(bchron_quants_row))
  
  colnames(bchron_quants) = c('depths', 'ylo', 'ymid', 'yhi')
  
  adjustcolor( "blue", alpha.f = 0.2)
  
  ggplot() +
    geom_ribbon(data = bchron_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#0000FF33") +
    geom_line(data = bchron_quants, aes(x = depths, y = ymid))
  
  # Now bchron and Bacon on one plot #
  
  colors = c("Bchron" = "blue", "Bacon" = "red")
  
  p <- ggplot() +
    geom_ribbon(data = bchron_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#0000FF33") +
    geom_line(data = bchron_quants, aes(x = depths, y = ymid, color = "Bchron"), size = 1.5) +
    geom_ribbon(data = wang_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#FF000033") +
    geom_line(data = wang_quants, aes(x = depths, y = ymid, color = "Bacon"), size = 1.5) +
    geom_point(data = geochron_quants, aes(x = depth-1, y = ymid), colour='#1F77B4', alpha=0.8) +
    geom_linerange(data = geochron_quants, aes(x = depth-1, ymin = ylo, ymax = yhi), colour='#1F77B4', alpha=0.8, lwd=1) +
    geom_point(data = geochron_bacon_quants, aes(x = depth+1, y = ymid), colour='#D62728', alpha=0.8) +
    geom_linerange(data = geochron_bacon_quants, aes(x = depth+1, ymin = ylo, ymax = yhi), colour='#D62728', alpha=0.8, lwd=1) +
    # geom_point(data = controls, aes(x = depth, y = age)) +
    labs(title = paste0(dsid, '; ',  wang_fc$handle[idx_dsid]), x = "Depths (cm)", y = "Age", color = "Legend") +
    scale_color_manual(values = colors)
  
  print(p)
  
  # ggsave(paste0('figures/age_depth_compare_', dsid, '.png'))
  
  ###
  
  
  wang_mean = data.frame(depths=wang_posts[,'depths'], age_w = rowMeans(wang_posts[,2:ncol(wang_posts)]))
  
  bchron_mean = data.frame(depths=bchron_posts[,'depths'], age_b = rowMeans(bchron_posts[,2:ncol(bchron_posts)]))
  
  age_means = merge(bchron_mean, wang_mean)
  
  diffs = rbind(diffs,
                data.frame(dsid = rep(dsid),
                           age_means))
  
}
dev.off()

# 
# ggplot(data = diffs) +
#   geom_point(aes(x = depths, y = age_b), color = 'blue', alpha = .2) +
#   geom_point(aes(x = depths, y = age_w), color = 'red', alpha = .2)


diffs$age_diff = diffs$age_b - diffs$age_w

ggplot(data=diffs) +
  geom_histogram(aes(x=age_diff, y=..density..), bins=100)

summary(diffs$age_diff)



dsid = datasetids[1]
idx_dsid = which(wang_fc$datasetid == dsid)

files = list.files(paste0('wang/Cores_full/', wang_fc$handle[idx_dsid]))                
idx_file = which(str_sub(files,-8,-1) == 'rout.csv')
wang_posts = read.csv(paste0('wang/Cores_full/', wang_fc$handle[idx_dsid], '/', files[idx_file]))
wang_depths =  scan(paste0('wang/Cores_full/', wang_fc$handle[idx_dsid], '/', wang_fc$handle[idx_dsid], '_depths.txt'))
wang_posts = data.frame(depths = wang_depths, wang_posts)
=======
>>>>>>> effea0c267fbe408395932b29c9af8a8f51d57b0

bchron_posts = read.csv(paste0('Cores/', dsid, '/', dsid, '_bchron_samples.csv'))

n_depths = nrow(bchron_posts)

<<<<<<< HEAD

olap = data.frame(datasetid = numeric(0),
                  depths = numeric(0),
                  overlap = numeric(0))

for (i in 1:n_depths){
  
=======
diffs = data.frame(dsid = numeric(0),
                   depths = numeric(0),
                   age_b  = numeric(0),
                   age_w  = numeric(0),
                   age_n  = numeric(0))

pdf('figures/age_depth_compare.pdf', width=10, height=6)
for (i in 132:N_datasetids){#N_datasetids){
  
  print(i)
  
  dsid = datasetids[i]
  
  idx_dsid = which(wang_fc$datasetid == dsid)
  
  neo_site = neo_dat[which(neo_dat$dataset_id == dsid)]
  
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
  geochron_quants = data.frame(depth = geochron$depth, geochron_quants)
  
  if (length(idx_dsid) == 0){
    print(paste0('No Bacon age-depth model for dataset id: ', dsid))
    next
  }
  
  geochron_bacon = read.csv(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid], '/', wang_fc$handle[idx_dsid], '.csv'))

  ageSds_bacon = geochron_bacon$error
  calCurves_bacon = rep(NA, nrow(geochron_bacon))
  calCurves_bacon[which(geochron_bacon$cc == 1)] = 'intcal20'
  calCurves_bacon[which(geochron_bacon$cc == 0)] = 'normal'
  
  geochron_bacon_cal <- BchronCalibrate(ages = geochron_bacon$age,
                                  ageSds = ageSds_bacon,
                                  calCurves = calCurves_bacon,
                                  allowOutside = TRUE)
  # # #  we want the weighted means from "calibrated"
  # wmean.date <- function(x) sum(x$ageGrid*x$densities / sum(x$densities))
  # control_young = wmean.date(cal)
  
  geochron_bacon_samples = sampleAges(geochron_bacon_cal)
  geochron_bacon_quants = t(apply(geochron_bacon_samples, 2, quantile, prob=c(0.025, 0.5, 0.975)))
  colnames(geochron_bacon_quants) = c('ylo', 'ymid', 'yhi')
  geochron_bacon_quants = data.frame(depth = geochron_bacon$depth, geochron_bacon_quants)
  
  files = list.files(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid]))                
  idx_file = which(str_sub(files,-8,-1) == 'rout.csv')
  wang_posts = read.csv(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid], '/', files[idx_file]))
  
  wang_depths =  scan(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid], '/', wang_fc$handle[idx_dsid], '_depths.txt'))
  
  # wang_posts = data.frame(depths = wang_depths, wang_posts)
  wang_posts = data.frame(depths = wang_depths, wang_posts)#[,1:100])
  
  wang_posts_long = melt(wang_posts, id.vars = "depths")
  colnames(wang_posts_long) = c('depths', 'iter', 'age')
  
  # wang_posts_long$iter 
  
  # ggplot() +
  #   geom_point(data = wang_posts, aes(x = depths, y = V1)) +
  #   geom_point(data = wang_posts, aes(x = depths, y = V2))
  # 
  # ggplot() +
  #   geom_point(data = wang_posts_long, aes(x = depths, y = age))
  # 
  # ggplot() +
  #   geom_line(data = wang_posts_long, aes(x = depths, y = age))
  # 
  # ggplot() +
  #   geom_line(data = wang_posts_long, aes(x = depths, y = age, group = iter))
  
  
  quantile(wang_posts$V1, c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  wang_quants_row = apply(wang_posts[,2:ncol(wang_posts)], 1, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  wang_quants = data.frame(depths = wang_posts[,1], t(wang_quants_row))
  colnames(wang_quants) = c('depths', 'ylo', 'ymid', 'yhi')
  
  adjustcolor( "red", alpha.f = 0.2)
  
  ggplot() +
    geom_ribbon(data = wang_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#FF000033") +
    geom_line(data = wang_quants, aes(x = depths, y = ymid))
  
  # ### Now do the same with bchron ### #
  
  bchron_posts = read.csv(paste0('Cores/', dsid, '/', dsid, '_bchron_samples.csv'))
  
  if (nrow(bchron_posts)==0){
    print(paste0('No posterior samples for dataset id: ', dsid))
    next
  }
  
  
  bchron_posts_long = melt(bchron_posts, id.vars = "depths")
  colnames(bchron_posts_long) = c('depths', 'iter', 'age')
  
  quantile(wang_posts$V1, c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  bchron_quants_row = apply(bchron_posts[,2:ncol(bchron_posts)], 1, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  bchron_quants = data.frame(depths = bchron_posts[,1], t(bchron_quants_row))
  
  colnames(bchron_quants) = c('depths', 'ylo', 'ymid', 'yhi')
  
  adjustcolor( "blue", alpha.f = 0.2)
  
  ggplot() +
    geom_ribbon(data = bchron_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#0000FF33") +
    geom_line(data = bchron_quants, aes(x = depths, y = ymid))
  
  # Now bchron and Bacon on one plot #
  
  colors = c("Bchron" = "blue", "Bacon" = "red")
  
  p <- ggplot() +
    geom_ribbon(data = bchron_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#0000FF33") +
    geom_line(data = bchron_quants, aes(x = depths, y = ymid, color = "Bchron"), size = 1.5) +
    geom_ribbon(data = wang_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#FF000033") +
    geom_line(data = wang_quants, aes(x = depths, y = ymid, color = "Bacon"), size = 1.5) +
    geom_point(data = geochron_quants, aes(x = depth-1, y = ymid), colour='#1F77B4', alpha=0.8) +
    geom_linerange(data = geochron_quants, aes(x = depth-1, ymin = ylo, ymax = yhi), colour='#1F77B4', alpha=0.8, lwd=1) +
    geom_point(data = geochron_bacon_quants, aes(x = depth+1, y = ymid), colour='#D62728', alpha=0.8) +
    geom_linerange(data = geochron_bacon_quants, aes(x = depth+1, ymin = ylo, ymax = yhi), colour='#D62728', alpha=0.8, lwd=1) +
    # geom_point(data = controls, aes(x = depth, y = age)) +
    labs(title = paste0(dsid, '; ',  wang_fc$handle[idx_dsid]), x = "Depths (cm)", y = "Age", color = "Legend") +
    scale_color_manual(values = colors)
  
  print(p)
  
  # ggsave(paste0('figures/age_depth_compare_', dsid, '.png'))
  
  ###
  
  
  wang_mean = data.frame(depths=wang_posts[,'depths'], age_w = rowMeans(wang_posts[,2:ncol(wang_posts)]))
  
  bchron_mean = data.frame(depths=bchron_posts[,'depths'], age_b = rowMeans(bchron_posts[,2:ncol(bchron_posts)]))
  
  neo_age = data.frame(depths=neo_site$depth, age_n = neo_site$age)
  
  agetype = neo_site$agetype[1]
  if (agetype == "Radiocarbon years BP"){
    n_neo = length(neo_age$depths)
    
    
    foo  <- BchronCalibrate(ages = neo_site$age,
                            ageSds = rep(200, n_neo),
                            calCurves = rep('intcal20', n_neo),
                            allowOutside = TRUE)
    
    goo = sampleAges(foo)
    
    neo_age$age_n = colMeans(goo)
    
  }
  
  
  age_means = merge(bchron_mean, wang_mean)
  age_means = merge(age_means, neo_age)
  diffs = rbind(diffs,
                data.frame(dsid = rep(dsid),
                           age_means))
  
}
dev.off()

# Make just two figures comparing models in one panel for paper #






# ggplot(data = diffs) +
#   geom_point(aes(x = depths, y = age_b), color = 'blue', alpha = .2) +
#   geom_point(aes(x = depths, y = age_w), color = 'red', alpha = .2)


diffs$age_diff = diffs$age_b - diffs$age_w

ggplot(data=diffs) +
  geom_histogram(aes(x=age_diff, y=..density..), bins=100)

summary(diffs$age_diff)



dsid = datasetids[1]
idx_dsid = which(wang_fc$datasetid == dsid)

files = list.files(paste0('wang/Cores_full/', wang_fc$handle[idx_dsid]))                
idx_file = which(str_sub(files,-8,-1) == 'rout.csv')
wang_posts = read.csv(paste0('wang/Cores_full/', wang_fc$handle[idx_dsid], '/', files[idx_file]))
wang_depths =  scan(paste0('wang/Cores_full/', wang_fc$handle[idx_dsid], '/', wang_fc$handle[idx_dsid], '_depths.txt'))
wang_posts = data.frame(depths = wang_depths, wang_posts)

bchron_posts = read.csv(paste0('Cores/', dsid, '/', dsid, '_bchron_samples.csv'))

n_depths = nrow(bchron_posts)


olap = data.frame(datasetid = numeric(0),
                  depths = numeric(0),
                  overlap = numeric(0))

for (i in 1:n_depths){
  
>>>>>>> effea0c267fbe408395932b29c9af8a8f51d57b0
  # print(i)
  
  depth = bchron_posts[i,1]
  
  print(depth)
  
  
  bchron_sub = bchron_posts[i,]
  wang_sub   = wang_posts[i,]
  
  if (any(is.na(wang_sub))) {
    next
  }
  
  if (any(is.na(bchron_sub))) {
    next
  }
  
  d_list = list(wang = as.numeric(wang_sub)[-1], bchron = as.numeric(bchron_sub)[-1])
  olap_index = overlap(d_list, plot=TRUE)$OV
  
  olap_site = data.frame(datasetid = dsid, depths = depth, overlap = olap_index)
  
  olap = rbind(olap,
               olap_site)
}








# ###Try overlap function###
# library(reshape2)


olap = data.frame(datasetid = numeric(0),
                  depths = numeric(0),
                  overlap = numeric(0))



bchron_sub = bchron_posts[6,]
wang_sub   = wang_posts[6,]

d_list = list(wang = as.numeric(wang_sub)[-1], bchron = as.numeric(bchron_sub)[-1])
olap_index = overlap(d_list, plot=TRUE)$OV

# 
# for (j in 1:N_datasetids){
#   d_list = list(as.numeric(wang_mean[j,2:ncol(wang_mean)]), as.numeric(bchron_mean[j,2:ncol(bchron_mean)]))
#   
#   olap_index = overlap(d_list, plot=TRUE)$OV
#   
#   olap = rbind(olap,
#                olap_site)
# 
# 
#   }
