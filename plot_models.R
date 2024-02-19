library(stringr)
library(ggplot2)
library(overlapping)
library(reshape2)
library(Bchron)
library(cowplot)
library(mgcv)
library(dplyr)



vers = 1.0

radio = read.csv('data/radiocarbon-dates-errors.csv')
mod_radio <- gam(error ~ s(age, k=15), data=radio, method='REML', family=Gamma(link="identity"))

chron_control_types <- read.csv("data/chroncontrol_types-edited.csv")

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

geochron_neo = read.csv('data/chroncontrol_summary_v2.csv')

## This doesn't work! But we need to clean up  geochron_neo more ##
#geochron_neo$age = na.omit(geochron_neo$age)

n_neo = length(geochron_neo$depth)
geochron_neo$error = (geochron_neo$limitolder - geochron_neo$limityounger)

diffs = data.frame(dsid = numeric(0),
                   depths = numeric(0),
                   age_mean_b  = numeric(0),
                   age_sd_b  = numeric(0),
                   age_mean_w  = numeric(0),
                   age_sd_b  = numeric(0),
                   age_mean_n  = numeric(0),
                   age_sd_b  = numeric(0))

# pdf('figures/age_depth_compare.pdf', width=10, height=6)
for (i  in 1:1){#N_datasetids){
  
  print(i)
  
  dsid = datasetids[i]
  
  if (i==20295){
    next
  }
  
  idx_dsid = which(wang_fc$datasetid == dsid)
  
  neo_site = neo_dat[which(neo_dat$dataset_id == dsid),]
  
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
  
  ##
  ## Bacon posterior samples
  ##
  
  # posteriors: pollen sample depths
  files = list.files(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid]))                
  idx_file = which(str_sub(files,-8,-1) == 'rout.csv')
  wang_posts = read.csv(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid], '/', files[idx_file]))
  wang_depths =  scan(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid], '/', wang_fc$handle[idx_dsid], '_depths.txt'))
  # wang_posts = data.frame(depths = wang_depths, wang_posts)
  wang_posts = data.frame(depths = wang_depths, wang_posts)#[,1:100])
  wang_posts_long = melt(wang_posts, id.vars = "depths")
  colnames(wang_posts_long) = c('depths', 'iter', 'age')
  
  wang_quants_row = apply(wang_posts[,2:ncol(wang_posts)], 1, 
                          function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  wang_quants = data.frame(depths = wang_posts[,1], t(wang_quants_row))
  colnames(wang_quants) = c('depths', 'ylo', 'ymid', 'yhi')
  
  
  # posteriors: geochron depths
  wang_geo_posts = read.csv(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid], "/", 
                          wang_fc$handle[idx_dsid], "_span_samples.csv"))
  
  wang_geo_posts_long = melt(wang_geo_posts, id.vars = "depths")
  colnames(wang_geo_posts_long) = c('depths', 'iter', 'age')
  
  wang_geo_quants_row = apply(wang_geo_posts[,2:ncol(wang_geo_posts)], 1, 
                          function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  wang_geo_quants = data.frame(depths = wang_geo_posts[,1], t(wang_geo_quants_row))
  colnames(wang_geo_quants) = c('depths', 'ylo', 'ymid', 'yhi')
  
  # posteriors: span across all depths
  wang_span_posts = read.csv(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid], "/", 
                           wang_fc$handle[idx_dsid], "_span_samples.csv"))
  
  wang_span_posts_long = melt(wang_span_posts, id.vars = "depths")
  colnames(wang_span_posts_long) = c('depths', 'iter', 'age')
  
  wang_span_quants_row = apply(wang_span_posts[,2:ncol(wang_span_posts)], 1, 
                              function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  wang_span_quants = data.frame(depths = wang_span_posts[,1], t(wang_span_quants_row))
  colnames(wang_span_quants) = c('depths', 'ylo', 'ymid', 'yhi')

  

  adjustcolor( "red", alpha.f = 0.2)
  
  ggplot() +
    geom_ribbon(data = wang_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#FF000033") +
    geom_line(data = wang_quants, aes(x = depths, y = ymid))
  
  ##
  ## Bchron posterior samples
  ##
  
  
  # posteriors: pollen sample depths
  bchron_posts = read.csv(paste0('Cores/', dsid, '/', dsid, '_bchron_samples.csv'))
  
  if (nrow(bchron_posts)==0){
    print(paste0('No posterior samples for dataset id: ', dsid))
    next
  }
  
  bchron_posts_long = melt(bchron_posts, id.vars = "depths")
  colnames(bchron_posts_long) = c('depths', 'iter', 'age')
  
  bchron_quants_row = apply(bchron_posts[,2:ncol(bchron_posts)], 1, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  bchron_quants = data.frame(depths = bchron_posts[,1], t(bchron_quants_row))
  colnames(bchron_quants) = c('depths', 'ylo', 'ymid', 'yhi')
  
  # posteriors: geochron depths
  bchron_geo_posts = read.csv(paste0('Cores/', dsid, '/', dsid, '_bchron_geo_samples.csv'))
  bchron_geo_posts = bchron_geo_posts[,-1]
  
  bchron_geo_posts_long = melt(bchron_geo_posts, id.vars = c("depths"))
  colnames(bchron_geo_posts_long) = c('depths', 'iter', 'age')
  
  bchron_geo_quants_row = apply(bchron_geo_posts[,2:ncol(bchron_geo_posts)], 1, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  bchron_geo_quants = data.frame(depths = bchron_geo_posts[,1], t(bchron_geo_quants_row))
  colnames(bchron_geo_quants) = c('depths', 'ylo', 'ymid', 'yhi')
  
  # posteriors: span across all depths
  bchron_span_posts = read.csv(paste0('Cores/', dsid, '/', dsid, '_bchron_span_samples.csv'))
  
  bchron_span_posts_long = melt(bchron_span_posts, id.vars = c("depths"))
  colnames(bchron_span_posts_long) = c('depths', 'iter', 'age')
  
  bchron_span_quants_row = apply(bchron_span_posts[,2:ncol(bchron_span_posts)], 1, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  bchron_span_quants = data.frame(depths = bchron_span_posts[,1], t(bchron_span_quants_row))
  colnames(bchron_span_quants) = c('depths', 'ylo', 'ymid', 'yhi')
  
  adjustcolor( "blue", alpha.f = 0.2)
  
  # geochron_neo_site = geochron_neo[which(geochron_neo$datasetid == dsid),]
  # ageSds_neo = geochron_neo_site$error
  # calCurves_neo = rep(NA, nrow(geochron_neo_site))
  # 
  # geochron_neo_cal  <- BchronCalibrate(ages = geochron_neo_site$age,
  #                                      ageSds = ageSds_neo,
  #                                      calCurves = rep('intcal20', n_neo),
  #                                      allowOutside = TRUE)
  # 
  # geochron_neo_samples = sampleAges(geochron_neo_cal)
  # geochron_neo_quants = t(apply(geochron_neo_samples,2,quantile, prob=c(0.025, 0.5, 0.975)))
  # colnames(geochron_neo_quants) = c('ylo', 'ymid', 'yhi')
  # geochron_neo_quants = data.frame(depth = geochron_neo_site$depth, geochron_neo_quants)
  
  
  neo_site = neo_dat[which(neo_dat$dataset_id == dsid),]
  neo_site = neo_site[which(!is.na(neo_site$age)),]
  
  if (nrow(neo_site) == 0){
    
    neo_quants = data.frame(depth = numeric(0),
                            ylo   = numeric(0),
                            ymid  = numeric(0),
                            yhi   = numeric(0))
    
    neo_mean = data.frame(depths = numeric(0),
                          age_mean_n = numeric(0),
                          age_sd_n = numeric(0))
    
    # neo_quants = data.frame(depth = NA,
    #                         ylo   = NA,
    #                         ymid  = NA,
    #                         yhi   = NA)
    
    # neo_mean = data.frame(depth = NA,
    #                       age_n = NA)
    
  } else if (neo_site$agetype[1] == "Radiocarbon years BP") {
    
    n_neo_site = nrow(neo_site)

    # ageSds_neo_site = rep(500, n_neo_site)
    ageSds_neo_site = rep(NA, n_neo_site)

    for (j in 1:n_neo_site){
      ageSds_neo_site[j] = round(mean(predict.gam(mod_radio, new_age=data.frame(age=neo_site$age[j]), type='response', se.fit=TRUE)$fit))

    }
    neo_cal  <- BchronCalibrate(ages = neo_site$age,
                                ageSds = ageSds_neo_site,
                                calCurves = rep('intcal20', n_neo_site),
                                allowOutside = TRUE)

    neo_samples = sampleAges(neo_cal)
    neo_quants = t(apply(neo_samples,2,quantile, prob=c(0.025, 0.5, 0.975)))
    colnames(neo_quants) = c('ylo', 'ymid', 'yhi')
    neo_quants = data.frame(depth = neo_site$depth, neo_quants)

    neo_mean = data.frame(depths = neo_site$depth,
                          age_mean_n = apply(neo_samples,2, mean, na.rm=TRUE),
                          age_sd_n = apply(neo_samples,2, sd, na.rm=TRUE))
    
    
    geochron_neo = read.csv(paste0('Cores/', dsid, '/', dsid, '.csv'))
    geochron_neo = geochron_neo[order(geochron_neo$depth),]

    ageSds_neo_site = rep(NA, nrow(geochron_neo))
    ageSds_neo_site = abs((geochron_neo$limitolder - geochron_neo$limityounger) / 4)
    # for (j in 1:n_neo_site){
    #   ageSds_neo_site[j] = round(mean(predict.gam(mod_radio, 
    #                                               new_age=data.frame(age=geochron_neo$age[j]), 
    #                                               type='response', 
    #                                               se.fit=TRUE)$fit))
    # 
    # }
    
    if (any(is.na(ageSds_neo_site))){
      ageSds_neo_site[which(is.na(ageSds_neo_site))] = 10
    }
    
    # geochron$error[match(geochron_neo$chroncontrolid, geochron$chroncontrolid)]
    
    calCurves_neo = rep(NA, nrow(geochron_neo))
    geochron_neo$cc = geochron$cc[match(geochron_neo$chroncontrolid, geochron$chroncontrolid)]
    
    if (dsid %in% c(1008, 1826)){
      geochron_neo$cc[which(geochron_neo$type == 'Extrapolated')] = 1
    } else if (dsid==1732){
      geochron_neo$cc[which(geochron_neo$type == 'Interpolated')] = 1
    } else if (dsid==2270){
      geochron_neo$cc[which(geochron_neo$type == 'Interpolated')] = 0
    } else if (any(is.na(geochron_neo$cc))) {
      geochron_neo$cc[which(geochron_neo$type == 'Biostratigraphic, pollen')] = 0
      geochron_neo$cc[which(geochron_neo$type == 'Sediment stratigraphic')] = 0
      geochron_neo$cc[which(geochron_neo$type == 'Introduced Pinus rise')] = 0
      geochron_neo$cc[which(geochron_neo$type == 'Guess')] = 0
    }
    
    calCurves_neo[which(geochron_neo$cc == 1)] = 'intcal20'
    calCurves_neo[which(geochron_neo$cc == 0)] = 'normal'
    
    geochron_cal_neo <- BchronCalibrate(ages = geochron_neo$age,
                                    ageSds = ageSds_neo_site,
                                    calCurves = calCurves_neo,
                                    allowOutside = TRUE)
    # # #  we want the weighted means from "calibrated"
    # wmean.date <- function(x) sum(x$ageGrid*x$densities / sum(x$densities))
    # control_young = wmean.date(cal)
    
    geochron_neo_samples = sampleAges(geochron_cal_neo)
    geochron_neo_quants = t(apply(geochron_neo_samples, 2, quantile, prob=c(0.025, 0.5, 0.975)))
    colnames(geochron_neo_quants) = c('ylo', 'ymid', 'yhi')
    geochron_neo_quants = data.frame(depth = geochron_neo$depth, geochron_neo_quants)
    
  } else {
    neo_quants = data.frame(depth = neo_site$depth,
                            ylo = neo_site$age,
                            ymid = neo_site$age,
                            yhi = neo_site$age)
    
    neo_mean = data.frame(depths = neo_site$depth,
                          age_mean_n = neo_site$age,
                          age_sd_n = 0)
    
    geochron_neo = read.csv(paste0('Cores/', dsid, '/', dsid, '.csv'))
    geochron_neo = geochron_neo[order(geochron_neo$depth),]
    
    ageSds_neo_site = rep(NA, nrow(geochron_neo))
    ageSds_neo_site = (geochron_neo$limitolder - geochron_neo$limityounger) / 4
    if (any(is.na(ageSds_neo_site))){
      ageSds_neo_site[which(is.na(ageSds_neo_site))] = 10
    }
    
    # for (j in 1:n_neo_site){
    #   ageSds_neo_site[j] = round(mean(predict.gam(mod_radio, 
    #                                               new_age=data.frame(age=geochron_neo$age[j]), 
    #                                               type='response', 
    #                                               se.fit=TRUE)$fit))
    # 
    # }
    
    # geochron$error[match(geochron_neo$chroncontrolid, geochron$chroncontrolid)]
    
    calCurves_neo = rep(NA, nrow(geochron_neo))
    geochron_neo$cc = geochron$cc[match(geochron_neo$chroncontrolid, geochron$chroncontrolid)]
    
    if (dsid==1008){
      geochron_neo$cc[which(geochron_neo$type == 'Extrapolated')] = 1
    } else if (any(is.na(geochron_neo$cc))) {
      geochron_neo$cc[which(geochron_neo$type == 'Biostratigraphic, pollen')] = 0
      geochron_neo$cc[which(geochron_neo$type == 'Sediment stratigraphic')] = 0
      geochron_neo$cc[which(geochron_neo$type == 'IDW-2d Picea decline')] = 0
      geochron_neo$cc[which(geochron_neo$type == 'Tsuga decline')] = 0
    }
    
    calCurves_neo[which(geochron_neo$cc == 1)] = 'intcal20'
    calCurves_neo[which(geochron_neo$cc == 0)] = 'normal'
    
    geochron_cal_neo <- BchronCalibrate(ages = geochron_neo$age,
                                        ageSds = ageSds_neo_site,
                                        calCurves = calCurves_neo,
                                        allowOutside = TRUE)
    # # #  we want the weighted means from "calibrated"
    # wmean.date <- function(x) sum(x$ageGrid*x$densities / sum(x$densities))
    # control_young = wmean.date(cal)
    
    geochron_neo_samples = sampleAges(geochron_cal_neo)
    geochron_neo_quants = t(apply(geochron_neo_samples, 2, quantile, prob=c(0.025, 0.5, 0.975)))
    colnames(geochron_neo_quants) = c('ylo', 'ymid', 'yhi')
    geochron_neo_quants = data.frame(depth = geochron_neo$depth, geochron_neo_quants)
    
    
  }
  
  # Now bchron, Bacon, and Neotoma? on one plot #
  
  colors = c("Bchron" = "blue", "Bacon" = "red", "Neotoma" = "orange")
  
 add_olap = readRDS('olap_dots.RDS')
  
  p <- ggplot() +
    geom_ribbon(data = bchron_span_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#0000FF33") +
    geom_line(data = bchron_span_quants, aes(x = depths, y = ymid, color = "Bchron"), size = 1.5) +
    geom_ribbon(data = wang_span_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#FF000033") +
    geom_line(data = wang_span_quants, aes(x = depths, y = ymid, color = "Bacon"), size = 1.5) + 
    geom_ribbon(data = neo_quants, aes(x = depth, ymin = ylo, ymax = yhi), fill = '#FFA500AA') +
    geom_line(data = neo_quants, aes(x = depth, y = ymid, colour = "Neotoma")) +
    geom_point(data = geochron_quants, aes(x = depth-1, y = ymid), colour='#1F77B4', alpha=0.8) +
    geom_linerange(data = geochron_quants, aes(x = depth-1, ymin = ylo, ymax = yhi), colour='#1F77B4', alpha=0.8, lwd=1) +
    geom_point(data = geochron_bacon_quants, aes(x = depth+1, y = ymid), colour='#D62728', alpha=0.8) +
    geom_linerange(data = geochron_bacon_quants, aes(x = depth+1, ymin = ylo, ymax = yhi), colour='#D62728', alpha=0.8, lwd=1) +
    # geom_point(data = geochron_neo_quants, aes(x = depth, y = ymid), colour='#FF8C00', alpha=0.8) +
    # geom_linerange(data = geochron_neo_quants, aes(x = depth, ymin = ylo, ymax = yhi), colour='#FF8C00', alpha=0.8, lwd=1) +
    # # geom_point(data = controls, aes(x = depth, y = age)) +
    labs(title = paste0(dsid, '; ',  wang_fc$handle[idx_dsid]), x = "Depths (cm)", y = "Age", color = "Legend") +
    scale_color_manual(values = colors)
  
  print(p)
  

  
  ggsave(paste0('figures/age_depth_compare_', dsid, '.png'))
  ggsave(paste0('figures/age_depth_compare_', dsid, '.pdf'))
  
  dsid = datasetids[i]
  print(i)
  
  neo_plot = ggplot() +
    geom_ribbon(data = neo_quants, aes(x = depth, ymin = ylo, ymax = yhi), fill = "#FFA500AA") +
    geom_line(data = neo_quants, aes(x = depth, y = ymid))+
    # geom_point(data = neo_quants, aes(x = depth, y = ymid), colour='#FF8C00', alpha=0.8) +
    geom_point(data = geochron_neo_quants, aes(x = depth, y = ymid), colour='#FF8C00', alpha=0.8) +
    geom_linerange(data = geochron_neo_quants, aes(x = depth, ymin = ylo, ymax = yhi), colour='#FF8C00', alpha=0.8, lwd=1)+
    # geom_linerange(data = neo_quants, aes(x = depth, ymin = ylo, ymax = yhi), colour='#FF8C00', alpha=0.8, lwd=1)+
    labs(title = "Neotoma Ages", x = "Depth (cm)", y = "Age Mean")+
    labs(title = paste0(dsid, '; ', wang_fc$handle[idx_dsid]))
  
  
  bchron_plot = ggplot() +
    geom_ribbon(data = bchron_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#0000FF33") +
    geom_line(data = bchron_quants, aes(x = depths, y = ymid))+
    geom_point(data = geochron_quants, aes(x = depth-1, y = ymid), colour='#1F77B4', alpha=0.8) +
    geom_linerange(data = geochron_quants, aes(x = depth-1, ymin = ylo, ymax = yhi), colour='#1F77B4', alpha=0.8, lwd=1)+
    labs(title = "Bchron Ages Calibrated" , x = "Depth (cm)", y = "Age Mean")
  
  
  bacon_plot = ggplot() +
    geom_ribbon(data = wang_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#FF000033") +
    geom_line(data = wang_quants, aes(x = depths, y = ymid))+
    geom_point(data = geochron_bacon_quants, aes(x = depth+1, y = ymid), colour='#D62728', alpha=0.8) +
    geom_linerange(data = geochron_bacon_quants, aes(x = depth+1, ymin = ylo, ymax = yhi), colour='#D62728', alpha=0.8, lwd=1)+
    labs(title = "Bacon Ages Calibrated", x = "Depth (cm)", y = "Age Mean")
  
  tri_plot = plot_grid(neo_plot, bchron_plot, bacon_plot, olap_dots) 
  
  print(tri_plot)
  
  
  ggsave(paste0('figures/three_plot_', dsid, '.png'))
  ggsave(paste0('figures/three_plot_', dsid, '.pdf'))
  
  
  ###
  #new figure
  #ggsave()
  
  wang_mean = data.frame(depths=wang_posts[,'depths'], 
                         age_mean_w = rowMeans(wang_posts[,2:ncol(wang_posts)]),
                         age_sd_w = apply(wang_posts[,2:ncol(wang_posts)], 1, sd, na.rm=TRUE))
  
  bchron_mean = data.frame(depths=bchron_posts[,'depths'], 
                           age_mean_b = rowMeans(bchron_posts[,2:ncol(bchron_posts)]),
                           age_sd_b = apply(bchron_posts[,2:ncol(bchron_posts)], 1, sd, na.rm=TRUE))
  
  age_means = merge(bchron_mean, wang_mean, by = 'depths')
  age_means = merge(age_means, neo_mean, by = 'depths', all = TRUE)
  
  diffs = rbind(diffs,
                data.frame(dsid = rep(dsid),
                           age_means))
  
}


fnames = list.files('figures', 'age_depth_compare_.*.pdf', recursive=TRUE)

fname_str = sapply(fnames, function(x) paste0('figures/', x))
fname_str = paste(fname_str, collapse = ' ')

sys_str = paste0("gs -sDEVICE=pdfwrite -o age_depth_compare_v", vers, ".pdf ", fname_str)

system(sys_str)
