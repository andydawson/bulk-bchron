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
                   name = character(0),
                   depths = numeric(0),
                   age_mean_b  = numeric(0),
                   age_sd_b  = numeric(0),
                   age_mean_w  = numeric(0),
                   age_sd_w  = numeric(0),
                   age_mean_n  = numeric(0),
                   age_sd_n = numeric(0))

geo_diffs = data.frame(dsid = numeric(0),
                   depths = numeric(0),
                   geo_age_mean_b  = numeric(0),
                   geo_age_sd_b  = numeric(0),
                   geo_age_mean_w  = numeric(0),
                   geo_age_sd_w = numeric(0),
                   geo_age_mean_n  = numeric(0),
                   geo_age_sd_n  = numeric(0))

# pdf('figures/age_depth_compare.pdf', width=10, height=6)



for (i  in 575:N_datasetids){#N_datasetids){#N_datasetids){
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
  
  files = list.files(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid]))                
  idx_file = which(str_sub(files,-8,-1) == 'rout.csv')
  wang_posts = read.csv(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid], '/', files[idx_file]))
  
  wang_depths =  scan(paste0('Cores_bacon/Cores_full/', wang_fc$handle[idx_dsid], '/', wang_fc$handle[idx_dsid], '_depths.txt'))
  
  # wang_posts = data.frame(depths = wang_depths, wang_posts)
  wang_posts = data.frame(depths = wang_depths, wang_posts)#[,1:100])
  
  wang_posts_long = melt(wang_posts, id.vars = "depths")
  colnames(wang_posts_long) = c('depths', 'iter', 'age')

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
  bchron_geo_posts = read.csv(paste0('Cores/', dsid, '/', dsid, '_bchron_geo_samples.csv'))
  bchron_geo_posts = bchron_geo_posts[,-1]
  
  if (nrow(bchron_posts)==0){
    print(paste0('No posterior samples for dataset id: ', dsid))
    next
  }
  
  bchron_posts_long = melt(bchron_posts, id.vars = "depths")
  colnames(bchron_posts_long) = c('depths', 'iter', 'age')
  
  bchron_geo_posts_long = melt(bchron_geo_posts, id.vars = c("depths"))
  colnames(bchron_geo_posts_long) = c('depths', 'iter', 'age')
  
  # bchron
  
  # quantile(bchron_posts$V1, c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  bchron_quants_row = apply(bchron_posts[,2:ncol(bchron_posts)], 1, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  bchron_quants = data.frame(depths = bchron_posts[,1], t(bchron_quants_row))
  colnames(bchron_quants) = c('depths', 'ylo', 'ymid', 'yhi')
  
  bchron_geo_quants_row = apply(bchron_geo_posts[,2:ncol(bchron_geo_posts)], 1, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  bchron_geo_quants = data.frame(depths = bchron_geo_posts[,1], t(bchron_geo_quants_row))
  colnames(bchron_geo_quants) = c('depths', 'ylo', 'ymid', 'yhi')
  
  adjustcolor( "blue", alpha.f = 0.2)
  
  bchron_all_quants = rbind(bchron_quants, bchron_geo_quants)
  
  ggplot() +
    geom_ribbon(data = bchron_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#0000FF33") +
    geom_line(data = bchron_quants, aes(x = depths, y = ymid))
  
  ggplot() +
    geom_ribbon(data = bchron_geo_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#0000FF33") +
    geom_line(data = bchron_geo_quants, aes(x = depths, y = ymid))
  
  # I'm not sure we want to bind these together
  # need to predict ages for all depths that cover samples and controls
  ggplot() +
    geom_ribbon(data = bchron_all_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#0000FF33") +
    geom_line(data = bchron_all_quants, aes(x = depths, y = ymid))
  
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
    
  } else {
    neo_quants = data.frame(depth = neo_site$depth,
                            ylo = neo_site$age,
                            ymid = neo_site$age,
                            yhi = neo_site$age)
    
    neo_mean = data.frame(depths = neo_site$depth,
                          age_mean_n = neo_site$age,
                          age_sd_n = 0)
    
    
  }
  
  # Now bchron, Bacon, and Neotoma? on one plot #
  
  colors = c("Bchron" = "blue", "Bacon" = "red", "Neotoma" = "orange")
  
  p <- ggplot() +
    geom_ribbon(data = bchron_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#0000FF33") +
    geom_line(data = bchron_quants, aes(x = depths, y = ymid, color = "Bchron"), size = 1.5) +
    geom_ribbon(data = wang_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "#FF000033") +
    geom_line(data = wang_quants, aes(x = depths, y = ymid, color = "Bacon"), size = 1.5) + 
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
    geom_point(data = neo_quants, aes(x = depth, y = ymid), colour='#FF8C00', alpha=0.8) +
    geom_linerange(data = neo_quants, aes(x = depth, ymin = ylo, ymax = yhi), colour='#FF8C00', alpha=0.8, lwd=1)+
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
  
  tri_plot = plot_grid(neo_plot, bchron_plot, bacon_plot) 
  
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
               data.frame(dsid = rep(dsid, nrow(age_means)),
                          name = rep(wang_fc$handle[idx_dsid], nrow(age_means)),
                           age_means))
  
 }
# dev.off()

diffs = diffs[order(diffs$dsid),]


##For creating subset of usable data. Don't overwrite data!!##
#compare_sites_ids = diffs[!duplicated(diffs$dsid),]
#compare_sites_ids$visual_check = 1


#write.csv(compare_sites_ids[, c('dsid', 'name','visual_check')], 'compare_sites_ids.csv', row.names = FALSE)


fnames = list.files('figures', 'age_depth_compare_.*.pdf', recursive=TRUE)

fname_str = sapply(fnames, function(x) paste0('figures/', x))
fname_str = paste(fname_str, collapse = ' ')

sys_str = paste0("gs -sDEVICE=pdfwrite -o age_depth_compare_v", vers, ".pdf ", fname_str)

system(sys_str)

# 

model_check = read.csv('compare_sites_ids.csv', stringsAsFactors = FALSE)
dsid_pass = model_check[which(model_check$pass ==1), 'dsid']

diffs = diffs[which(diffs$dsid %in% dsid),]



diffs$diff_bb = abs(diffs$age_b - diffs$age_w)
diffs$diff_bn = abs(diffs$age_b - diffs$age_n)

med_diff = diffs[which(diffs$diff_bb >= 7000 & diffs$diff_bb <= 25000),]

which(diffs$diff_bb >= 25000)
lrg_diff = diffs[which(diffs$diff_bb >= 25000),]

sml_diff = diffs[which(diffs$diff_bb <=7000),]

low_diff = diffs[which(diffs$diff_bb >= 0 & diffs$diff_bb <=1),]


#use these sites
#use_diff = diffs[which(diffs$diff_bb <= 10000 & diffs$diff_bb >= 0),]

ggplot(data = use_diff) +
  geom_histogram(aes(x=diff_bb, y=..density..), bins=100)

ggplot(data = sml_diff) +
  geom_histogram(aes(x=diff_bb, y=..density..), bins=100)

site_diffs = diffs %>% group_by(dsid) %>% summarize(site_diff_bb = mean(diff_bb),  site_diff_bn = mean(diff_bn))



summary(diffs$diff_bb)

  

plot_idx = data.frame(radio$datasetid)


# ggplot(data = diffs) +
#   geom_point(aes(x = depths, y = age_b), color = 'blue', alpha = .2) +
#   geom_point(aes(x = depths, y = age_w), color = 'red', alpha = .2)






#Compare age means and SDs

bacon_comp = plot (diffs$age_mean_w, diffs$age_sd_w)

bchron_comp = plot (diffs$age_mean_b, diffs$age_sd_b)

ggplot(data=diffs) + geom_point(aes(x=age_mean_b, y=age_sd_b, colour= factor(dsid)))

ggplot(data=diffs) + geom_point(aes(x=age_mean_w, y=age_sd_w, colour= factor(dsid)))



ggplot(data=subset(diffs, dsid==15356)) + geom_point(aes(x=age_mean_b, y=age_sd_b), colour='blue') + geom_point(aes(x=age_mean_w, y=age_sd_w), colour='black')+
  geom_line(aes(x=age_mean_b, y=age_sd_b), colour='blue') + geom_line(aes(x=age_mean_w, y=age_sd_w), colour='black')


legend_1 = c("Bchron" = "blue", "Bacon" = "black")

ggplot(data=subset(diffs, dsid==1136)) + geom_point(aes(x=age_mean_b, y=age_sd_b, colour="Bchron")) + geom_point(aes(x=age_mean_w, y=age_sd_w, colour="Bacon")) +
  geom_line(aes(x=age_mean_b, y=age_sd_b, colour="Bchron")) + geom_line(aes(x=age_mean_w, y=age_sd_w, colour="Bacon")) +
  labs(title = "Age Mean vs. SD: 1136", x = "age mean", y = "age sd", color = "Legend") +
  scale_colour_manual(values = legend_1)


ggplot(data=diffs) + geom_point(aes(x=age_sd_w, y=age_sd_b)) +
  geom_abline(slope = 1, intercept = 0) +
  coord_equal() +
  xlim(c(0,800)) + ylim(c(0,800)) 


ggplot(data=subset(diffs, dsid ==14680)) + geom_point(aes(x=age_sd_w, y=age_sd_b)) +
  geom_abline(slope = 1, intercept = 0) +
  coord_equal() +
  xlim(c(0,800)) + ylim(c(0,800)) 


goo = diffs$age_sd_b - diffs$age_sd_w
bacon_bigger = which(goo<0)
length(bacon_bigger)
#I suspect this number is so large because of the errors with the memory strength 

bchron_bigger = which (goo>0)
length(bchron_bigger)

same = which(goo == 0)


dsid = datasetids[i]
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



pdf('figures/olap_figs.pdf', width=10, height=6)
for (i in 1:n_depths){
  
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
dev.off()








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
