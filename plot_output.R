library(stringr)
library(ggplot2)
library(overlapping)
library(reshape2)
library(Bchron)
library(cowplot)
library(mgcv)
library(dplyr)


vers = 1.0

##########################################################################################################
## read and summarize biotic velocities
##########################################################################################################

# radio = read.csv('data/radiocarbon-dates-errors.csv')
# mod_radio <- gam(error ~ s(age, k=15), data=radio, method='REML', family=Gamma(link="identity"))
# 
# chron_control_types <- read.csv("data/chroncontrol_types-edited.csv")

wang_fc = read.csv('Cores_bacon/SiteInfo_fullcore.csv', stringsAsFactors = FALSE)
wang_cores  = list.files('Cores_bacon/Cores_full')
wang_ncores = length(wang_cores)

site_meta = read.csv('bchron_report_v9.0.csv', stringsAsFactors = FALSE)
neo_dat = read.csv("data/pollen_north_america_lc6k2_v8.0.csv")

site_meta = read.csv('bchron_report_v9.0.csv', stringsAsFactors = FALSE)
site_meta = site_meta[which(site_meta$success == 1),]

datasetids = unique(site_meta$datasetid)
N_datasetids = length(datasetids)

which(wang_fc$datasetid %in% site_meta$datasetid)

geochron_neo = read.csv('data/chroncontrol_summary_v2.csv')

n_neo = length(geochron_neo$depth)
geochron_neo$error = (geochron_neo$limitolder - geochron_neo$limityounger)

##########################################################################################################
## read in mean age differences and age SDs
##########################################################################################################

# age estimate differences for pollen sample depths
diffs = readRDS('diffs.RDS')

# age estimate differences for chron control depths
geo_diffs = readRDS('geo_diffs.RDS')

foo = diffs
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

model_check = read.csv('compare_sites_ids.csv', stringsAsFactors = FALSE, header = TRUE)
dsid_pass = model_check[which(model_check$visual_check ==1), 'dsid']

diffs = diffs[which(diffs$dsid %in% dsid_pass),]
geo_diffs = geo_diffs[which(geo_diffs$dsid %in% dsid_pass),]


diffs$diff_bb = (diffs$age_mean_b - diffs$age_mean_w)
diffs$diff_bn = (diffs$age_mean_w - diffs$age_mean_n)

geo_diffs$geo_diff_bb = (geo_diffs$geo_age_mean_b - geo_diffs$geo_age_mean_w)

#geo_diffs$geo_diff_bn = (geo_diffs$geo_age_mean_w - geo_diffs$geo_age_mean_n)


# med_diff = diffs[which(diffs$diff_bb >= 7000 & diffs$diff_bb <= 25000),]
# 
# which(diffs$diff_bb >= 5000)
# lrg_diff = diffs[which(diffs$diff_bb >= 25000),]
# 
# sml_diff = diffs[which(diffs$diff_bb <=7000),]
# 
# low_diff = diffs[which(diffs$diff_bb >= 0 & diffs$diff_bb <=1),]


#use these sites
#use_diff = diffs[which(diffs$diff_bb <= 10000 & diffs$diff_bb >= 0),]

use_diff = diffs

ggplot(data = use_diff) +
  geom_histogram(aes(x=diff_bb, y=after_stat(density)), bins=100)

#ggplot(data = sml_diff) +
#geom_histogram(aes(x=diff_bb, y=..density..), bins=100)

site_diffs = diffs %>% group_by(dsid) %>% summarize(site_diff_bb = mean(diff_bb),  site_diff_bn = mean(diff_bn))



summary(diffs$diff_bb)



plot_idx = data.frame(radio$datasetid)


# ggplot(data = diffs) +
#   geom_point(aes(x = depths, y = age_b), color = 'blue', alpha = .2) +
#   geom_point(aes(x = depths, y = age_w), color = 'red', alpha = .2)






#Compare age means and SDs
diffs = subset(diffs, dsid != 14104)

diffs[which((!is.na(diffs$age_sd_w))&(diffs$age_sd_w<30)),]


bacon_comp = plot(diffs$age_mean_w, diffs$age_sd_w)

bchron_comp = plot(diffs$age_mean_b, diffs$age_sd_b)

ggplot(data=diffs) + geom_point(aes(x=age_mean_b, y=age_sd_b, colour= factor(dsid)))

ggplot(data=diffs) + geom_point(aes(x=age_mean_w, y=age_sd_w, colour= factor(dsid)))

legend_1 = c("Bchron" = "blue", "Bacon" = "black")



# ggplot(data=subset(diffs, dsid==15356)) + geom_point(aes(x=age_mean_b, y=age_sd_b), colour='Bchron') + geom_point(aes(x=age_mean_w, y=age_sd_w), colour='black')+
#   geom_line(aes(x=age_mean_b, y=age_sd_b), colour='Bchron') + geom_line(aes(x=age_mean_w, y=age_sd_w), colour='Bacon')
# 


ggplot(data=subset(diffs, dsid==983)) + geom_point(aes(x=age_mean_b, y=age_sd_b, colour="Bchron")) + geom_point(aes(x=age_mean_w, y=age_sd_w, colour="Bacon")) +
  geom_line(aes(x=age_mean_b, y=age_sd_b, colour="Bchron")) + geom_line(aes(x=age_mean_w, y=age_sd_w, colour="Bacon")) +
  labs(title = "Age Mean vs. SD: 983", x = "age mean", y = "age sd", color = "Legend") +
  scale_colour_manual(values = legend_1)

fit_sd = lm(age_sd_b ~ age_sd_w, data=diffs)
summary(fit_sd)
fit_intercept = coef(fit_sd)[1]
fit_slope = coef(fit_sd)[2]

ggplot(data=diffs) + geom_point(aes(x=age_sd_w, y=age_sd_b)) +
  geom_abline(slope = 1, intercept = 0) +
  coord_equal() +
  #xlim(c(0,800)) + ylim(c(0,800)) +
  geom_abline(slope = fit_slope, intercept = fit_intercept, colour = 'red')

#diffs[diffs$age_sd_w>2000,'age_sd_w']

diffs[which((!is.na(diffs$age_sd_w))&(diffs$age_sd_w>2000)),]

fit = lm(age_sd_b ~ age_sd_w, data=diffs)

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

##The same plots but now with geo_diffs data
geo_diffs = subset(geo_diffs, dsid != 14104)

geo_bacon_comp = plot(geo_diffs$geo_age_mean_w, geo_diffs$geo_age_sd_w)
geo_bchron_comp = plot(geo_diffs$geo_age_mean_b, geo_diffs$geo_age_sd_b)

ggplot(data=subset(geo_diffs, dsid==1000)) + geom_point(aes(x=geo_age_mean_b, y=geo_age_sd_b, colour="Bchron")) + geom_point(aes(x=geo_age_mean_w, y=geo_age_sd_w, colour="Bacon")) +
  geom_line(aes(x=geo_age_mean_b, y=geo_age_sd_b, colour="Bchron")) + geom_line(aes(x=geo_age_mean_w, y=geo_age_sd_w, colour="Bacon")) +
  labs(title = "Geo Age Mean vs. SD: 1000", x = "age mean", y = "age sd", color = "Legend") +
  scale_colour_manual(values = legend_1)


ggplot(data = geo_diffs) +
  geom_histogram(aes(x=geo_diff_bb, y=after_stat(density)), bins=100)

ggplot(data = geo_diffs) +
  geom_point (aes(x=geo_age_mean_b, y=geo_diff_bb)) +
  labs(title = "Age Estimate Differences Between Models at Chronological Control", x = "Age Mean", y = "Difference")

ggplot(data = diffs) +
  geom_point (aes(x=age_mean_b, y=diff_bb)) +
  labs(title = "Age Estimate Differences Between Models", x = "Age Mean", y = "Difference")

geo_fit_sd = lm(geo_age_sd_b ~ geo_age_sd_w, data=geo_diffs)
summary(geo_fit_sd)
geo_fit_intercept = coef(geo_fit_sd)[1]
geo_fit_slope = coef(geo_fit_sd)[2]

ggplot(data=geo_diffs) + geom_point(aes(x=geo_age_sd_w, y=geo_age_sd_b)) +
  geom_abline(slope = 1, intercept = 0) +
  coord_equal() +
  #xlim(c(0,800)) + ylim(c(0,800)) +
  geom_abline(slope = geo_fit_slope, intercept = geo_fit_intercept, colour = 'red')

geo_bw_diff = geo_diffs$geo_age_sd_b - geo_diffs$geo_age_sd_w

geo_bchron_bigger = which (geo_bw_diff>0)
length(geo_bchron_bigger)

geo_bacon_bigger = which (geo_bw_diff<0)
length(geo_bacon_bigger)


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
