library(stringr)
library(ggplot2)



wang_fc = read.csv('wang/SiteInfo_fullcore.csv', stringsAsFactors = FALSE)

wang_cores  = list.files('wang/Cores_full')
wang_ncores = length(wang_cores)

site_meta = read.csv('bchron_report_v9.0.csv', stringsAsFactors = FALSE)


# there are 444 values for datasetid in wang_fc
which(wang_fc$datasetid %in% site_meta$datasetid)


site_meta = read.csv('bchron_report_v9.0.csv', stringsAsFactors = FALSE)
site_meta = site_meta[which(site_meta$success == 1),]

datasetids = unique(site_meta$datasetid)
N_datasetids = length(datasetids)

which(wang_fc$datasetid %in% site_meta$datasetid)



diffs = data.frame(dsid = numeric(0),
                   depths = numeric(0),
                   age_b  = numeric(0),
                   age_w  = numeric(0))

for (i in 1:N_datasetids){
  
  dsid = datasetids[i]
  
  idx_dsid = which(wang_fc$datasetid == dsid)
  
  if (length(idx_dsid) == 0){
    print(paste0('No Bacon age-depth model for dataset id: ', dsid))
    next
  }
  
  files = list.files(paste0('wang/Cores_full/', wang_fc$handle[idx_dsid]))                
  idx_file = which(str_sub(files,-8,-1) == 'rout.csv')
  wang_posts = read.csv(paste0('wang/Cores_full/', wang_fc$handle[idx_dsid], '/', files[idx_file]))
  
  wang_depths =  scan(paste0('wang/Cores_full/', wang_fc$handle[idx_dsid], '/', wang_fc$handle[idx_dsid], '_depths.txt'))
  
  wang_posts = data.frame(depths = wang_depths, wang_posts)
  wang_mean = data.frame(depths=wang_posts[,'depths'], age_w = rowMeans(wang_posts[,2:ncol(wang_posts)]))
  
  bchron_posts = read.csv(paste0('Cores/', dsid, '/', dsid, '_bchron_samples.csv'))
  
  if (nrow(bchron_posts)==0){
    print(paste0('No posterior samples for dataset id: ', dsid))
    next
  }
  
  bchron_mean = data.frame(depths=bchron_posts[,'depths'], age_b = rowMeans(bchron_posts[,2:ncol(bchron_posts)]))
  
  age_means = merge(bchron_mean, wang_mean)
  
  diffs = rbind(diffs,
                data.frame(dsid = rep(dsid),
                           age_means))

}


ggplot(data = diffs) +
  geom_point(aes(x = depths, y = age_b), color = 'blue', alpha = .2) +
  geom_point(aes(x = depths, y = age_w), color = 'red', alpha = .2)










# ###Try overlap function###
# library(reshape2)
# library(overlapping)

# olap = data.frame(datasetid = numeric(0),
#                   depths = numeric(0),
#                   overlap = numeric(0))
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
