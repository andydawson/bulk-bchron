library(stringr)
library(ggplot2)
library(overlapping)
library(reshape2)

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
  
  # wang_posts = data.frame(depths = wang_depths, wang_posts)
  wang_posts = data.frame(depths = wang_depths, wang_posts[,2:101])
  
  wang_posts_long = melt(wang_posts, id.vars = "depths")
  colnames(wang_posts_long) = c('depths', 'iter', 'age')
  
  # wang_posts_long$iter 
  
  ggplot() +
    geom_point(data = wang_posts, aes(x = depths, y = V1)) +
    geom_point(data = wang_posts, aes(x = depths, y = V2))
  
  ggplot() +
    geom_point(data = wang_posts_long, aes(x = depths, y = age))
  
  ggplot() +
    geom_line(data = wang_posts_long, aes(x = depths, y = age))
  
  ggplot() +
    geom_line(data = wang_posts_long, aes(x = depths, y = age, group = iter))
  
  
  quantile(wang_posts$V1, c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  wang_quants_row = apply(wang_posts[,2:ncol(wang_posts)], 1, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  wang_quants = data.frame(depths = wang_posts[,1], t(wang_quants_row))
  colnames(wang_quants) = c('depths', 'ylo', 'ymid', 'yhi')
  
  
  ggplot() +
    geom_ribbon(data = wang_quants, aes(x = depths, ymin = ylo, ymax = yhi), fill = "pink") +
    geom_line(data = wang_quants, aes(x = depths, y = ymid))
    
  
  
  ###
  
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
