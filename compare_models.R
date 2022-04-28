library(stringr)

wang_fc = read.csv('wang/SiteInfo_fullcore.csv', stringsAsFactors = FALSE)

wang_cores  = list.files('wang/Cores_full')
wang_ncores = length(wang_cores)


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
