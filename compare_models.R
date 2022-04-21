wang_fc = read.csv('wang/SiteInfo_fullcore.csv', stringsAsFactors = FALSE)

site_meta = read.csv('bchron_report_v8.0.csv', stringsAsFactors = FALSE)

which(wang_fc$datasetid %in% site_meta$datasetid)



wang_cores  = list.files('wang/Cores_full')
wang_ncores = length(wang_cores)



dsid = 48879

foo = read.csv(paste0('Cores/', dsid, '/', dsid, '_bchron_samples.csv'), stringsAsFactors = FALSE)

# for (i in 1:)
  
  which(wang_fc$handle == wang_cores[1])


# }