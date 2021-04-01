load('data/pollen_north_america_p25.rdata')


# pollen_north_america


core.ids = list.files('Cores')
ncores = length(core.ids)

bchron_report = read.csv(paste0('bchron_report_v3.0.csv'))

bchron_report$success[which(grepl("unreliable controls", bchron_report$reason, fixed=TRUE))] = 1


for (i in 1:ncores) {
  print(paste0("Core ", i, " of ", ncores, " cores"))
  idx_report = which(bchron_report$datasetid == core.ids[i])
  
  idx_pollen = which(pollen_north_america$dataset_id == core.ids[i])
  
  if (bchron_report$success[idx_report]==1){
    
    samps = read.csv(paste0('Cores/', core.ids[i], '/', core.ids[i], '_bchron_samples.csv'))
    age_mean = apply(samps[, 2:ncol(samps)], 1, median)
    
    # idx_in = which(!is.na(match(pollen_north_america$depth[idx_pollen], samps$depths)))
    # idx_out = 
    
    pollen_north_america$age[idx_pollen] = age_mean[match(pollen_north_america$depth[idx_pollen], samps$depths)]
    
    if (!all(pollen_north_america$depth[idx_pollen]  %in% samps$depths)){
      print(paste0("Problem with dataset ", pollen_north_america$sitename[idx_pollen][1]))
    } #else { 
      #pollen_north_america$age[idx_pollen] = age_mean 
    #}
    
  } else {
    
    next
  }
}

write.csv(pollen_north_america, 'data/pollen_north_america_v2.0.csv', row.names = FALSE)
