library(plyr)
library(ggplot2)
library(Bchron)

version = '7.0'

# type can be LC6K or ABI
type = "LC6K"

if (type == "LC6K"){
  #load('data/pollen_north_america_lc6k.rdata')
  pollen_north_america <- read.csv("data/pollen_north_america_lc6k.csv")
} else if (type == "ABI"){
  #load('data/pollen_north_america_p25.rdata')
  pollen_north_america <- read.csv("data/pollen_north_america_abi.csv")
}

pollen_north_america = pollen_north_america[which(!is.na(pollen_north_america$agetype)),]

site_meta = read.csv('data/sites_north_america.csv')

extrap = 1000

cc = read.csv('data/chroncontrol_summary_pollen_full.csv')

ncontrols = ddply(cc, .(datasetid), nrow)

id_keep = ncontrols$datasetid[which(ncontrols$V1>1)]

cc = cc[which(cc$datasetid %in% id_keep),]

cc[which(cc$chroncontrolid == 105848),'age'] = -10
cc[which(cc$chroncontrolid == 105849),'age'] = 12
cc[which(cc$chroncontrolid == 105848),'limityounger'] = 1950 - cc[which(cc$chroncontrolid == 105848),'limityounger']
cc[which(cc$chroncontrolid == 105849),'limityounger'] = 1950 - cc[which(cc$chroncontrolid == 105849),'limityounger']

cc[which(cc$chroncontrolid == 104870),'type'] = 'Tephra'

cc[which(cc$limityounger < (-70)), 'limityounger'] = -70
cc[which((cc$type=="Lead-210")&(cc$limitolder>200)), 'limitolder'] = 200

cc[which(is.na(cc$age)),'age'] = (cc[which(is.na(cc$age)),'limityounger'] + cc[which(is.na(cc$age)),'limitolder']) / 2

cc = cc[which(!is.na(cc$age)),]

rc_types = cc$type %in% c('Radiocarbon', 
                          'Radiocarbon, reservoir correction', 
                          'Radiocarbon, average of two or more dates', 
                          'Radiocarbon, infinite',
                          'Radiocarbon, calibrated from calendar years',
                          'Radiocarbon, calibrated')
cc$geoagetype[which(rc_types)] = 'Radiocarbon years BP'

cc[which(cc$chroncontrolid == 37261),'geoagetype'] = 'Years BP'
cc[which(cc$chroncontrolid == 50647),'geoagetype'] = 'Years BP'

ar_types = cc$type %in% c('Ambrosia rise')
cc$geoagetype[which(ar_types)] = NA

cc$geoagetype[which(is.na(cc$geoagetype))] = "Years BP"


# translate any dates in radiocarbon years to calendar years
# radio.years <- ((cc$type %in% c("Radiocarbon years BP", "Calibrated radiocarbon years BP"))&
#   (cc$age > 95 ) &
#   (cc$age < 50193))
# sryears <- sum(radio.years, na.rm = TRUE)
# 
# # BChronCalibrate is in the BChron package:
# calibrated <- BchronCalibrate(compiled.cores$age[radio.years],
#                               ageSds = rep(100, sryears),
#                               calCurves = rep("intcal20", sryears))
# #  we want the weighted means from "calibrated"
# wmean.date <- function(x) sum(x$ageGrid*x$densities / sum(x$densities))
# compiled.cores$age[radio.years] <- sapply(calibrated, wmean.date)
# # saveRDS(compiled.cores, 'data/compiled_cores_P25_all_times.RDS')


# pollen_north_america
core.ids = sort(unique(pollen_north_america$dataset_id))
ncores = length(core.ids)

# core.ids = list.files('Cores')
# ncores = length(core.ids)

bchron_report = read.csv(paste0('bchron_report_v', version, '.csv'))

bchron_report$success[which(grepl("unreliable controls", bchron_report$reason, fixed=TRUE))] = 1

pollen_north_america = data.frame(checked = rep(0, length=nrow(pollen_north_america)), 
                                  keep = rep(0, length=nrow(pollen_north_america)), 
                                  pollen_north_america)

count = 0 
for (i in 1:ncores) {
  #print(paste0("Core ", i, " of ", ncores, " cores"))
  
  
  idx_pollen = which(pollen_north_america$dataset_id == core.ids[i])
  
  if (core.ids[i] %in% bchron_report$datasetid){
    
    idx_report = which(bchron_report$datasetid == core.ids[i])
    
    if (bchron_report$success[idx_report]==1){
      
      samps = read.csv(paste0('Cores/', core.ids[i], '/', core.ids[i], '_bchron_samples.csv'))
      if (nrow(samps)==0){
        pollen_north_america$checked[idx_pollen] = 1
        next
      }
      
      age_mean = apply(samps[, 2:ncol(samps)], 1, median)
      
      # idx_in = which(!is.na(match(pollen_north_america$depth[idx_pollen], samps$depths)))
      # idx_out = 
      
      for (k in 1:length(idx_pollen)){
        if (length(which(samps$depth == pollen_north_america$depth[idx_pollen[k]]))!=0){
          
          if (length(age_mean[which(samps$depth == pollen_north_america$depth[idx_pollen[k]])]) != 1){
            print(paste0('Dataset ', i))
          }
          pollen_north_america$age[idx_pollen[k]] = age_mean[which(samps$depth == pollen_north_america$depth[idx_pollen[k]])][1]
          pollen_north_america$keep[idx_pollen[k]] = 1
          pollen_north_america$checked[idx_pollen[k]] = 1
        } else {
          pollen_north_america$checked[idx_pollen[k]] = 1
        }
      }
      # if (!all(pollen_north_america$depth[idx_pollen]  %in% samps$depths)){
      #   print(paste0("Problem with dataset ", pollen_north_america$sitename[idx_pollen][1]))
      # } #else { 
      # #pollen_north_america$age[idx_pollen] = age_mean 
      # #}
      
    } else if (bchron_report$success[idx_report]==0) {
      
      pollen_north_america$checked[idx_pollen] = rep(1, length(idx_pollen))
      
    } 
    
  } else if (pollen_north_america$agetype[idx_pollen][1] == 'Varve years BP'){
    
    pollen_north_america$keep[idx_pollen] = rep(1, length(idx_pollen))
    pollen_north_america$checked[idx_pollen] = rep(1, length(idx_pollen))
    
  } else if (pollen_north_america$dataset_id[idx_pollen][1] == '15320') { 
    
    pollen_north_america$checked[idx_pollen] = rep(1, length(idx_pollen))
    
    next 
    
  } else {
  
    control = cc[which(cc$datasetid == pollen_north_america$dataset_id[idx_pollen[1]]),]
    
    if ((nrow(control)==0)&(length(idx_pollen)==1)){
      pollen_north_america$keep[idx_pollen] = rep(1, length(idx_pollen))
      pollen_north_america$checked[idx_pollen] = rep(1, length(idx_pollen))
      next
    }
    if ((nrow(control) == 0)&(length(idx_pollen)>1)) {
      count = count + 1
      control_young = pollen_north_america$ageboundyounger[idx_pollen][1]
      control_old = pollen_north_america$ageboundolder[idx_pollen][1]
      
      idx_keep = which((pollen_north_america$age[idx_pollen]<control_old)&(pollen_north_america$age[idx_pollen]>control_young))
      pollen_north_america$keep[idx_pollen[idx_keep]] = rep(1, length(idx_keep))
      
      pollen_north_america$checked[idx_pollen] = rep(1, length(idx_pollen))
      
      next
      
    } else if (all(is.na(control$depth))) {
      control_min_depth = control[which.min(control$age),]
      control_max_depth = control[which.max(control$age),]
    } else {
      control_min_depth = control[which.min(control$depth),]
      control_max_depth = control[which.max(control$depth),]
    }
    
    
    if (control_min_depth$geoagetype == 'Radiocarbon years BP'){
      
      calibrated <- BchronCalibrate(control_min_depth$age,
                                    ageSds = 100,
                                    calCurves = "intcal20")
      # #  we want the weighted means from "calibrated"
      wmean.date <- function(x) sum(x$ageGrid*x$densities / sum(x$densities))
      control_young = wmean.date(calibrated)
      
    } else if (control_min_depth$geoagetype == 'Calendar years AD/BC'){
      
      control_young = 1950 - control_min_depth$age
    
    } else {
      
      control_young = control_min_depth$age
      
    }
    
    
    if (control_max_depth$geoagetype == 'Radiocarbon years BP'){
      
      calibrated <- BchronCalibrate(control_max_depth$age,
                                    ageSds = 100,
                                    calCurves = "intcal20")
      # #  we want the weighted means from "calibrated"
      wmean.date <- function(x) sum(x$ageGrid*x$densities / sum(x$densities))
      control_old = wmean.date(calibrated$Date1)
      
    } else if (control_max_depth$geoagetype == 'Calendar years AD/BC'){
      
      control_old = 1950 - control_max_depth$age
      
    } else {
      
      control_old = control_max_depth$age
      
    }
    
    
    control_old   = control_old + extrap
    control_young = ifelse((control_young - extrap) < -70, -70, control_young - extrap)

    idx_keep = which((pollen_north_america$age[idx_pollen]<control_old)&(pollen_north_america$age[idx_pollen]>control_young))

    pollen_north_america$keep[idx_pollen[idx_keep]] = rep(1, length(idx_keep))
    pollen_north_america$checked[idx_pollen] = rep(1, length(idx_pollen))
    
  }

}

sum(pollen_north_america$checked==1)/length(pollen_north_america$checked)
sum(pollen_north_america$keep==1)/length(pollen_north_america$keep)


if (type == "LC6K"){
  write.csv(pollen_north_america, paste0('data/pollen_north_america_LC6K_v', version, '.csv'), row.names = FALSE)
} else if (type == "ABI"){
  write.csv(pollen_north_america, paste0('data/pollen_north_america_ABI_v', version, '.csv'), row.names = FALSE)
}

