#  Note, put this code in a directory with a Bacon.R file and a Cores directory.
require(raster)
require(fields)
require(sp)
require(Bchron)
library(ggplot2)
library(parallel)
library(mgcv)
# 
# source('R/config.r')
# source('R/utils/helpers.r')

chron.control.types <- read.csv("chroncontrol_types-edited.csv")


radio = read.csv('data/radiocarbon-dates-errors.csv')
mod_radio <- gam(error ~ s(age, k=15), data=radio, method='REML', family=Gamma(link="identity"))

lead = read.csv('data/lead-dates-errors.csv')
mod_lead <- gam(error ~ s(age, k=4), data=lead, method='REML', family=Gamma(link="log"))

extrap = 1000

options(show.error.messages = TRUE)
#options(show.error.messages = FALSE)

version='200 - 250'

get_sitename <- function(core.id) {
  geochron <- try(read.table(paste0('Cores/', core.id, '/', core.id, '.csv'), sep=',', header=TRUE))
  if(is(geochron, "try-error")) {
    return(NA)
  }
  geochron$sitename[1]
}


fix_geochron_errors <- function(geochron){
  
  # fix some errors in data
  if (any(geochron$chroncontrolid== 105848)){
    geochron$age[which(geochron$chroncontrolid== 105848)] = -10
    geochron$limityounger[which(geochron$chroncontrolid== 105848)] = 1950 - geochron[which(geochron$chroncontrolid == 105848),'limityounger']
  }
  if (any(geochron$chroncontrolid== 105849)){
    geochron$age[which(geochron$chroncontrolid== 105849)] = -12
    geochron$limityounger[which(geochron$chroncontrolid== 105849)] = 1950 - geochron[which(geochron$chroncontrolid == 105849),'limityounger']
  }
  if (any(geochron$chroncontrolid== 104870)){
    geochron$type[which(geochron$chroncontrolid== 104870)] = 'Tephra'
  }
  
  # if (any(geochron$limityounger < (-70))){
  #   geochron$limityounger[which(geochron$limityounger < (-70))] = -70
  # }
  
 if (any((!is.na(geochron$limityounger))&(geochron$limityounger < (-70)))){
       geochron$limityounger[which(geochron$limityounger < (-70))] = -70
     }
  
  if (any((geochron$type == 'Lead-210') & (geochron$limitolder > 200))){
    geochron$limitolder[which((geochron$type == 'Lead-210') & (geochron$limitolder > 200))] = 200
  }
  if (any(is.na(geochron$age))){
    geochron$age[which(is.na(geochron$age))] = (geochron$limityounger[which(is.na(geochron$age))] + geochron$limitolder[which(is.na(geochron$age))])/2
  }
  
  return(geochron)
}

set_bio_flag <- function(chron.control.meta, keep, idx){
  
  bio.flag = c()
  if (any(keep != 1)) {
    for (i in 1:length(idx)){
      if (chron.control.meta$chron.control.type[idx[i]] == "Biostratigraphic, pollen") {
        keep[i] = 1 
        bio.flag = c(bio.flag, "Biostratigraphic, pollen")
      } else if (chron.control.meta$chron.control.type[idx[i]] == "Ambrosia rise") {
        keep[i] = 1 
        bio.flag = c(bio.flag, "Ambrosia rise")
      } else if (chron.control.meta$chron.control.type[idx[i]] == "European settlement horizon") {
        keep[i] = 1 
        bio.flag = c(bio.flag, "European settlement horizon, Pre-EuroAmerican settlement horizon")
      } else if (chron.control.meta$chron.control.type[idx[i]] == "Tsuga decline") {
        keep[i] = 1 
        bio.flag = c(bio.flag, "Tsuga decline")
      }
    }
  }
  
  return(bio.flag)
}


set_missing_control_errors <- function(geochron, mod_radio, mod_lead){
  
  if (any(is.na(geochron$error)|(geochron$error == 0))){
    error.na = which(is.na(geochron$error)|(geochron$error == 0))
    
    for (i in 1:length(error.na)){
      if (geochron$type[error.na[i]] %in% c('Radiocarbon', 'Radiocarbon, average of two or more dates', 'Radiocarbon, reservoir correction')) {
        geochron$error[error.na[i]] = round(mean(predict.gam(mod_radio, new_age=data.frame(age=geochron$age[error.na[i]]), type='response', se.fit=TRUE)$fit))
      } else if (geochron$type[error.na[i]] %in% c('Core top')) {
        geochron$error[error.na[i]] = 32#.27
        #} else if (geochron$type[error.na[i]] %in% c('Radiocarbon, average of two or more dates', 'Radiocarbon, reservoir correction')){
        #  geochron$error[error.na[i]] = 50
        # >>>>>>> 904b764fa902585d4b72f3c53a83ce22e4271314
      } else if (geochron$type[error.na[i]] %in% c('Annual laminations (varves)')) {
        geochron$error[error.na[i]] = (geochron$age[error.na[i]]+70)*0.05
      } else if (geochron$type[error.na[i]] %in% c('Biostratigraphic, pollen', 'Tsuga decline')) {
        geochron$error[error.na[i]] = 250
      } else if (geochron$type[error.na[i]] %in% c('Ambrosia rise', 'European settlement horizon')) {
        geochron$error[error.na[i]] = 50
      } else if (geochron$type[error.na[i]] %in% c('Tephra')) {
        geochron$error[error.na[i]] = 334
      } else if (geochron$type[error.na[i]] %in% c('Lead-210')) {
        geochron$error[error.na[i]] = round(mean(predict.gam(mod_lead, new_age=data.frame(age=geochron$radio[error.na[i]]), type='response', se.fit=TRUE)$fit))
      } else if (geochron$type[error.na[i]] %in% c('Other dating methods')) {
        geochron$error[error.na[i]] = 350
      } else {
        stop(paste0('an na chron control error could not be assigned for ', geochron$type[error.na[i]]))
      }
    }
  }
  
  return(geochron)
  
}


do_core_bchron <- function(core.id, chron.control.meta, mod_radio, mod_lead, extrap) {

  # read in the geochron table for core.id
  geochron = read.table(paste0('Cores/', core.id, '/', core.id, '.csv'), sep=',', header=TRUE)
  
  # remove duplicated rows (rare)
  geochron <- unique(geochron)
  
  # fix some specific errors
  geochron = fix_geochron_errors(geochron)
  
  # STOP if only one control
  if (nrow(geochron) == 1) {
    stop("only 1 control before filtering")
  }  
  
  # STOP if no controls
  if (nrow(geochron) < 1) {
    stop("no controls before filtering")
  }
  
  # match core.id control types to master list
  # STOP if any types not in master list
  idx <- match(geochron$type, chron.control.meta$chron.control.type, nomatch=NA)
  if (any(is.na(idx))) {
    stop(paste0(nrow(geochron), " controls; types not in master"))
  }
  
  # define an error column
  # half the distance between the limitolder and limityounger
  geochron$error = abs(geochron$limitolder-geochron$limityounger) / 2
  
  
  # idx should index all rows
  # above statement: STOP if any don't match master control type list
  # add the calibrated column to geochron file
  # make a keep vector; we will modify this according to our rules for inclusion
  geochron$cc = chron.control.meta$cc[idx]
  keep        = chron.control.meta$keep[idx]
  
  # if (is.na(keep)) {
  #   stop(paste0(nrow(geochron), " controls; nbut"))
  # }
  
  bio.flag = set_bio_flag(chron.control.meta, keep, idx)

  age.flag = NA
  if (any(is.na(geochron$age))) {
    age.flag = "had one or more NA ages, but ran anyway"
  }
  
  # remove rows from geochron that have NA ages
  # STOP if fewer than 2 ages left
  geochron = geochron[!is.na(geochron$age),]
  if (nrow(geochron) < 2) {
    stop("less than 2 non-NA ages")
  }
  
  # discard unreliable controls
  # only keep rows with keep flag of 1
  geochron = geochron[which(keep==1),]
  

  # only use site if 3 or more non-biostrat dates
  # or if biostrat date is ambrosia rise, use
  types = chron.control.meta$chron.control.type[match(geochron$type, chron.control.meta$chron.control.type, nomatch=NA)]
  if (any(types %in% c('Biostratigraphic, pollen', 'Tsuga decline'))){
    # if (sum(!(types %in% c('Biostratigraphic, pollen', 'Tsuga decline')))>= 3){
    #   stop("Biostrat control with fewer than 3 controls of other types")
    # }
    if (sum(!(types %in% c('Biostratigraphic, pollen', 'Tsuga decline')))>= 2){
      stop("Biostrat control with fewer than 2 controls of other types")
    }
  }
  
  # for any NA or 0 control errors, set them to informed value
  # positive error required for Bchron
  geochron = set_missing_control_errors(geochron)
  
  if (nrow(geochron) <= 1) {
    stop("one or no controls after filtering")
  }
  
  if (any(is.na(geochron$cc))) {
    stop('one or more cc values are NA')
  }
  
  if (any(geochron$error == 0)){
    stop('geochron error still set to zero')
  }

  #if (is.na(error)) next
  
  depths = scan(paste0('Cores/', core.id, '/', core.id, '_depths.txt'))
  if (diff(range(depths)) == 0) {
    stop('diff(range(depths)) is zero')
  }

  calCurves = rep(NA, nrow(geochron))
  calCurves[which(geochron$cc == 0)] = "normal"
  calCurves[which(geochron$cc == 1)] = "intcal20"

  if (any(is.na(calCurves))) {
    stop('one or more calcurve values are NA')
  }
  
  #geochron$error[which(geochron$error == 0)] = 1

  out = Bchronology(ages    = geochron$age,
                     ageSds = geochron$error, 
                     calCurves = calCurves,
                     positions = geochron$depth, 
                     positionThicknesses = rep(4,nrow(geochron)),
                     ids = geochron$chroncontrolid, 
                     predictPositions = depths)

  summary(out, type='convergence', na.rm=TRUE)
  
  fname = paste0('Cores/', core.id, '/', core.id, '_bchron.pdf')
  p <-  plot(out, ageScale = c('bp')) +
    labs(title =  paste0(geochron$sitename[1], '; Dataset ID ', core.id),
         x = 'Age (cal years BP)',
         y = 'Depth (cm)')
  ggsave(fname, plot=p)
    

  
  predict_geo = predict(out, 
                        newPositions = geochron$depth)
  
  post_geo = data.frame(labid=geochron$chroncontrolid, depths=geochron$depth, t(predict_geo))
  write.table(post_geo, paste0('.', '/Cores/', core.id, '/', 
                           core.id, '_bchron_geo_samples.csv'), sep=',', col.names = TRUE, row.names = FALSE)
 
  post_geo_means = rowMeans(t(predict_geo))
  post_geo_old = max(post_geo_means) + extrap
  post_geo_young = ifelse((min(post_geo_means) - extrap)<(-70), -70, min(post_geo_means) - extrap)
  
  post_sample_means = rowMeans(t(out$thetaPredict))
  idx_reliable = which((post_sample_means < post_geo_old) & (post_sample_means > post_geo_young))
  
  print(paste0('Core ' , core.id, ';', length(idx_reliable)/length(post_sample_means)))
  
  post_sample = data.frame(depths=depths[idx_reliable], matrix(t(out$thetaPredict)[idx_reliable,], nrow=length(idx_reliable)))
  
  write.table(post_sample, paste0('.', '/Cores/', core.id, '/', 
                           core.id, '_bchron_samples.csv'), sep=',', col.names = TRUE, row.names = FALSE)
   
  if (length(bio.flag) >= 1){
    stop(paste0('unreliable controls used: ', paste(bio.flag, collapse='; ' )))
  }

  if (!is.na(age.flag)) {
    stop(age.flag)
  }

  TRUE
}

# debugging
#do_core_bchron(1000, chron.control.types, mod)
#stop("debug stop")

mc.cores = 4 # processors

core.ids = list.files('Cores')
ncores = length(core.ids)

# core.ids = core.ids[1:4]
# ncores = length(core.ids)


# bchron.reports = mclapply(core.ids, function(core.id) {
#   message(core.id)
#   print(core.id)
#   bchron.report = data.frame(datasetid = core.id, sitename=get_sitename(core.id), success = 1, reason=NA)
#   success = try(do_core_bchron(core.id, chron.control.types, mod_radio, mod_lead, extrap))
#   if(is(success, "try-error")) {
#     bchron.report$success = 0
#     bchron.report$reason = geterrmessage()
#   }
#   bchron.report
# }, mc.cores=mc.cores)

bchron.report = data.frame(datasetid = numeric(0), sitename=character(0), success = numeric(0), reason=character(0))

# do in chunks cause busted
for (i in 201:250) {

  if (i==1429){next}
  
  core.id = core.ids[i]
  
  print(i)
  print(core.id)
  
  bchron.report.site = data.frame(datasetid = core.id, sitename=get_sitename(core.id), success = 1, reason=NA)
  success = try(do_core_bchron(core.id, chron.control.types, mod_radio, mod_lead, extrap))
  if(is(success, "try-error")) {
    bchron.report.site$success = 0
    bchron.report.site$reason = geterrmessage()
  }
  bchron.report = rbind(bchron.report, bchron.report.site)
  write.csv(bchron.report, paste0('bchron_report_v', version, '.csv'), row.names=FALSE)
}



 # bchron.report = do.call(rbind, bchron.reports)
write.csv(bchron.report, paste0('bchron_report_v', version, '.csv'), row.names=FALSE)

fnames = list.files('Cores', '*_bchron.pdf', recursive=TRUE)

fname_str = sapply(fnames, function(x) paste0('Cores/', x))
fname_str = paste(fname_str, collapse = ' ')

sys_str = paste0("gs -sDEVICE=pdfwrite -o bchron_plots_v", version, ".pdf ", fname_str)
system(sys_str)

