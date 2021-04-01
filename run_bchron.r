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
mod <- gam(error ~ s(age, k=15), data=radio, method='REML', family=Gamma(link="identity"))


options(show.error.messages = TRUE)
#options(show.error.messages = FALSE)

version='3.0'

get_sitename <- function(core.id) {
  geochron <- try(read.table(paste0('Cores/', core.id, '/', core.id, '.csv'), sep=',', header=TRUE))
  if(is(geochron, "try-error")) {
    return(NA)
  }
  geochron$sitename[1]
}

do_core_bchron <- function(core.id, chron.control.meta, mod) {

  geochron = read.table(paste0('Cores/', core.id, '/', core.id, '.csv'), sep=',', header=TRUE)
  geochron <- unique(geochron)
  
  geochron$error = abs(geochron$limitolder-geochron$limityounger) / 2
  
  if (nrow(geochron) == 1) {
    stop("only 1 control before filtering")
  }  
  
  if (nrow(geochron) < 1) {
    stop("no controls before filtering")
  }
  
  idx <- match(geochron$type, chron.control.meta$chron.control.type, nomatch=NA)
  if (any(is.na(idx))) {
    stop(paste0(nrow(geochron), " controls; types not in master"))
  }

  geochron$cc <- chron.control.meta$cc[idx]
  keep = chron.control.meta$keep[idx]
  
  # if (is.na(keep)) {
  #   stop(paste0(nrow(geochron), " controls; nbut"))
  # }
  
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
        bio.flag = c(bio.flag, "European settlement horizon")
      } else if (chron.control.meta$chron.control.type[idx[i]] == "Tsuga decline") {
        keep[i] = 1 
        bio.flag = c(bio.flag, "Tsuga decline")
    }
    }
  }

  age.flag = NA
  if (any(is.na(geochron$age))) {
    age.flag = "had one or more NA ages, but ran anyway"
  }
  geochron = geochron[!is.na(geochron$age),]
  if (nrow(geochron) < 2) {
    stop("less than 2 non-NA ages")
  }
  
  geochron = geochron[which(keep==1),]
  
  if (any(is.na(geochron$error)|(geochron$error == 0))){
    error.na = which(is.na(geochron$error)|(geochron$error == 0))
    
    for (i in 1:length(error.na)){
      if (geochron$type[error.na[i]] %in% c('Radiocarbon')) {
        geochron$error[error.na[i]] = round(mean(predict.gam(mod, new_age=data.frame(age=geochron$age[error.na[i]]), type='response', se.fit=TRUE)$fit))
      } else if (geochron$type[error.na[i]] %in% c('Core top')) {
        geochron$error[error.na[i]] = 33#.27
      } else if (geochron$type[error.na[i]] %in% c('Radiocarbon, average of two or more dates', 'Radiocarbon, reservoir correction')){
        geochron$error[error.na[i]] = 50
      } else if (geochron$type[error.na[i]] %in% c('Annual laminations (varves)')) {
        geochron$error[error.na[i]] = 10
      } else if (geochron$type[error.na[i]] %in% c('Biostratigraphic, pollen', 'Tsuga decline')) {
        geochron$error[error.na[i]] = 100
      } else if (geochron$type[error.na[i]] %in% c('Ambrosia rise', 'European settlement horizon')) {
        geochron$error[error.na[i]] = 50
      } else if (geochron$type[error.na[i]] %in% c('Tephra')) {
        geochron$error[error.na[i]] = 332
      } else if (geochron$type[error.na[i]] %in% c('Lead-210')) {
        geochron$error[error.na[i]] = 40 # fix this
      } else if (geochron$type[error.na[i]] %in% c('Other dating methods')) {
        geochron$error[error.na[i]] = 350
      } else {
        stop(paste0('an na chron control error could not be assigned for ', geochron$type[error.na[i]]))
      }
    }
  }
  
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
  
  #geochron$error[which(geochron$error == 0)] = 1

  out = Bchronology(ages=geochron$age,
                     ageSds=geochron$error, 
                     calCurves=calCurves,
                     positions=geochron$depth, 
                     positionThicknesses=rep(4,nrow(geochron)),
                     ids=geochron$chroncontrolid, 
                     predictPositions=depths,
                     jitterPositions=TRUE)

  summary(out, type='convergence', na.rm=TRUE)
  
  fname = paste0('Cores/', core.id, '/', core.id, '_bchron.pdf')
  p <-  plot(out, ageScale = c('bp')) +
    labs(title =  paste0(geochron$sitename[1], '; Dataset ID ', core.id),
         x = 'Age (cal years BP)',
         y = 'Depth (cm)')
  ggsave(fname, plot=p)
    
  post = data.frame(depths=depths, t(out$thetaPredict))
  write.table(post, paste0('.', '/Cores/', core.id, '/', 
                           core.id, '_bchron_samples.csv'), sep=',', col.names = TRUE, row.names = FALSE)
  
  predict_geo = predict(out, 
                        newPositions = geochron$depth)
  
  post = data.frame(labid=geochron$chroncontrolid, depths=geochron$depth, t(predict_geo))
  write.table(post, paste0('.', '/Cores/', core.id, '/', 
                           core.id, '_bchron_geo_samples.csv'), sep=',', col.names = TRUE, row.names = FALSE)
  
  if (length(bio.flag) >= 1){
    stop(paste0('unreliable controls used: ', paste(bio.flag, collapse='; ' )))
  }

  if (!is.na(age.flag)) {
    stop(age.flag)
  }

  TRUE
}

# debugging
#do_core_bchron(981, chron.control.types, mod)
#stop("debug stop")

mc.cores = 4 # processors

core.ids = list.files('Cores')
ncores = length(core.ids)

bchron.reports = mclapply(core.ids, function(core.id) {
  message(core.id)
  bchron.report = data.frame(datasetid = core.id, sitename=get_sitename(core.id), success = 1, reason=NA)
  success = try(do_core_bchron(core.id, chron.control.types, mod))
  if(is(success, "try-error")) {
    bchron.report$success = 0
    bchron.report$reason = geterrmessage()
  }
  bchron.report
}, mc.cores=mc.cores)

bchron.report = do.call(rbind, bchron.reports)
write.csv(bchron.report, paste0('bchron_report_v', version, '.csv'), row.names=FALSE)

fnames = list.files('Cores', '*_bchron.pdf', recursive=TRUE)

fname_str = sapply(fnames, function(x) paste0('Cores/', x))
fname_str = paste(fname_str, collapse = ' ')

sys_str = paste0("gs -sDEVICE=pdfwrite -o bchron_plots_v", version, ".pdf ", fname_str)
system(sys_str)
