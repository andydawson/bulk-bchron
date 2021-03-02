#  Note, put this code in a directory with a Bacon.R file and a Cores directory.
require(raster)
require(fields)
require(sp)
require(Bchron)
library(ggplot2)
library(parallel)
# 
# source('R/config.r')
# source('R/utils/helpers.r')

options(show.error.messages = TRUE)
#options(show.error.messages = FALSE)

version='1.0'

get_sitename <- function(core.id) {
  geochron <- try(read.table(paste0('Cores/', core.id, '/', core.id, '.csv'), sep=',', header=TRUE))
  if(is(geochron, "try-error")) {
    return(NA)
  }
  geochron$sitename[1]
}

do_core_bchron <- function(core.id) {

  geochron = read.table(paste0('Cores/', core.id, '/', core.id, '.csv'), sep=',', header=TRUE)
  depths   = scan(paste0('Cores/', core.id, '/', core.id, '_depths.txt'))
  
  stopifnot(nrow(geochron) >= 2)
  
  calCurves = rep(NA, nrow(geochron))
  calCurves[which(geochron$cc == 0)] = "normal"
  calCurves[which(geochron$cc == 1)] = "intcal20"
  
  geochron$error[which(geochron$error == 0)] = 1
  
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

  TRUE
}





mc.cores = 4 # processors

core.ids = list.files('Cores')
ncores = length(core.ids)

bchron.reports = mclapply(core.ids, function(core.id) {
  message(core.id)
  bchron.report = data.frame(datasetid = core.id, sitename=get_sitename(core.id), success = 1)
  success = try(do_core_bchron(core.id))
  if(is(success, "try-error")) {
    bchron.report$success = 0
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



