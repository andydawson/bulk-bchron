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

extrap = 1000

options(show.error.messages = TRUE)
#options(show.error.messages = FALSE)

chron_control_types <- read.csv("data/chroncontrol_types-edited.csv")

radio = read.csv('data/radiocarbon-dates-errors.csv')
mod_radio <- gam(error ~ s(age, k=15), data=radio, method='REML', family=Gamma(link="identity"))

lead = read.csv('data/lead-dates-errors.csv')
mod_lead <- gam(error ~ s(age, k=4), data=lead, method='REML', family=Gamma(link="log"))

version='10.0'

get_sitename <- function(core.id) {
  geochron <- try(read.table(paste0('Cores/', core.id, '/', core.id, '_prepared.csv'), sep=',', header=TRUE))
  if(is(geochron, "try-error")) {
    return(NA)
  }
  geochron$sitename[1]
}


# do_core_bchron(core.id, chron_control_types, mod_radio, mod_lead, extrap)

do_core_bchron <- function(core.id, chron_control_types, mod_radio, mod_lead, extrap) {
  
  # read in the geochron table for core.id
  geochron = read.table(paste0('Cores/', core.id, '/', core.id, '_prepared.csv'), sep=',', header=TRUE)
  
  # remove duplicated rows (rare)
  geochron <- unique(geochron)
  
  geochron <- geochron[which(!is.na(geochron$depth)),]
  
  # STOP if only one control
  if (nrow(geochron) <= 1) {
    stop("only 1 control before filtering")
  }  

  
  # STOP if no controls
  if (nrow(geochron) < 1) {
    stop("no controls before filtering")
  }
  

  # discard unreliable controls
  # only keep rows with keep flag of 1
  geochron = geochron[which(geochron$keep==1),]
  
  # remove rows from geochron that have NA ages
  # STOP if fewer than 2 ages left
  geochron = geochron[!is.na(geochron$age),]
  if (nrow(geochron) < 2) {
    stop("less than 2 non-NA ages")
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
  geochron = geochron[order(geochron$depth),]
  
  out = Bchronology(ages    = geochron$age,
                    ageSds = geochron$error, 
                    calCurves = calCurves,
                    positions = geochron$depth, 
                    positionThicknesses = rep(4,nrow(geochron)),
                    ids = geochron$chroncontrolid, 
                    predictPositions = depths,
                    allowOutside = TRUE)
  
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
  
  TRUE
}

# debugging
#do_core_bchron(1000, chron_control_types, mod)
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
#   success = try(do_core_bchron(core.id, chron_control_types, mod_radio, mod_lead, extrap))
#   if(is(success, "try-error")) {
#     bchron.report$success = 0
#     bchron.report$reason = geterrmessage()
#   }
#   bchron.report
# }, mc.cores=mc.cores)

bchron.report = data.frame(datasetid = numeric(0), sitename=character(0), success = numeric(0), reason=character(0))

# do in chunks cause busted
for (i in 1023:ncores) {
  
  # check these
  # if (i==1429){next}
  if (i==801){next} 
  if (i==813){next}   
  if (i==823){next} 
  if (i==852){next}
  if (i==853){next} 
  if (i==1022){next} 
  
  core.id = core.ids[i]
  
  print(i)
  print(core.id)
  
  bchron.report.site = data.frame(datasetid = core.id, sitename=get_sitename(core.id), success = 1, reason=NA)
  success = try(do_core_bchron(core.id, chron_control_types, mod_radio, mod_lead, extrap))
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

