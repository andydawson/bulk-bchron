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

chron_control_types <- read.csv("data/chroncontrol_types-edited.csv")

radio = read.csv('data/radiocarbon-dates-errors.csv')
mod_radio <- gam(error ~ s(age, k=15), data=radio, method='REML', family=Gamma(link="identity"))

lead = read.csv('data/lead-dates-errors.csv')
mod_lead <- gam(error ~ s(age, k=4), data=lead, method='REML', family=Gamma(link="log"))
# 
# extrap = 1000

options(show.error.messages = TRUE)
#options(show.error.messages = FALSE)

version='10.0'

get_sitename <- function(core.id) {
  geochron <- try(read.table(paste0('Cores/', core.id, '/', core.id, '.csv'), sep=',', header=TRUE))
  if(is(geochron, "try-error")) {
    return(NA)
  }
  geochron$sitename[1]
}


fix_geochron_errors <- function(geochron){
  
  # fix some errors in data
  if (any(is.na(geochron$chroncontrolid))){
    geochron$chroncontrolid[is.na(geochron$chroncontrolid)] = -999
  }
  if (any(is.na(geochron$type))){
    geochron$type[is.na(geochron$type)] = -999
  }
  if (geochron$datasetid[1] == "774"){
    geochron$age = 1950 - geochron$age
  }
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
  
  if ((all(!(is.na(geochron$limityounger))))&(any(geochron$limityounger < (-70)))){
    geochron$limityounger[which(geochron$limityounger < (-70))] = -70
  }
  if (any((geochron$type == 'Lead-210') & (!is.na(geochron$limitolder)) & (geochron$limitolder > 200))){
    geochron$limitolder[which((geochron$type == 'Lead-210') & (geochron$limitolder > 200))] = 200
  }
  if (any(is.na(geochron$age))){
    geochron$age[which(is.na(geochron$age))] = (geochron$limityounger[which(is.na(geochron$age))] + geochron$limitolder[which(is.na(geochron$age))])/2
  }
  
  return(geochron)
}

set_bio_flag <- function(chron_control_types, keep, idx){
  
  bio.flag = c()
  if (any(keep != 1)) {
    for (i in 1:length(idx)){
      if (chron_control_types$chron.control.type[idx[i]] == "Biostratigraphic, pollen") {
        keep[i] = 1 
        bio.flag = c(bio.flag, "Biostratigraphic, pollen")
      } else if (chron_control_types$chron.control.type[idx[i]] == "Ambrosia rise") {
        keep[i] = 1 
        bio.flag = c(bio.flag, "Ambrosia rise")
      } else if (chron_control_types$chron.control.type[idx[i]] == "European settlement horizon") {
        keep[i] = 1 
        bio.flag = c(bio.flag, "European settlement horizon, Pre-EuroAmerican settlement horizon")
      } else if (chron_control_types$chron.control.type[idx[i]] == "Tsuga decline") {
        keep[i] = 1 
        bio.flag = c(bio.flag, "Tsuga decline")
      } else if (chron_control_types$chron.control.type[idx[i]] == "Casuarina rise"){
        keep[i] = 1 
        bio.flag = c(bio.flag, "Casuarina rise")
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
      } else if (geochron$type[error.na[i]] %in% c('Core top', 'Collection date', 'Core top, estimated')) {
        geochron$error[error.na[i]] = 32#.27
        #} else if (geochron$type[error.na[i]] %in% c('Radiocarbon, average of two or more dates', 'Radiocarbon, reservoir correction')){
        #  geochron$error[error.na[i]] = 50
        # >>>>>>> 904b764fa902585d4b72f3c53a83ce22e4271314
      } else if (geochron$type[error.na[i]] %in% c('Annual laminations (varves)')) {
        geochron$error[error.na[i]] = (geochron$age[error.na[i]]+70)*0.05
      } else if (geochron$type[error.na[i]] %in% c('Biostratigraphic, pollen', 'Tsuga decline', 'Casuarina rise', 'Salsola rise', 'Biostratigraphic')) {
        geochron$error[error.na[i]] = 250
      } else if (geochron$type[error.na[i]] %in% c('Ambrosia rise', 'European settlement horizon', 'Historical fire')) {
        geochron$error[error.na[i]] = 50
      } else if (geochron$type[error.na[i]] %in% c('Tephra')) {
        geochron$error[error.na[i]] = 334
      } else if (geochron$type[error.na[i]] %in% c('Lead-210')) {
        geochron$error[error.na[i]] = round(mean(predict.gam(mod_lead, new_age=data.frame(age=geochron$radio[error.na[i]]), type='response', se.fit=TRUE)$fit))
      } else if (geochron$type[error.na[i]] %in% c('Other dating methods')) {
        geochron$error[error.na[i]] = 350
      } else {
        # stop(paste0('an na chron control error could not be assigned for ', geochron$type[error.na[i]]))
        print(paste0('an na chron control error could not be assigned for ', geochron$type[error.na[i]]))
      }
    }
  }
  
  return(geochron)
  
}

prepare_geochron <- function(core.id, chron_control_types, mod_radio, mod_lead, extrap) {
  
  # read in the geochron table for core.id
  geochron = read.table(paste0('Cores/', core.id, '/', core.id, '.csv'), sep=',', header=TRUE)
  
  # remove duplicated rows (rare)
  geochron <- unique(geochron)
  
  geochron <- geochron[which(!is.na(geochron$depth)),]
  
  if (nrow(geochron) < 1) {
    stop("no controls")
  }  
  
  # STOP if only one control
  if (nrow(geochron) == 1) {
    # stop("only 1 control before filtering")
    geochron[, 'keep'] = 0
    print("only 1 control before filtering")
  }  
  
  # fix some specific errors
  geochron = fix_geochron_errors(geochron)
  
  
  # STOP if no controls
  if (nrow(geochron) < 1) {
    stop("no controls before filtering")
  }
  
  # match core.id control types to master list
  # STOP if any types not in master list
  idx <- match(geochron$type, chron_control_types$chron.control.type, nomatch=NA)
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
  geochron$cc = chron_control_types$cc[idx]
  geochron$keep  = chron_control_types$keep[idx]
  
  # if (is.na(keep)) {
  #   stop(paste0(nrow(geochron), " controls; nbut"))
  # }
  
  bio.flag = set_bio_flag(chron_control_types, geochron$keep, idx)
  
  age.flag = NA
  if (any(is.na(geochron$age))) {
    age.flag = "had one or more NA ages, but ran anyway"
  }
  
  # remove rows from geochron that have NA ages
  geochron = geochron[!is.na(geochron$age),]
  
  # discard unreliable controls
  # only keep rows with keep flag of 1
  
  idx_keep = which(geochron$keep==1)
  # geochron = geochron[which(keep==1),]
  # STOP if fewer than 2 ages left

  if (nrow(geochron[idx_keep,]) < 2) {
    geochron[, 'keep'] = 0
    # stop("less than 2 non-NA ages")
    print("less than 2 non-NA ages")
  }
  
  # only use site if 3 or more non-biostrat dates
  # or if biostrat date is ambrosia rise, use
  types = chron_control_types$chron.control.type[match(geochron[idx_keep,]$type, chron_control_types$chron.control.type, nomatch=NA)]
  if (any(types %in% c('Biostratigraphic, pollen', 'Tsuga decline'))){
    # if (sum(!(types %in% c('Biostratigraphic, pollen', 'Tsuga decline')))>= 3){
    #   stop("Biostrat control with fewer than 3 controls of other types")
    # }
    if (sum(!(types %in% c('Biostratigraphic, pollen', 'Tsuga decline'))) < 2){
      geochron[, 'keep'] = 0
      # stop("Biostrat control with fewer than 2 controls of other types")
      print("Biostrat control with fewer than 2 controls of other types")
    }
  }
  
  # for any NA or 0 control errors, set them to informed value
  # positive error required for Bchron
  geochron = set_missing_control_errors(geochron, mod_radio, mod_lead)
  
  if (nrow(geochron) <= 1) {
    geochron[, 'keep'] = 0
    # stop("one or no controls after filtering")
    print("one or no controls after filtering")
  }
  
  if (any(is.na(geochron$cc))) {
    geochron[, 'keep'] = 0
    # stop('one or more cc values are NA')
    print('one or more cc values are NA')
  }
  
  if (any((geochron$error == 0)|is.na(geochron$error))){
    geochron[, 'keep'] = 0
    # stop('geochron error still set to zero')
    print('geochron error still set to zero')
  }
  
  write.csv(geochron, paste0('Cores/', core.id, '/', core.id, '_prepared.csv'), row.names = FALSE) 
 
}

# debugging

core.ids = list.files('Cores')
ncores = length(core.ids)

for (i in 1:ncores) {

  core.id = core.ids[i]
  
  print(i)
  # print(core.id)
  
  success = try(prepare_geochron(core.id, chron_control_types, mod_radio, mod_lead, extrap))
}

