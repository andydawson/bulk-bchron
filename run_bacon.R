library(rbacon)
library(dplyr)
library(readr)

#' @title Get age posteriors
#' @description Using the output files from Bacon get the full posterior at depths.
#' @param handle The site handle.
#' 
# samples[j,] = bacon_geo_posts(d=geochron$depth[j], b.depths=info$depths, out=outputs, thick=site_params$bestthick)

# Bacon.Age.d(geochron$depth[j], set=info, its=output, BCAD=set$BCAD)
Bacon.Age.d <- function(d, set=info, its=info$output, BCAD=set$BCAD)
{ 
  # intercept
  elbows <- cbind(its[,1])
  
  # accumulation rates
  if ((ncol(its)-2) > 1){
    accs <- its[,2:(ncol(its)-1)]
    # accs <- its[,1:(ncol(its)-1)]
  } else {
    accs <- as.matrix(its[,2:(ncol(its)-1)], ncol=1)
  }
  
  for(i in 1:ncol(accs))
    elbows <- cbind(elbows, elbows[,ncol(elbows)] + (set$thick * accs[,i]))
  
  # if(d %in% set$depths){
  #   ages <- elbows[,which(set$depths == d)] 
  # } else {
  #   maxd <- max(which(set$depths < d))
  #   ages <- elbows[,maxd] + ((d-set$depths[maxd]) * accs[,maxd])
  # }
  # if(!is.na(set$hiatus.depths)[1])
  #   for(hi in set$hiatus.depths)
  #     #        {
  #     #          below <- min(which(set$d > hi), set$K-1)+1
  #     #          above <- max(which(set$d < hi))
  #     #          if(d > set$d[above] && d < set$d[below])
  #     #            {
  #     #              start.hiatus <- elbows[,below] - (its[,1+below] * (set$d[below] - hi))
  #     #              end.hiatus <- elbows[,above] + (its[,above] * (hi - set$d[above]))
  #     #              ok <- which(end.hiatus < start.hiatus)
  #     #              if(d < hi)
  #     #                ages[ok] <- elbows[ok,above] + (its[ok,above] * (d - set$d[above])) else
  #     #                ages[ok] <- elbows[ok,below] - (its[ok,1+below] * (set$d[below] - d))
  # #            }
  # #        }
  #   {
  #     above <- min(which(set$depths > hi))
  #     below <- max(which(set$depths < hi))
  #     
  #     if (d > set$depths[below] && d < set$depths[above])
  #     {
  #       start.hiatus <- elbows[,below] + (accs[,below] * (hi - set$depths[below]))
  #       end.hiatus   <- elbows[,above] - (accs[,above] * (set$depths[above] - hi))
  #       ok <- which(end.hiatus > start.hiatus)
  #       if(d < hi){
  #         ages[ok] <- elbows[ok,below] + (accs[ok,below] * (d - set$depths[below])) 
  #       } else {
  #         ages[ok] <- elbows[ok,above] - (accs[ok,above] * (set$depths[above] - d))
  #       }
  #     }
  #   }
  if(d %in% set$elbows){
    ages <- elbows[,which(set$elbows == d)] 
  } else {
    maxd <- max(which(set$elbows < d))
    ages <- elbows[,maxd] + ((d-set$elbows[maxd]) * accs[,maxd])
  }
  if(!is.na(set$hiatus.depths)[1])
    for(hi in set$hiatus.depths)
      #        {
      #          below <- min(which(set$d > hi), set$K-1)+1
      #          above <- max(which(set$d < hi))
      #          if(d > set$d[above] && d < set$d[below])
      #            {
      #              start.hiatus <- elbows[,below] - (its[,1+below] * (set$d[below] - hi))
      #              end.hiatus <- elbows[,above] + (its[,above] * (hi - set$d[above]))
      #              ok <- which(end.hiatus < start.hiatus)
      #              if(d < hi)
      #                ages[ok] <- elbows[ok,above] + (its[ok,above] * (d - set$d[above])) else
      #                ages[ok] <- elbows[ok,below] - (its[ok,1+below] * (set$d[below] - d))
  #            }
  #        }
    {
      above <- min(which(set$elbows > hi))
      below <- max(which(set$elbows < hi))
      
      if (d > set$elbows[below] && d < set$elbows[above])
      {
        start.hiatus <- elbows[,below] + (accs[,below] * (hi - set$elbows[below]))
        end.hiatus   <- elbows[,above] - (accs[,above] * (set$elbows[above] - hi))
        ok <- which(end.hiatus > start.hiatus)
        if(d < hi){
          ages[ok] <- elbows[ok,below] + (accs[ok,below] * (d - set$elbows[below])) 
        } else {
          ages[ok] <- elbows[ok,above] - (accs[ok,above] * (set$elbows[above] - d))
        }
      }
    }
  if(BCAD) ages <- 1950 - ages
  ages
}

bacon_age_posts <- function(site_params, core_path) {
  
  handle = site_params$handle
  
  thick = site_params$bestthick
  
  depth_file <- paste0(core_path, "/",
                       handle, "/", handle, "_depths.txt")
  
  settings_file <- paste0(core_path, "/",
                          handle, "/", handle, "_settings.txt")
  
  out_files <- list.files(paste0(core_path, "/", handle),
                          pattern = ".out$",
                          full.names = TRUE)
  
  # assertthat::assert_that(length(out_files) > 0,
  #                         msg = list.files(paste0(settings$core_path, "/", handle)))
  
  depth <- suppressMessages(readr::read_csv(depth_file,
                                            col_names = FALSE)) %>%
    as.data.frame()
  
  bacon_settings <- suppressMessages(readr::read_csv(settings_file,
                                                     col_names = FALSE, comment = "#")) %>%
    as.data.frame()
  
  if (length(out_files) > 1) {
    message("Multiple Bacon output files exist.")
  }
  
  for (k in 1:length(out_files)) {
    
    outer <- suppressMessages(readr::read_delim(out_files[k],
                                                col_names = FALSE, delim = " ")) %>%
      as.data.frame
    
    outer = outer[sample(nrow(outer), 1000),]
    
    # We can do this match because we know how Bacon writes out files.
    sections <- stringr::str_match(out_files[k], "(?:_)([0-9]*)\\.")[2] %>%
      as.numeric()
    
    if (ncol(outer) == (sections + 3)) {
      posteriors <- matrix(NA, nrow = nrow(depth), ncol = nrow(outer))
      
      for (j in 1:nrow(outer)) {
        # x <- seq(from = bacon_settings[1, 1],
        #          to = bacon_settings[2, 1],
        #          # length.out = sections)
        #          length.out = sections + 1)
        x <- seq(from = bacon_settings[1, 1],
                 by = thick,
                 # length.out = sections)
                 length.out = sections)
        y <- c(outer[j, 1],
               outer[j, 1] +
                 cumsum( ( diff(x) * outer[j, 2:(ncol(outer) - 2)]) %>%
                           as.numeric))
        posteriors[, j] <- bacon_extrap(x,
                                        y = y,
                                        xout = depth %>% unlist())
      }
      
      posterior_file <- paste0(core_path, "/",
                               handle, "/", handle, "_", sections, "_posteriorout.csv")
      
      readr::write_csv(posterior_file,
                       x = as.data.frame(posteriors))
    }
    
  }
  if (length(out_files) == 1) {
    return(as.data.frame(posteriors))
  }
  
}


#' @title Extrapolation for Bacon posterior estimation.
#' @description When sample depths are below the last chronological control in a Bacon 
#' record we get NAs in the interpolation.  This function attempts to manage that issue.
#' It gets called in the \code{bacon_age_posts()} function.
#' @return A numeric vector.
bacon_extrap <- function(x, y, xout) {
  
  out_depth <- xout %>% unlist
  
  if (max(out_depth) > max(x)) {
    
    slope <- diff(tail(y, n = 2)) / diff(tail(x, n=2))
    x <- c(x, max(out_depth))
    y <- c(y, tail(y, n = 1) + tail(diff(x), n = 1) * slope)
  }
  outputs <- round(approx(x = x, y = y,
                          xout = xout %>% unlist())$y, 0)
  
  return(outputs)
  
}


settings = read.csv('Cores_bacon/SiteInfo_fullcore.csv', stringsAsFactors = FALSE)

ids = unique(settings$datasetid)
N_ids = length(ids)


for (i in 1:N_ids){
  print(i)
  
  site_params = settings[i, ]
  
  id = site_params$datasetid[1]
  
  core_path = paste0('Cores_bacon/Cores_full/')
  bacon_chrons <- paste0('Cores_bacon/Cores_full/', site_params$handle, '/', site_params$handle, '.csv')
  bacon_depths <- paste0('Cores_bacon/Cores_full/', site_params$handle, '/', site_params$handle, '_depths.txt')
  
  # find hiatus depth
  geochron <- suppressMessages(readr::read_csv(bacon_chrons))
  depths<- scan(bacon_depths)
  
  if (any(diff(geochron$depth) < 0)) {
    geochron <- geochron %>% arrange(depth)
    readr::write_csv(geochron, bacon_chrons)
  }
  
  gcol <- which(tolower(colnames(geochron)) == 'labid')
  
  sett_layer <- stringr::str_detect(unlist(geochron[,gcol]), "sett")
  
  if (any(sett_layer) & nrow(geochron) > 2) {
    
    # determine which bacon parameters to input if preset is the last sample
    if (which(sett_layer) == nrow(geochron)) {
      hiatus.depth <- NA
      acc.mean.val <- site_params$acc.mean.mod
      acc.shape.val <- site_params$acc.shape.mod
      site_params$hiatus <- 0
    } else if (which(sett_layer) == 1) {
      # if preset is the first sample
      hiatus.depth <- NA
      acc.mean.val <- site_params$bestacc
      acc.shape.val <- site_params$acc.shape.old
      site_params$hiatus <- 0
    } else {
      hiatus.depth <- geochron$depth[sett_layer]
      acc.mean.val <- c(site_params$acc.mean.mod,
                        site_params$bestacc)
      acc.shape.val <- c(site_params$acc.shape.mod,
                         site_params$acc.shape.old)
      site_params$hiatus <- 1
    }
    
  } else if (any(sett_layer) & nrow(geochron) == 2) {
    # if preset and only two geochron samples, use modern priors
    hiatus.depth <- NA
    acc.mean.val <- site_params$acc.mean.mod
    acc.shape.val <- site_params$acc.shape.mod
    site_params$hiatus <- 0
  } else if (!any(sett_layer)) {
    # if no preset then use historical priors
    hiatus.depth <- NA
    acc.mean.val <- site_params$bestacc
    acc.shape.val <- site_params$acc.shape.old
    site_params$hiatus <- 0
  }
  
  
  bestthick = as.numeric(site_params$bestthick)
  
  acc.mean = acc.mean.val
  
  mem.strength = 10
  
  if (id %in% c(13029, 15750, 2598)) {
    mem.strength = 10
  } else {
    mem.strength = 10
  }
  
  if (id %in% c(14539, 14645)){
    acc.mean = c(3.02, 30)
    bestthick = 5
  } else if (id %in% c(14991)){
    acc.mean = c(3.02, 40)
    bestthick = 5
  } else if (id %in% c(15035)){
    acc.mean = 10
  } else if  (id %in% c(15138)){
    acc.mean = c(3.02, 50)
    mem.strength = 10
    bestthick = 50
  } else if (id %in% c(15185)){
    acc.mean = 10
  } else if (id %in% c(15197)){
    acc.mean = 10
  } else if (id %in% c(1523)){
    acc.mean = 10
    bestthick = 60
  }  else if (id %in% c(15302)){
    acc.mean = 10
    bestthick = 40
    # bestthick = 60
  }  else if (id %in% c(15738)){
    acc.mean = 10
    bestthick = 40
    # bestthick = 60
  }
  
  # if (id %in% c(15750, 1136)) {
  #   acc.mean = 5
  # } else if (id %in% c(1144)) {
  #   acc.mean = 100
  # } else if (id %in% c(12001)) {
  #   acc.mean = 50
  # } else {
  #   acc.mean = 20
  # }


  acc.shape.val
  as.numeric(site_params$mem.strength)
  as.numeric(site_params$mem.mean)
 # 
  
  
 # source('Bacon.R')
  out <- try(Bacon(core = site_params$handle,
                   coredir = core_path,
                   acc.mean = acc.mean,
                   # acc.shape = as.numeric(acc.shape.val),
                   mem.strength = mem.strength,
                   # mem.mean = as.numeric(site_params$mem.mean),
                   thick = bestthick,#as.numeric(site_params$bestthick),
                   ask = FALSE,
                   suggest = FALSE,
                   depths.file = TRUE,
                   postbomb = 1,
                   hiatus.max = 10,
                   hiatus.depths = hiatus.depth,
                   ssize = 5000,
                   burnin = 1000,
                   plot.pdf=TRUE))
  
  # agedepth(info)
  # iters = 1000
  

  
  write.table(info$elbows, paste0('Cores_bacon/Cores_full/', site_params$handle, '/', site_params$handle, '_bacon_depths.txt'), row.names = FALSE, col.names = FALSE)
  
  outfile = paste0(core_path,  site_params$handle, '/', site_params$handle, '_', length(info$elbows), '.out')
  output = read.table(outfile)
  
  
  # outputs <- bacon_age_posts(site_params$handle, core_path)
  
  iters = nrow(info$output)
  posteriors <- matrix(NA, nrow = length(depths), ncol = iters)
  
  for (d in 1:length(depths)){
    posteriors[d,] = Bacon.Age.d(depths[d], set=info, its=output, BCAD=info$BCAD)
  }
  
  posterior_file <- paste0(core_path, "/",
                           site_params$handle, "/", site_params$handle, "_", length(info$elbows), "_posteriorout.csv")
  
  # readr::write_csv(posterior_file,
  #                  x = as.data.frame(posteriors))
  
  write.table(posteriors, posterior_file, sep=',', col.names = TRUE, row.names = FALSE)
    
  # 
  # out_files <- list.files(paste0(core_path, "/", handle),
  #                         pattern = ".out$",
  #                         full.names = TRUE)
  # 
  # assertthat::assert_that(length(out_files) > 0,
  #                         msg = list.files(paste0(settings$core_path, "/", handle)))
  # 
  # outer <- suppressMessages(readr::read_delim(out_files[k],
  #                                             col_names = FALSE, delim = " ")) %>%
  #   as.data.frame
  
  samples = matrix(0, nrow = length(geochron$depth), ncol = nrow(output))
  colnames(samples) = paste0('iter', rep(1:nrow(output)))
  for (j in 1:length(geochron$depth)){    
    # print(j)
    samples[j,] = Bacon.Age.d(geochron$depth[j], set=info, its=output, BCAD=info$BCAD)
    # samples[j,] = bacon_geo_posts(d=geochron$depth[j], b.depths=depths, out=outer, thick=site_params$bestthick)
  }
  
  # samples = matrix(0, nrow = length(geochron$depth), ncol = nrow(info$output))
  # colnames(samples) = paste0('iter', rep(1:nrow(info$output)))
  # for (j in 1:length(geochron$depth)){    
  #   # print(j)
  #   samples[j,] = Bacon.Age.d(geochron$depth[j], set=info, its=info$output, BCAD=set$BCAD)
  #   # samples[j,] = bacon_geo_posts(d=geochron$depth[j], b.depths=depths, out=outer, thick=site_params$bestthick)
  # }
  
  post = data.frame(depths=geochron$depth, samples)
  
  write.table(post, paste0('Cores_bacon/Cores_full/', site_params$handle, "/", 
                           site_params$handle, "_geo_samples.csv"), sep=',', col.names = TRUE, row.names = FALSE)
  
  
  depths_all = sort(c(depths, geochron$depth))
  depth_min  = min(depths_all)
  depth_max = max(depths_all)
  depths_span = seq(depth_min, depth_max, length=50)
  
  samples_span = matrix(0, nrow = length(depths_span), ncol = nrow(output))
  colnames(samples_span) = paste0('iter', rep(1:nrow(output)))
  for (j in 1:length(depths_span)){    
    # print(j)
    samples_span[j,] = Bacon.Age.d(depths_span[j], set=info, its=output, BCAD=info$BCAD)
    # samples[j,] = bacon_geo_posts(d=geochron$depth[j], b.depths=depths, out=outer, thick=site_params$bestthick)
  }
  
  # samples = matrix(0, nrow = length(geochron$depth), ncol = nrow(info$output))
  # colnames(samples) = paste0('iter', rep(1:nrow(info$output)))
  # for (j in 1:length(geochron$depth)){    
  #   # print(j)
  #   samples[j,] = Bacon.Age.d(geochron$depth[j], set=info, its=info$output, BCAD=set$BCAD)
  #   # samples[j,] = bacon_geo_posts(d=geochron$depth[j], b.depths=depths, out=outer, thick=site_params$bestthick)
  # }
  
  post_span = data.frame(depths=depths_span, samples_span)
  
  write.table(post_span, paste0('Cores_bacon/Cores_full/', site_params$handle, "/", 
                           site_params$handle, "_span_samples.csv"), sep=',', col.names = TRUE, row.names = FALSE)
  
}







