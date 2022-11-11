library(maps)
library(mapdata)
library(ggplot2)

# Maps #
dat =  read.csv('data/chroncontrol_summary_pollen_full.csv', stringsAsFactors = FALSE)
dat =  dat[which(dat$longitude<100),]

bacon_dat = read.csv('wang/SiteInfo_fullcore.csv')

foo = read.csv('bchron_report_v9.0.csv')
foo = foo[which(foo$success == 1),]

match(foo$datasetid, dat$datasetid)
dat[match(foo$datasetid, dat$datasetid), 'latitude']

foo$latitude = dat[match(foo$datasetid, dat$datasetid), 'latitude']

foo$longitude = dat[match(foo$datasetid, dat$datasetid), 'longitude']



colors = c("All Sites" = "blue", "Sites With Reliable Controls" = "green", "Bacon Sites" = "red")

#  #
world = map_data('world')
na = world[which(world$region %in% c('Canada', 'USA', 'Mexico')),]
na = na[which((na$long < (-50))&(na$long >(-200))),]
ggplot() +
  geom_polygon(data= na, aes(x=long, y=lat, group=group)) +
  geom_point(data = dat, aes(x =longitude, y =latitude, color = "All Sites"), size = 0.5) +
  geom_point(data = foo, aes(x =longitude, y =latitude, color = "Sites With Reliable Controls"), size = 0.5) +
  geom_point(data = bacon_dat, aes(x =lon, y =lat, color = "Bacon Sites"), size = 0.5) +
  labs (title = "Map of Pollen Core Collection Sites", hjust=0, 
        caption = "A map of North America showing each site where a pollen core was taken") +
  theme(plot.caption = element_text(hjust=0, face="italic", size=7),
        plot.title = element_text(hjust=0.5)) +
  xlab ("Longitude") +
  ylab ("Latitude")

ggsave('figures/full_map.png')

world = map_data('world')
na = world[which(world$region %in% c('Canada', 'USA', 'Mexico')),]
na = na[which((na$long < (-50))&(na$long >(-200))),]
ggplot() +
  geom_polygon(data= na, aes(x=long, y=lat, group=group)) +
  geom_point(data = foo, aes(x =longitude, y =latitude), color = "Dodger Blue 1", size = 0.5) +
  labs (title = "Map of Pollen Core Collection Sites", hjust=0, 
        caption = "A map of North America showing each site where a pollen core was taken") +
  theme(plot.caption = element_text(hjust=0, face="italic", size=7),
        plot.title = element_text(hjust=0.5)) +
  xlab ("Longitude") +
  ylab ("Latitude")

ggsave('figures/bchron_only_map.png')

# Just Bacon sites #

world = map_data('world')
na = world[which(world$region %in% c('Canada', 'USA', 'Mexico')),]
na = na[which((na$long < (-50))&(na$long >(-200))),]
ggplot() +
  geom_polygon(data= na, aes(x=long, y=lat, group=group)) +
  geom_point(data = bacon_dat, aes(x =lon, y =lat), color = "red", size = 0.5) +
  labs (title = "Map of Bacon Pollen Core Collection Sites", hjust=0, 
        caption = "A map of North America showing each site where both Bchron and Bacon were used") +
  theme(plot.caption = element_text(hjust=0, face="italic", size=7),
        plot.title = element_text(hjust=0.5)) +
  xlab ("Longitude") +
  ylab ("Latitude")

ggsave('figures/bacon_only_map.png')

# ### Tables ### #
# Frequency of each control type #
chron_control_types <- read.csv("chroncontrol_types-edited.csv")
all_sites = read.csv("data/chroncontrol_summary_pollen_full.csv")

control_freq = table(all_sites$type)


# frequency of sites for each method #
# bchron
length(unique(dat$datasetid))
# reliable bchron
length()
# bacon
length(unique(bacon_dat$siteID))




