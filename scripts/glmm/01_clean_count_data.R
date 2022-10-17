################################################################################
################# West Coast Env Health GLMM Analysis ##########################
#################      Step 1: Clean Count Data       ##########################
################################################################################


##################### 
# install.packages("_")
library(data.table)
library(sf)
library(sp)
library(igraph)
library(dplyr)
library(here)
library(readr)
library(terra)

# read in data
all_data <- read_csv("data/raw_data_from_uwin/yasmine_west_coast_data.csv")
length(unique(all_data$locationAbbr))
class(all_data)

# check data
head(all_data)
nrow(all_data)
unique(all_data$commonName)

# read in detection history data set
detections <- read_csv("data/raw_data_from_uwin/initial_data_yasmine.csv")

# remove all sites with 0 days running from detection data
detections <- detections[detections$J != 0, ]
(unique(detections$Site))

# 77 total sites with data 

# collapse data so each instance of city/location/season/species is a new column

counts <- all_data %>% group_by(commonName, locationID, city, season) %>% summarize(count = n())
nrow(counts)
nrow(all_data)

colnames(all_data)
colnames(counts)

# take site info
all_data <- all_data %>%
  dplyr::select(locationID, locationAbbr, utmEast, utmNorth, utmZone) %>% distinct()
head(all_data)

# join counts to site data 
joined_data <- left_join(counts, all_data, by = "locationID")
nrow(joined_data)
length(unique(joined_data$locationAbbr))

joined_data
# now we have counts for each site/season combo 

# let's change UTM to lat long


# first make all of our zones just 10 or 11
joined_data$utmZone[joined_data$utmZone == "11S"] <- 11
joined_data$utmZone[joined_data$utmZone == "10N"] <- 10
joined_data$utmZone[joined_data$utmZone == "11"] <- 11
joined_data$utmZone[joined_data$utmZone == "11N"] <- 11
unique(joined_data$utmZone)

# split data based on utm zone
data_10 <- subset(joined_data, utmZone == "10")
data_11 <- subset(joined_data, utmZone == "11")

# get lat long for all points in zone 10
utm1 <- data.frame(x=data_10$utmEast, y=data_10$utmNorth) 
coordinates(utm1) <- ~x+y 
class(utm1)
proj4string(utm1) <- CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +ellps=WGS84") 
utm2 <- spTransform(utm1,CRS("+proj=longlat +datum=WGS84"))
utm2

data_10$long <- utm2$x
data_10$lat <- utm2$y

# get lat long for all points in zone 11
utm1 <- data.frame(x=data_11$utmEast, y=data_11$utmNorth) 
coordinates(utm1) <- ~x+y 
class(utm1)
proj4string(utm1) <- CRS("+proj=utm +zone=11 +datum=WGS84 +units=m +ellps=WGS84") 
utm2 <- spTransform(utm1,CRS("+proj=longlat +datum=WGS84"))
utm2

data_11$long <- utm2$x
data_11$lat <- utm2$y

# bind back into 1 data set
new_dat <- rbind(data_10,data_11)
head(new_dat)

# change all city into lowercase
new_dat$city <- tolower(new_dat$city)

# change colnames
colnames(new_dat) <- c("Species", "locationID", "city", "season",
                       "count", "Site", "utmEast", "utmNorth",
                       "utmZone", "long", "lat")

################################################################################
## removing/merging sites 

# use mason's functions to collapse sites by location 
source("scripts/mason_site_collapse_code/qaqc_sites.R")
source("scripts/mason_site_collapse_code/long_to_zone.R")
source("scripts/mason_site_collapse_code/fix_site_names.R")

new_dat_1 <- qaqc_sites(x = new_dat, cities="city", sites = "Site",
                        my_coords=c("long","lat"), my_crs=4326)

# this function will tell us what sites need to be merged or removed based
# on proximity to the next site(s)
# the "remove" may not be totally accurate for this data set because it is 
# made to keep the site with more data in an occupancy data set, but 
# we have counts. for now we will keep these separate 

fix_site_names(new_dat_1)

# we will also remove sites that don't have enough data in detection data
# the detection data only includes sites that have been fully tagged for a season

# first see which sites are not present in count data
sites_det_only <- setdiff(unique(detections$Site),unique(new_dat$Site))

# now see which sites are not present in detection data 
sites_count_only <- setdiff(unique(new_dat$Site),unique(detections$Site))

# manually decide which sites to remove from count data after comparing 

# merge  H01-SPS1 and H01-SPS2; H01-HNP1 and H01-HNP2; H01-CFR2 and H01-CFR1
# this is based on the fix_site_names function after cross-checking detections

new_dat$Site[new_dat$Site == "H01-SPS2"] <- "H01-SPS1" # keep
new_dat$Site[new_dat$Site == "H01-HNP2"] <- "H01-HNP1" # keep

# remove all sites that don't exist in detection dataset 

new_dat_cut <- subset(new_dat, Site!="KP" & Site!="H9" & Site != "P-MCP"
                      & Site != "H01-ARB1" & Site != "H01-GAT4"
                      & Site != "H01-CFR2" & Site != "H01-AST1" & Site != "H01-LTC1"
                      & Site != "ED08" & Site != "PCNYN02" & Site != "H01-TCR1"
                      & Site != "ED02" & Site != "PCNYN01" & Site != "LCW01"
                      & Site != "LCW02" & Site != "H01-GSC1" & Site != "H01-GSC2"
                      & Site !=  "H01-ROS2" & Site != "H01-SPS2" & Site != "H01-ELY1"
                      & Site != "H01-BGC1" & Site != "ED03" & Site != "H01-CFR1"
                      & Site != "H01-WTP2" & Site != "H01-HNP2" & Site != "LWNT02"
                      & Site != "LWNT01" & Site != "H01-SAB1" & Site != "CHG"
                      & Site != "GKSK" & Site != "UPWE")
new_dat_cut
new_dat
# now remove all seasons of a site that don't exist in detection data set 

# take closer look at detection data set 
uniquesiteseason <- unique(detections[c("Site", "Season")])
uniques <- arrange(uniquesiteseason, Site, Season)
table(uniquesiteseason$Site)

# same with counts 
uniquesiteseason2 <- unique(new_dat_cut[c("Site", "season")])
table(uniquesiteseason2$Site)

# let's remove oct 21 since it doesn't exist in detections data and is probably
# incomplete 

new_dat <- subset(new_dat_cut, season != "OC21")

# merge pasadena and long beach data sets
# skipping this for now but keeping in code just in case
# unique(new_dat$city)
# new_dat$city[new_dat$city == "paca"] <- "mela"
# new_dat$city[new_dat$city == "lbca"] <- "mela"

################################################################################

# merge all deer to "Deer"
new_dat$Species[new_dat$Species == "White-tailed deer"] <- "Deer"
new_dat$Species[new_dat$Species == "Mule deer"] <- "Deer"

# merge all rabbit to "Rabbit"
new_dat$Species[new_dat$Species == "Eastern cottontail rabbit"] <- "Rabbit"
new_dat$Species[new_dat$Species == "Desert cottontail rabbit"] <- "Rabbit"
new_dat$Species[new_dat$Species == "Rabbit (cannot ID)"] <- "Rabbit"

# merge all squirrel to "Squirrel"
# going to keep separate squirrel spp for now 
# new_dat$Species[new_dat$Species == "California Ground Squirrel"] <- "Squirrel"
# new_dat$Species[new_dat$Species == "Douglas squirrel"] <- "Squirrel"
# new_dat$Species[new_dat$Species == "Fox squirrel"] <- "Squirrel"
# new_dat$Species[new_dat$Species == "Western gray squirrel"] <- "Squirrel"

################################################################################

# save data

write_csv(new_dat, "data/count_data_fall20-sum21.csv")


