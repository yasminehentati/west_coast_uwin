################################################################################
################# West Coast Env Health GLMM Analysis ##########################
#################      Step 1: Clean Count Data       ##########################
################################################################################


##################### 
# install.packages("_")
library(data.table)
library(sf)
library(igraph)
library(dplyr)
library(here)
library(readr)

# read in data
all_data <- read_csv("yasmine_west_coast_data.csv")
length(unique(all_data$locationAbbr))
?read_csv
class(all_data)

# check data
head(all_data)
nrow(all_data)
unique(all_data$commonName)

# read in detection history data set
detections <- read_csv("initial_data_yasmine.csv")

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
  select(locationID, locationAbbr, utmEast, utmNorth, utmZone) %>% distinct()
head(all_data)

# join counts to site data 
joined_data <- left_join(counts, all_data, by = "locationID")
nrow(joined_data)
length(unique(joined_data$locationAbbr))

View(joined_data)

# now we have counts for each site/season combo 

# let's change UTM to lat long
library(terra)

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
# now use mason's functions to collapse sites by location 
source("mason_site_collapse/qaqc_sites.R")
source("mason_site_collapse/long_to_zone.R")
source("mason_site_collapse/fix_site_names.R")

new_dat_1 <- qaqc_sites(x = new_dat, cities="city", sites = "Site",
                        my_coords=c("long","lat"), my_crs=4326)

fix_site_names(new_dat_1)

# merge  H01-SPS1 and H01-SPS2; H01-HNP1 and H01-HNP2; H01-CFR2 and H01-CFR1

new_dat$Site[new_dat$Site == "H01-SPS2"] <- "H01-SPS1"
new_dat$Site[new_dat$Site == "H01-HNP2"] <- "H01-HNP1"
new_dat$Site[new_dat$Site == "H01-CFR2"] <- "H01-CFR2"

# merge pasadena and long beach data sets
unique(new_dat$city)
new_dat$city[new_dat$city == "paca"] <- "mela"
new_dat$city[new_dat$city == "lbca"] <- "mela"



# merge all deer to "Deer"
new_dat$Species[new_dat$Species == "White-tailed deer"] <- "Deer"
new_dat$Species[new_dat$Species == "Mule deer"] <- "Deer"


# merge all rabbit to "Rabbit"
new_dat$Species[new_dat$Species == "Eastern cottontail rabbit"] <- "Rabbit"
new_dat$Species[new_dat$Species == "Desert cottontail rabbit"] <- "Rabbit"
new_dat$Species[new_dat$Species == "Rabbit (cannot ID)"] <- "Rabbit"

# merge all squirrel to "Squirrel"
new_dat$Species[new_dat$Species == "California Ground Squirrel"] <- "Squirrel"
new_dat$Species[new_dat$Species == "Douglas squirrel"] <- "Squirrel"
new_dat$Species[new_dat$Species == "Fox squirrel"] <- "Squirrel"
new_dat$Species[new_dat$Species == "Western gray squirrel"] <- "Squirrel"

unique(new_dat$Species)


# separate data set for each species 


