################################################################################
################# West Coast Env Health GLMM Analysis ##########################
################# Step 2B: Calculate Species Richness ##########################
#################   Yasmine Hentati yhentati@uw.edu   ##########################
################################################################################

# packages
library(readr)
library(here)
library(tidyverse)
library(vegan)

# read in detection data
counts <- read_csv("data/count_data_fall20-sum21.csv")

# add together all occurrences of each spp per site 
counts_sum <- counts %>% group_by(Species, Site, season) %>% summarise(ySum = sum(count))

# pivot into wide format 
# wide_counts <- counts_sum %>% pivot_wider(names_from = "Species", values_from = "ySum")

# replace all NAs with 0s
# wide_counts %>% mutate_all(~replace(., is.na(.), 0))



# now add in sites that have 0 observations of anything and SR of 0 
# going to use detection data to do this 

# read in detection history data set
detections <- read_csv("data/raw_data_from_uwin/initial_data_yasmine.csv")

# summarize Ys for each site 

det_1 <- detections %>% group_by(Site, Season) %>% summarize(ySum = sum(Y))

# now only keep rows where Y is 0 - this means nothing was detected for 
# that site/season combination

det_0 <- det_1 %>% subset(ySum == 0)

# now we have all site/season combos where nothing was detected 
intersect(det_0$Site, counts_sum$site)

# now we need to add a new row to the richness data
# for each site that already exists but had no detections in a season












# get richness for each site/season
plr<-counts_sum %>%
  group_by(locationID, season) %>%
  summarise(species.richness=n()) %>%
  arrange(-species.richness)


# average richness across seasons for each site 
# skipping this for now because our sites have different amounts of seasons 
# plr_avg <- plr %>% group_by(Site) %>% summarise(avg_richness = mean(species.richness))

## merge with other site data 

# change counts data to only have 1 row per site and remove variables we don't care about

counts_new <- counts %>% select(-c(Species, count, season)) %>% distinct(locationID, .keep_all=TRUE)

# merge site metadata with richness
spp_rich_dat<- plr %>% left_join(counts_new)
spp_rich_dat



# save as csv
write_csv(spp_rich_dat, "data/avg_spp_rich_fall20-sum21.csv")



