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

dets_all <- read_csv("data/raw_data_from_uwin/initial_data_yasmine.csv")
glimpse(dets_all)

# remove all sites with 0 days running from detection data
dets <- dets_all[dets_all$J != 0, ]
(unique(dets$Site))

# add together all occurrences of each spp per site 
dets_sum <- dets %>% group_by(Species, Site, Season) %>% summarise(ySum = sum(Y))

dets_sum

# pivot into wide format 
# wide_counts <- counts_sum %>% pivot_wider(names_from = "Species", values_from = "ySum")

# replace all NAs with 0s
# wide_counts %>% mutate_all(~replace(., is.na(.), 0))

# get richness for each site/season
plr<- dets_sum %>%
  group_by(Site, Season) %>%
  summarise(species.richness=n()) %>%
  arrange(-species.richness)

dets_sum$ySum
# average richness across seasons for each site 
# skipping this for now because our sites have different amounts of seasons 
# plr_avg <- plr %>% group_by(Site) %>% summarise(avg_richness = mean(species.richness))

## merge with other site data 

# change counts data to only have 1 row per site and remove variables we don't care about

dets_new <- dets_all %>% select(-c(Species, Y, Season)) %>% distinct(Site, .keep_all=TRUE)

# merge site metadata with richness
spp_rich_avg <- plr %>% left_join(dets_all)
spp_rich_avg

# now add in sites that have 0 observations of anything and SR of 0 
# going to use detection data to do this 

# save as csv
write_csv(spp_rich_avg, "data/avg_spp_rich_fall20-sum21.csv")



