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

# read in data

counts <- read_csv("data/count_data_fall20-sum21.csv")

# add together all occurrences of each spp per site 
counts_sum <- counts %>% group_by(Species, locationID, season) %>% summarise(ySum = sum(count))

# pivot into wide format 
# wide_counts <- counts_sum %>% pivot_wider(names_from = "Species", values_from = "ySum")

# replace all NAs with 0s
# wide_counts %>% mutate_all(~replace(., is.na(.), 0))

# get richness for each site/season
plr<-counts_sum %>%
  group_by(locationID, season) %>%
  summarise(species.richness=n()) %>%
  arrange(-species.richness)

# average richness across seasons for each site 
plr_avg <- plr %>% group_by(locationID) %>% summarise(avg_richness = mean(species.richness))

## merge with other site data 

# change counts data to only have 1 row per site and remove variables we don't care about

counts_new <- counts %>% select(-c(Species, count, season)) %>% distinct(locationID, .keep_all=TRUE)

# merge site metadata with richness
spp_rich_avg <- plr_avg %>% left_join(counts_new)
spp_rich_avg

# now add in sites that have 0 observations of anything and SR of 0 
# going to use detection data to do this 

# save as csv
write_csv(spp_rich_avg, "data/avg_spp_rich_fall20-sum21.csv")



