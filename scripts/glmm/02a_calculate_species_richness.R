################################################################################
################# West Coast Env Health GLMM Analysis  #########################
################# Step 2B: Calculate Species Richness, #########################
#################     Diversity, & Abundance           #########################
#################   Yasmine Hentati yhentati@uw.edu    #########################
################################################################################

# packages
library(readr)
library(here)
library(tidyverse)
library(vegan)
library(ggplot2)
library(viridis)
library(cowplot)
library(dplyr)

# read in count data
counts <- read_csv("data/count_data_fall20-sum21.csv")

# add together all occurrences of each spp per site (for grouped spp) 
counts_sum <- counts %>% group_by(Species, Site, season) %>% summarise(ySum = sum(count))
counts_sum


# pivot into wide format 
wide_counts <- counts_sum %>% pivot_wider(names_from = "Species", values_from = "ySum")
wide_counts

# replace all NAs with 0s
wide_counts <- wide_counts %>% mutate_all(~replace(., is.na(.), 0))
wide_counts

# now add in sites that have 0 observations of anything and SR of 0 
# going to use detection data to do this 

# read in detection history data set
detections <- read_csv("data/raw_data_from_uwin/initial_data_yasmine.csv")

# remove all rows where J is 0 
detections <- detections[detections$J != 0, ]
(unique(detections$Site))

# summarize Ys for each site 

det_1 <- detections %>% group_by(Site, Season) %>% summarize(ySum = sum(Y))

# now only keep rows where Y is 0 - this means nothing was detected for 
# that site/season combination

det_0 <- det_1 %>% subset(ySum == 0)
det_0

colnames(det_0) <- c("Site", "season", "ySum")

# make det_0 look like wide_counts by adding columns 
det_0 <- det_0 %>% dplyr::select(-ySum) %>% add_column("Black bear" = 0, "Bobcat" = 0, "Brush Rabbit" = 0,
                    "California Ground Squirrel" = 0, "Coyote" = 0, "Deer" = 0,
                    "Domestic cat" = 0, "Domestic dog" = 0, "Douglas squirrel" = 0, 
                    "Eastern gray squirrel" = 0, "Elk" = 0, "Fox squirrel" = 0, 
                    "Gray fox" = 0, "Mountain lion" = 0, "Rabbit" = 0,"Raccoon" = 0, 
                    "Red fox" = 0, "Stoat" = 0, "Striped Skunk" = 0, "Virginia opossum" = 0,
                    "Weasel (cannot ID)" = 0, "Western Chipmunks" = 0, "Western gray squirrel" = 0)


# now we need to add a new row to the data
# for each site that already exists but had no detections in a season

# bind both dfs 
counts_all <- rbind(wide_counts, det_0)
counts_all

# keep first unique site/season combo
counts_all <- counts_all %>% distinct(Site, season, .keep_all = TRUE)

# merge sites that are the same site with diff code
counts_all$Site[counts_all$Site == "H01-FIR1"] <- "H01-FIR2" # keep

# pivot back to long format 
counts_long <- pivot_longer(counts_all, cols = 3:25, names_to = "Species")

################################################################################

# add metadata - some new sites from detections so need to bind those 
counts_meta <- counts %>% dplyr::select(-c("Species", "season", "locationID", "count", "utmEast",
                                    "utmNorth", "utmZone"))
det_meta <- detections %>% dplyr::select(-c("Species", "Season", "Crs", "Y", "J"))

# give same column names 
colnames(counts_meta) <- c("City", "Site", "Long", "Lat")

# bind together 
metadata <- rbind(counts_meta, det_meta)
metadata <- metadata %>% distinct(Site, .keep_all = TRUE)

# change "mela" to "paca" or "lbca" 
metadata[76,1] = "paca"
metadata$City[metadata$City == "mela"] <- "lbca"

# merge metadata with long count data
counts_long <- counts_long %>% left_join(metadata, by = "Site")


## stuck - need to figure out how to remove 0 rows while keeping those sites 
# if 0 for all

# get richness for each site/season
sr<- counts_sum %>%
    group_by(Site, season) %>%
  summarise(species.richness=n()) %>%
  arrange(-species.richness)


# average richness across seasons for each site 
# skipping this for now because our sites have different amounts of seasons 
# plr_avg <- plr %>% group_by(Site) %>% summarise(avg_richness = mean(species.richness))
# hist(plr_avg$avg_richness)

# save as csv
write_csv(sr, "data/spp_rich_fall20-sum21.csv")


################################################################################

## calculating diversity 


div <- counts_long %>%
  group_by(Site, season) %>%
  filter(value>0) %>%
  summarise(N=sum(value),
            shannon.di=-sum((value/sum(value))*log(value/sum(value))),
            simpson.di=1-sum((value/sum(value))^2),
            inv.simpson.di=1/sum((value/sum(value))^2)) %>%
  arrange(-shannon.di)

# merge with metadata
div <- div %>% left_join(metadata, by = "Site")

# save as csv
# write_csv(div, "data/diversity-fall20-sum21.csv")

div <- read_csv("data/diversity-fall20-sum21.csv")


## get species abundance data sets 
# going to use counts_all from before

colnames(counts_all)

# pivot back to long format 
counts_long <- pivot_longer(counts_all, cols = 3:25, names_to = "Species")

# add metadata - some new sites from detections so need to bind those 
counts_meta <- counts %>% dplyr::select(-c("Species", "season", "locationID", "count", "utmEast",
                                    "utmNorth", "utmZone"))
det_meta <- detections %>% dplyr::select(-c("Species", "Season", "Crs", "Y", "J"))

# give same column names 
colnames(counts_meta) <- c("City", "Site", "Long", "Lat")

# bind together 
metadata <- rbind(counts_meta, det_meta)
metadata <- metadata %>% distinct(Site, .keep_all = TRUE)

# change "mela" to "paca" or "lbca" 
metadata[76,1] = "paca"
metadata$City[metadata$City == "mela"] <- "lbca"

# merge metadata with long count data
counts_long <- counts_long %>% left_join(metadata, by = "Site")

counts_long$City
# pivot back to wide 
counts_wide_all <- counts_long %>% pivot_wider(names_from = "Species", 
                                               values_from = "value")
colnames(counts_wide_all)

# create csv
# write_csv(counts_wide_all, here("data", "counts_cleaned.csv"))
# write_csv(counts_long, here("data", "counts_cleaned_long.csv"))



################################################################################
## same as above but using vegan

# with this method, we will get rid of seasons altogether 
counts_year <- counts_long  %>% dplyr::select(-season) %>% 
  group_by(Site, Species) %>% summarize(value = sum(value))
glimpse(counts_year)

# pivot to wide 
counts_wide_yr <- counts_year %>% pivot_wider(names_from = "Species", 
                                               values_from = "value")
glimpse(counts_wide_yr)

# calculate richness
richness <- specnumber(counts_wide_yr) 
counts_wide_yr$richness <- richness
counts_wide_yr

# shannon's div
div_vegan <- counts_year %>%
  group_by(Site) %>% 
  summarise(N=sum(value),
            shannon.di=diversity(value, index = "shannon", MARGIN = 2),
            simpson.di=diversity(value, index = "simpson", MARGIN = 2),
            inv.simpson.di=diversity(value, index = "invsimpson", MARGIN = 2)) %>%
  arrange(-shannon.di)

# bind together
rich_vegan <- counts_wide_yr %>% dplyr::select(Site, richness)

vegandf <- left_join(div_vegan, rich_vegan, by = "Site")

# pielou's evenness
evenness <- vegandf$shannon.di/log(vegandf$richness)
vegandf$evenness <- evenness

head(vegandf)

# bind back with metadata 
vegandf <- vegandf %>% left_join(metadata, by = "Site")
vegandf

write_csv(vegandf, "vegan_sites.csv")

################################################################################

