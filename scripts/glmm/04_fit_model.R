################################################################################
################# West Coast Env Health GLMM Analysis ##########################
#################      Step 4: Build Model       ###############################
################################################################################


library(readr)
library(here)


# read in data 
LA_dat <- read_csv(here("data", "LA_counts_covs_10-20-22.csv"))
SF_dat <- read_csv(here("data", "SF_counts_covs_10-20-22.csv"))
WA_dat <- read_csv(here("data", "WA_counts_covs_10-20-22.csv"))


mod <- glmer(counts ~ urbanization + income + env_health
           + (1|season) + (1|city) + (1|site), data=coyotecounts)

# other option: make season/city fixed effects and do model selection
mod <- glmer(counts ~ urbanization + income + env_health + 
               season + city + (1|site), data = coyotecounts)