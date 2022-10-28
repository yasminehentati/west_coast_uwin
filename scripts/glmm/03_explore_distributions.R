################################################################################
################# West Coast Env Health GLMM Analysis ##########################
#################      Step 3: Explore Data      ###############################
################################################################################


# packages
library(car)
library(MASS)
library(lme4)
library(readr)
library(here)
library(tidyverse)

# read in count and covariate data 
all_dat <- read_csv(here("data", "all_data_counts_covs_10-27-22.csv"))

# read in richness/diversity data 
all_dat_veg <- read_csv(here("data", "all_data_vegan_covs_10-27-22.csv"))

colnames(all_dat_veg)
hist(all_dat_veg$shannon.di)
hist(all_dat_veg$richness)


boxplot()

data2$shannondi <- data2$shannon.di+0.00001

all_dat$coyote <- all_dat$Coyote+0.00001
# look at distribution of data 

qqp(data2$richness, "norm")
qqp(data2$richness, "lnorm")


# qqp requires estimates of the parameters of the negative binomial, Poisson
# and gamma distributions. You can generate estimates using the fitdistr
# function. Save the output and extract the estimates of each parameter as I
# have shown below.

nbinom <- fitdistr(all_dat$Raccoon, "Negative Binomial")
qqp(all_dat$Raccoon, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

poisson <- fitdistr(all_dat$Raccoon, "Poisson")
qqp(all_dat$Raccoon, "pois",  lambda = poisson$estimate[[1]])

gamma <- fitdistr(data3$coyote, "gamma")
qqp(data3$coyote, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

# read data - div
data <- read_csv("data/diversity-fall20-sum21.csv")
glimpse(data)
data$shannon.di.t <- data$shannon.di + 1 
# look at distribution of data 

qqp(data$shannon.di.t, "norm")
qqp(data$shannon.di.t, "lnorm")


# qqp requires estimates of the parameters of the negative binomial, Poisson
# and gamma distributions. You can generate estimates using the fitdistr
# function. Save the output and extract the estimates of each parameter as I
# have shown below.
nbinom <- fitdistr(data$shannon.di.t, "Negative Binomial")
qqp(data$shannon.di.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

poisson <- fitdistr(data$shannon.di.t, "Poisson")
qqp(data$shannon.di.t, "pois", lambda = poisson$estimate[[1]])

gamma <- fitdistr(data$shannon.di.t, "gamma")
qqp(data$shannon.di.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

?qqp
