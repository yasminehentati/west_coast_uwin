


# packages
library(car)
library(MASS)
library(lme4)
library(readr)
library(here)
library(tidyverse)

# read data
data <- read_csv("data/spp_rich_fall20-sum21.csv")
glimpse(data)
data$species.richness

# look at distribution of data 

qqp(data$species.richness, "norm")
qqp(data$species.richness, "lnorm")


# qqp requires estimates of the parameters of the negative binomial, Poisson
# and gamma distributions. You can generate estimates using the fitdistr
# function. Save the output and extract the estimates of each parameter as I
# have shown below.
nbinom <- fitdistr(data$species.richness, "Negative Binomial")
qqp(data$species.richness, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

poisson <- fitdistr(data$species.richness, "Poisson")
qqp(data$species.richness, "pois",  lambda = poisson$estimate[[1]])

gamma <- fitdistr(data$species.richness, "gamma")
qqp(data$species.richness, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

?qqp

mod <-lmer(yield ~irrigation*density*fertilizer +(1|block/irrigation/density), data=crops)

# read data - avg
data <- read_csv("data/avg_spp_rich_fall20-sum21.csv")
glimpse(data)

# look at distribution of data 

qqp(data$avg_richness, "norm")
qqp(data$avg_richness, "lnorm")


# qqp requires estimates of the parameters of the negative binomial, Poisson
# and gamma distributions. You can generate estimates using the fitdistr
# function. Save the output and extract the estimates of each parameter as I
# have shown below.
nbinom <- fitdistr(data$avg_richness, "Negative Binomial")
qqp(data$avg_richness, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

poisson <- fitdistr(data$avg_richness, "Poisson")
qqp(data$avg_richness, "pois", lambda = poisson$estimate[[1]])

gamma <- fitdistr(data$avg_richness, "gamma")
qqp(data$avg_richness, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

?qqp

mod <-lmer(yield ~irrigation*density*fertilizer +(1|block/irrigation/density), data=crops)