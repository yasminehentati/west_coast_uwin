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

# read data
data <- read_csv("data/spp_rich_fall20-sum21.csv")
data2 <- read_csv("vegan_sites.csv")
data3 <- read_csv("data/counts_cleaned.csv")
data3$Coyote
glimpse(data2)

hist(data2$shannon.di)
hist(data$species.richness)
data$species.richness

data2$shannondi <- data2$shannon.di+0.00001

data3$coyote <- data3$Coyote+0.00001
# look at distribution of data 

qqp(data2$richness, "norm")
qqp(data2$richness, "lnorm")


# qqp requires estimates of the parameters of the negative binomial, Poisson
# and gamma distributions. You can generate estimates using the fitdistr
# function. Save the output and extract the estimates of each parameter as I
# have shown below.

nbinom <- fitdistr(data3$coyote, "Negative Binomial")
qqp(data3$coyote, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

poisson <- fitdistr(data3$coyote, "Poisson")
qqp(data3$coyote, "pois",  lambda = poisson$estimate[[1]])

gamma <- fitdistr(data3$coyote, "gamma")
qqp(data3$coyote, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

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
