


# packages
library(car)
library(MASS)
library(lme4)
library(readr)
library(here)
library(tidyverse)

# read data
data <- read_csv("data/avg_spp_rich_fall20-sum21.csv")
glimpse(data)

# look at distribution of data 
data$avg_richness <- recog$Aggression + 1
qqp(recog$Aggression.t, "norm")

mod <-lmer(yield ~irrigation*density*fertilizer +(1|block/irrigation/density), data=crops)