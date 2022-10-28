

library(rmarkdown)
library(MuMIn)
library(corrplot)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)
library(data.table)
library(extrafont)
library(ggeffects)
library(jtools)
library(RColorBrewer)
library(cowplot)
library(grid)
library(png)

# tacoma 

# read in count and covariate data 
all_dat <- read_csv(here("data", "all_data_counts_covs_10-27-22.csv")) %>% 
  filter(City == "oaca") %>% group_by(Site) 

?summarize

all_dat
head(all_dat)

glm1 <- glm(Coyote ~ rank_buff + urb_pca + med_inc, 
    family = "poisson",
    data = all_dat)

summary(glm1)
?glm

## plot 

newdat_rank <- data.frame(rank_buff = 2:10, veg_cover = mean(all_dat$prop_veg),
                          median_income = mean(all_dat$med_inc), samplingOccasions = 1)

pred_rank <- predict(glm1, newdata = newdat_rank, se.fit = T)
newdat_rank$fit.link <- pred_rank$fit
newdat_rank$se.link <- pred_rank$se.fit

newdat_rank$fit <- exp(pred_rank$fit)
newdat_rank$lcl <- exp(newdat_rank$fit.link - 1.96*newdat_rank$se.link)
newdat_rank$ucl <- exp(newdat_rank$fit.link + 1.96*newdat_rank$se.link)
