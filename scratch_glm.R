

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
library(pscl)


#########################################################################
# tacoma - coyotes

# read in count and covariate data 
all_dat <- read_csv(here("data", "all_data_counts_covs_10-30-22_J.csv")) %>% 
  filter(City == "tawa")  %>% mutate_at(vars("prop_veg",
                                             "imp_surf", "rank_buff",
                                             "huden2010", "med_inc", 
                                             "urb_pca"), scale)

glm1 <- glm(Coyote  ~ urb_pca + season, 
            family = "poisson", offset = log(J), data = all_dat_tawa)

summary(glm1)
glmm1 <- glmmTMB(Coyote ~ rank_buff + urb_pca + med_inc + season + (1|Site), 
    family = "poisson", offset = log(J),
    data = all_dat)

summary(glmm1)

zeroinfl(Coyote ~ rank_buff + urb_pca + med_inc + season | ## Predictor for the Poisson process
           rank_buff + urb_pca + med_inc + season, ## Predictor for the Bernoulli process;
         dist = 'poisson', offset = log(J),
         data = all_dat)

glm1 <- glmer(Coyote ~ rank_buff + prop_veg + med_inc + (1|season), 
            family = "poisson", offset = log(J),
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


#########################################################################

# tacoma - richness

# read in richness/diversity data 
all_dat_veg <- read_csv(here("data", "all_data_vegan_covs_10-27-22.csv"))  %>% 
  filter(City == "tawa")  %>% mutate_at(vars("prop_veg",
                                             "imp_surf", "rank_buff",
                                             "huden2010", "med_inc", 
                                             "urb_pca"), scale)


all_dat_veg
head(all_dat)

glm1 <- glm(shannon.di ~ rank_buff + huden2010 + med_inc,
            family = "gaussian", 
            data = all_dat_veg)

summary(glm1)


## plot 

newdat_rank <- data.frame(rank_buff = 2:10, veg_cover = mean(all_dat$prop_veg),
                          median_income = mean(all_dat$med_inc), samplingOccasions = 1)

pred_rank <- predict(glm1, newdata = newdat_rank, se.fit = T)
newdat_rank$fit.link <- pred_rank$fit
newdat_rank$se.link <- pred_rank$se.fit

newdat_rank$fit <- exp(pred_rank$fit)
newdat_rank$lcl <- exp(newdat_rank$fit.link - 1.96*newdat_rank$se.link)
newdat_rank$ucl <- exp(newdat_rank$fit.link + 1.96*newdat_rank$se.link)





# sjplot

?plot_model
modlbca <- plot_model(zinbmodlbca, title = "Long Beach, CA",
                      rm.terms = "season [JU21]",
                      show.values = TRUE, show.p = TRUE) +   
  scale_color_sjplot("circus") 

modlbca + scale_y_log10(limits = c(0.1, 10))

modpaca <- plot_model(glmmpaca, type = "pred", title = "Pasadena, CA",
                      terms = "med_inc", axis.title = "",
                      show.values = TRUE, show.p = TRUE) +   
  scale_color_sjplot("circus")
modpaca + scale_y_log10(limits = c(0.1, 10))



pr1 <- ggpredict(modtawa)

# plot multiple models
mods <- plot_models(list(zipmodtawa, glmmpaca, zipmodlbca), show.values = T, grid = TRUE, 
                    rm.terms = "season [JA21, JU21]")
?plot_models
mods + scale_y_log10(limits = c(0.3, 2.3)) +  theme_sjplot2() + 
  scale_color_sjplot("simply")


tab_model(mod1,
          show.reflvl = T, 
          show.intercept = F)


ggpredict(mod1, c("rank_buff", "urb_pca", "City")) %>% plot()
?plot_model
install.packages("effects")
library(effects)
confint(glmm1)
modtawa <- plot_model(mod1,  type = "pred", title = "Tacoma, WA",
                      terms = c("urb_pca"), axis.title = "",
                      show.values = TRUE, show.p = TRUE,) +   
  scale_color_sjplot("circus") 


plot1 <- ggpredict(mod1, terms = c("med_inc", "City [tawa]"), type = "re")
confint(zipmod)
plot(plot1)

ggpredict(zipmod, terms = c("rank_buff", "med_inc", "City")) %>% plot()

modtawa <- plot_model(zipmod,  type = "eff", title = "Tacoma, WA",
                      terms = c("med_inc", "City"), axis.title = "") +   
  scale_color_sjplot("circus") 

?ggpredict
modlbca <- plot_model(zipmodlbca, title = "Long Beach, CA",
                      rm.terms = "season [JU21]", axis.title = "",
                      show.values = TRUE, show.p = TRUE) +   
  scale_color_sjplot("circus") 
modlbca + scale_y_log10(limits = c(0.1, 10))

modpaca <- plot_model(glmmpaca, type = "pred", title = "Pasadena, CA",
                      terms = "med_inc", axis.title = "",
                      show.values = TRUE, show.p = TRUE) +   
  scale_color_sjplot("circus")
modpaca + scale_y_log10(limits = c(0.1, 10))



pr1 <- ggpredict(modtawa)

# plot multiple models
mods <- plot_models(list(zipmodtawa, glmmpaca, zipmodlbca), show.values = T, grid = TRUE, 
                    rm.terms = "season [JA21, JU21]"))
?plot_models
mods +  +  theme_sjplot2() + 
  scale_color_sjplot("simply")


tab_model(mod1,
          show.reflvl = T, 
          show.intercept = F)

