################################################################################
################# West Coast Env Health GLMM Analysis ##########################
#################      Step 4: Build Model       ###############################
################################################################################



library(readr)
library(MuMIn)
library(car)
library(here)
library(lme4)
library(dplyr)
library(magrittr)
# remotes::install_github("palday/coefplot2",
       #                  subdir = "pkg")
library(coefplot2)

library(corrplot)
library(glmmTMB)

library(aods3)

library(DHARMa)
library(sjPlot)
library(ggplot2)
library(ggeffects)

# read in count and covariate data 
all_dat <- read_csv(here("data", "all_data_counts_covs_10-30-22_J.csv")) 


all_dat$rank_buff
# read in richness/diversity data 
all_dat_veg <- read_csv(here("data", "all_data_vegan_covs_10-27-22.csv")) 


# make things factors
all_dat$City <- as.factor(all_dat$City)
all_dat$season <- as.factor(all_dat$season)
all_dat$Site <- as.factor(all_dat$Site)



# make things factors
all_dat_veg$City <- as.factor(all_dat_veg$City)

# explore data 
colnames(all_dat)



# scale covariates
colnames(all_dat)
all_dat <-  data.frame(all_dat %>% mutate_at(vars("prop_veg",
                                                  "imp_surf", "rank_buff",
                                                  "huden2010", "med_inc", 
                                                  "urb_pca"), scale))

all_dat_veg <-  data.frame(all_dat_veg %>% mutate_at(vars("prop_veg",
                                                  "imp_surf", "rank_buff",
                                                  "huden2010", "med_inc", 
                                                  "urb_pca"), scale))

# get city specific data set 
all_dat_tawa  <- all_dat %>% dplyr::filter(City == "tawa")
all_dat_paca  <- all_dat %>% dplyr::filter(City == "paca")
all_dat_lbca  <- all_dat %>% dplyr::filter(City == "lbca")
all_dat_oaca  <- all_dat %>% dplyr::filter(City == "oaca")


# look at rank distributions
ggplot(all_dat, aes(x = rank_buff)) +
  geom_histogram(fill = "white", colour = "black") +
  facet_grid(City ~ .)

# look at correlations 
cordat <- all_dat %>% dplyr::select(c(urb_pca, rank_buff, med_inc)) 
dev.off()
corr <- cor(cordat, method = "spearman")
corrplot(corr, method = "number")

head(cordat) # all below 0.7 - look ok 

################################################################################
# coyote counts model 

# check for zero inflation for the data set we're interested in 
100*sum(all_dat$Coyote== 0)/nrow(all_dat)

# 39% 0s for the coyotes - we'll try poisson and check for overdispersion 

# universal model with random effects 
mod1 <- glmmTMB(Coyote ~ urb_pca + rank_buff + season + (1|City) +
                (1|Site), 
                data = all_dat, family = "poisson")
summary(mod1)
confint(mod1)

# check for overdispersion
E2 <- resid(mod1, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(mod1))   
sum(E2^2) / (N - p)

# data is overdispersed 

# use negative binomial instead
# can't use lme4 for this distribution in glmms though


mod2 <- glmmTMB(Coyote ~ med_inc + urb_pca + rank_buff + season + (1|City) +
                  (1|Site), offset = log(J),
               data = all_dat, family = "nbinom2")
summary(mod2)

# check for underdispersion
E2 <- resid(mod2, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(mod2))   
sum(E2^2) / (N - p)


# if underdispersed - try zero-inflated poisson
# here we can choose different predictors for count and overdisp - 
# because we might expect different variables to drive presence/absence
# than the ones that drive total # of individuals (count) 

dat2 <- na.omit(all_dat)

zipmod <- glmmTMB(Coyote ~  urb_pca + rank_buff + (1|season) + (1|City) +
                    (1|Site),
                   data = all_dat,
                   ziformula=~1,
                   family="poisson")

summary(zipmod)
E2 <- resid(zipmod, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(zipmod))   
sum(E2^2) / (N - p)

#slightly overdispersed 
# try ZINB 


zinbmod <- glmmTMB(Coyote ~ med_inc + urb_pca + rank_buff + season + (1|City) +
                     (1|Site), 
                  data = all_dat,offset = log(J),
                  ziformula=~1,
                  family="nbinom2")

summary(zinbmod)
E2 <- resid(zinbmod, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(zinbmod))   
sum(E2^2) / (N - p)


# modeldiagnostics 
simOut <- simulateResiduals(fittedModel = mod2, plot = T)
testZeroInflation(mod1)
testOutliers(zipmod)
testDispersion(mod1)
plotResiduals(simOut, form = all_dat$Site)
plot(simOut, quantreg = T)
par(mfrow = c(1,2))
plotResiduals(simOut, all_dat$City)
plotResiduals(simOut, all_dat$Environment2)

# compare AICs 
AIC(mod1, zipmod, mod2, zinbmod)

################################################################################
# city specific 

glmm1 <- glmmTMB(Coyote ~ urb_pca + rank_buff + med_inc + season +  (1|Site), 
                 family = "poisson", offset = log(J),
                 data = all_dat_tawa)



# check for overdispersion
summary(glmm1)
E2 <- resid(glmm1, type = "pearson")
N  <- nrow(all_dat_tawa)
p  <- length(coef(glmm1))   
sum(E2^2) / (N - p)

# nbinom

mod2 <- glmmTMB(Coyote ~ rank_buff + urb_pca + med_inc + season + (1|Site), 
                offset = log(J),
                data = all_dat_tawa, family = "nbinom2")

summary(mod2)

# check for underdispersion
summary(mod2)
E2 <- resid(mod2, type = "pearson")
N  <- nrow(all_dat_tawa)
p  <- length(coef(mod2))   
sum(E2^2) / (N - p)


# if underdispersed - try zero-inflated poisson
# here we can choose different predictors for count and overdisp - 
# because we might expect different variables to drive presence/absence
# than the ones that drive total # of individuals (count) 

zipmodtawa <- glmmTMB(Coyote ~ rank_buff + med_inc + urb_pca + season + (1|Site), 
                #    offset = log(J),
                  data = all_dat_tawa, 
                  ziformula=~1,
                  family="poisson")

summary(zipmodtawa)

# check for underdispersion
summary(zipmodtawa )
E2 <- resid(zipmodtawa , type = "pearson")
N  <- nrow(all_dat_tawa)
p  <- length(coef(zipmodtawa ))   
sum(E2^2) / (N - p)

################################################################################

all_dat_paca <- read_csv(here("data", "all_data_counts_covs_10-30-22_J.csv")) %>% 
  filter(City == "paca")  %>% mutate_at(vars("prop_veg",
                                             "imp_surf", "rank_buff",
                                             "huden2010", "med_inc", 
                                             "urb_pca"), scale)


glmmpaca <- glmmTMB(Coyote ~ rank_buff + urb_pca + med_inc + season + (1|Site), 
                 family = "poisson", offset = log(J),
                 data = all_dat_paca)

summary(glmmpaca)

# check for overdispersion
summary(glmmpaca)
E2 <- resid(glmmpaca, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(glmmpaca))   
sum(E2^2) / (N - p)


# nbinom

mod2 <- glmmTMB(Coyote ~ rank_buff + urb_pca + season + (1|Site), 
                offset = log(J),
                data = all_dat_paca, family = "nbinom2")

summary(mod2)

# check for underdispersion
summary(mod2)
E2 <- resid(mod2, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(mod2))   
sum(E2^2) / (N - p)


# if underdispersed - try zero-inflated poisson
# here we can choose different predictors for count and overdisp - 
# because we might expect different variables to drive presence/absence
# than the ones that drive total # of individuals (count) 

dat2 <- na.omit(all_dat)

zipmodpaca <- glmmTMB(Coyote ~ rank_buff + urb_pca + med_inc + season + (1|Site), 
                  offset = log(J),
                  data = all_dat_paca, 
                  ziformula=~1,
                  family="poisson")

summary(zipmodpaca)

# check for underdispersion
summary(zipmodpaca)
E2 <- resid(zipmodpaca, type = "pearson")
N  <- nrow(all_dat_paca)
p  <- length(coef(zipmodpaca))   
sum(E2^2) / (N - p)


###############################################################################

all_dat_lbca <- read_csv(here("data", "all_data_counts_covs_10-30-22_J.csv")) %>% 
  filter(City == "lbca")  %>% mutate_at(vars("prop_veg",
                                             "imp_surf", "rank_buff",
                                             "huden2010", "med_inc", 
                                             "urb_pca"), scale)


glmm1 <- glmmTMB(Coyote ~ rank_buff + urb_pca + med_inc + season + (1|Site), 
                 family = "poisson", offset = log(J),
                 data = all_dat_lbca)

summary(glmm1)

# check for overdispersion
summary(glmm1)
E2 <- resid(glmm1, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(glmm1))   
sum(E2^2) / (N - p)


# nbinom

mod2 <- glmmTMB(Coyote ~ rank_buff + urb_pca + season + (1|Site), 
                offset = log(J),
                data = all_dat, family = "nbinom2")

summary(mod2)

# check for underdispersion
summary(mod2)
E2 <- resid(mod2, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(mod2))   
sum(E2^2) / (N - p)


# if underdispersed - try zero-inflated poisson
# here we can choose different predictors for count and overdisp - 
# because we might expect different variables to drive presence/absence
# than the ones that drive total # of individuals (count) 

dat2 <- na.omit(all_dat)

zipmodlbca <- glmmTMB(Coyote ~ rank_buff + urb_pca + med_inc + season + (1|Site), 
                  offset = log(J),
                  data = all_dat_lbca, 
                  ziformula=~1,
                  family="poisson")

# check for underdispersion
summary(mod2)
E2 <- resid(mod2, type = "pearson")
N  <- nrow(all_dat_lbca)
p  <- length(coef(mod2))   
sum(E2^2) / (N - p)



##################################################################################
# plot 

mydf <- ggpredict(zipmod, terms = "urb_pca")
summary(zipmod)

ggplot(mydf, aes(x, predicted)) +
  geom_line() +
  labs(
       x ="urbanization_pca", y = "counts") + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)


mydf2 <- ggpredict(zipmod, terms = "rank_buff")

ggplot(mydf2, aes(x, predicted)) +
  geom_line() +
  labs(
       x ="environmental health index", y = "counts") + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)


mydf3 <- ggpredict(zipmod, terms = c("med_inc", "rank_buff"))
confint(zipmod)
?ggpredict
ggplot(mydf3, aes(x, predicted)) +
  geom_line() +
  labs(
    x ="median income", y = "counts") + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)
summary(zipmod)
