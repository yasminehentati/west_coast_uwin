################################################################################
################# West Coast Env Health GLMM Analysis ##########################
#################      Step 4: Build Model       ###############################
################################################################################



library(readr)
library(here)
library(lme4)
library(dplyr)
library(magrittr)
# remotes::install_github("palday/coefplot2",
       #                  subdir = "pkg")
library(coefplot2)
library(glmmTMB)
library(aods3)
library(DHARMa)
install.packages("sjPlot")


# read in count and covariate data 
all_dat <- read_csv(here("data", "all_data_counts_covs_10-27-22.csv")) 

View(all_dat)
# read in richness/diversity data 
all_dat_veg <- read_csv(here("data", "all_data_vegan_covs_10-27-22.csv")) 


# make things factors
all_dat$City <- as.factor(all_dat$City)
all_dat$season <- as.factor(all_dat$season)
all_dat$Site <- as.factor(all_dat$Site)

# explore data 
colnames(all_dat)



# scale covariates
colnames(all_dat)
all_dat <-  data.frame(all_dat %>% mutate_at(vars("prop_veg",
                                                  "imp_surf", "rank_buff",
                                                  "huden2010", "med_inc", 
                                                  "urb_pca"), scale))
View(all_dat)
all_dat_veg <-  data.frame(all_dat_veg %>% mutate_at(vars("prop_veg",
                                                  "imp_surf", "rank_buff",
                                                  "huden2010", "med_inc", 
                                                  "urb_pca"), scale))


ggplot(all_dat, aes(x = rank_buff)) +
  geom_histogram(fill = "white", colour = "black") +
  facet_grid(City ~ .)

################################################################################
# coyote counts model 

# check for zero inflation for the data set we're interested in 
100*sum(all_dat$Coyote== 0)/nrow(all_dat)

# 39% 0s for the coyotes - we'll try poisson and check for overdispersion 

# universal model with random effects 
mod1 <- glmmTMB(Coyote ~ urb_pca + rank_buff + med_inc + (1|season),
                data = all_dat, family = "poisson")

summary(mod1)

# check for overdispersion
E2 <- resid(mod1, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(mod1))   
sum(E2^2) / (N - p)

# data is overdispersed 

# use negative binomial instead
# can't use lme4 for this distribution in glmms though


mod2 <- glmmTMB(Coyote ~ urb_pca + rank_buff + med_inc + City + season +
                  (1|City:Site), 
               data = all_dat, family = "nbinom2")


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

zipmod <- glmmTMB(Coyote ~ urb_pca + rank_buff + med_inc + (1|Site) + (1|season) +
                    (1|City), 
                   data = all_dat,
                   ziformula=~1,
                   family="poisson")



summary(zipmod)
E2 <- resid(zipmod, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(zipmod))   
sum(E2^2) / (N - p)

# only slightly overdispersed 

# modeldiagnostics 
simOut <- simulateResiduals(fittedModel = zipmod, plot = T)
testZeroInflation(mod1)
testOutliers(zipmod)
testDispersion(mod1)
plotResiduals(simOut, form = all_dat$Site)
plot(simOut, quantreg = T)
par(mfrow = c(1,2))
plotResiduals(simOut, all_dat$City)
plotResiduals(simOut, all_dat$Environment2)

# compare AICs 
AIC(mod1, zipmod, mod2)

# compare our 2 zero-inflated models 
# lrtest(zipmod, M4)



################################################################################
# raccoon counts model 

# check for zero inflation for the data set we're interested in 
100*sum(all_dat$Raccoon== 0)/nrow(all_dat)

# 56% 0s for the Raccoon - we'll try poisson and check for overdispersion 

# universal model with random effects 
mod1 <- glmer(Raccoon ~ urb_pca + rank_buff + med_inc  + (1|City) + 
                (1|season) + (1|Site), 
              data=all_dat, family = "poisson", 
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

# check for overdispersion
summary(mod1)
E2 <- resid(mod1, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(mod1))   
sum(E2^2) / (N - p)

# data is overdispersed 

# other option: make season/city fixed effects and do model selection
# mod <- glmer(Coyote ~ prop_veg + rank_buff + 
#   season + City + (1|Site), data = all_dat, family = poisson)
# 


# use negative binomial instead
# can't use lme4 for this distribution in glmms though


mod2 <- glmmTMB(Raccoon ~ urb_pca + rank_buff + med_inc  + (1|City) + 
                  (1|season) + (1|Site), 
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

# model didn't converge so lessened random effects 
zipmod <- glmmTMB(Raccoon ~ urb_pca + rank_buff + med_inc  + City + 
                    season + (1|Site), 
                  data = dat2,
                  ziformula=~1,
                  family="poisson")



summary(zipmod)
E2 <- resid(zipmod, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(zipmod))   
sum(E2^2) / (N - p)

# only  slightly overdispersed 

# compare AICs 
AIC(mod1, zipmod, mod2)

# if there's overdispersion of ZIP - use a ZINB 

M4 <- zeroinfl(Counts ~ Canopy.std |
                 Canopy.std,
               dist = 'negbin',
               data = OBFL)

# check dispersion
E2 <- resid(M4, type = "pearson")
N  <- nrow(OBFL)
p  <- length(coef(M4)) + 1 # '+1' is due to theta
sum(E2^2) / (N - p)


# compare our 2 zero-inflated models 
lrtest(M3, M4)
summary(M4)





################################################################################
# species richness

all_dat_veg_ca <- all_dat_veg %>% filter(City != "tawa")
all_dat_veg_ca
# universal model with random effects )
mod1 <- glmmTMB(richness ~ urb_pca + rank_buff + med_inc + (1|City),
data=all_dat_veg_ca, family = "poisson")

summary(mod1)


# modeldiagnostics 
simOut <- simulateResiduals(fittedModel = zipmod, plot = T)
testZeroInflation(mod1)
testOutliers(zipmod)
testDispersion(mod1)
plotResiduals(simOut, form = all_dat$Site)
plot(simOut, quantreg = T)
par(mfrow = c(1,2))
plotResiduals(simOut, all_dat$City)
plotResiduals(simOut, all_dat$Environment2)

# check for overdispersion
summary(mod1)
E2 <- resid(mod1, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(mod1))   
sum(E2^2) / (N - p)

# data is underdispersed 


# if underdispersed - try zero-inflated poisson
# here we can choose different predictors for count and overdisp - 
# because we might expect different variables to drive presence/absence
# than the ones that drive total # of individuals (count) 

dat2 <- na.omit(all_dat)

# model didn't converge so lessened random effects 
zipmod <- glmmTMB(Raccoon ~ urb_pca + rank_buff + med_inc  + City + 
                    season + (1|Site), 
                  data = dat2,
                  ziformula=~1,
                  family="poisson")



summary(zipmod)
E2 <- resid(zipmod, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(zipmod))   
sum(E2^2) / (N - p)

# only  slightly overdispersed 

# compare AICs 
AIC(mod1, zipmod, mod2)

# if there's overdispersion of ZIP - use a ZINB 

M4 <- zeroinfl(Counts ~ Canopy.std |
                 Canopy.std,
               dist = 'negbin',
               data = OBFL)

# check dispersion
E2 <- resid(M4, type = "pearson")
N  <- nrow(OBFL)
p  <- length(coef(M4)) + 1 # '+1' is due to theta
sum(E2^2) / (N - p)


# compare our 2 zero-inflated models 
lrtest(M3, M4)
summary(M4)





################################################################################
# shannon diversity


# universal model with random effects 
mod1 <- glmmTMB(shannon.di ~ urb_pca + rank_buff + med_inc + (1|City) + (1|Site),
              data=all_dat_veg, family = "gaussian")

summary(mod1)

simulateResiduals(mod1,plot=T)



################################################################################
# deer counts model 

# check for zero inflation for the data set we're interested in 
100*sum(all_dat$Deer== 0)/nrow(all_dat)

# 56% 0s for the Deer - we'll try poisson and check for overdispersion 

# universal model with random effects 
mod1 <- glmmTMB(Deer ~ urb_pca + rank_buff + med_inc + (1|City) + (1|season)+
                  (1|Site), 
                data = all_dat, family = "poisson")

# check for overdispersion
summary(mod1)
E2 <- resid(mod1, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(mod1))   
sum(E2^2) / (N - p)

# data is overdispersed 

# other option: make season/city fixed effects and do model selection
# mod <- glmer(Coyote ~ prop_veg + rank_buff + 
#   season + City + (1|Site), data = all_dat, family = poisson)
# 


# use negative binomial instead


mod2 <- glmmTMB(Deer ~ urb_pca + rank_buff + med_inc  + (1|City) + 
                  (1|season) + (1|Site), 
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

# model didn't converge so lessened random effects 
zipmod <- glmmTMB(Deer ~ urb_pca + rank_buff + med_inc  + City + 
                    season + (1|Site), 
                  data = dat2,
                  ziformula=~1,
                  family="poisson")



summary(zipmod)
E2 <- resid(zipmod, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(zipmod))   
sum(E2^2) / (N - p)

# only  slightly overdispersed 

# compare AICs 
AIC(mod1, zipmod, mod2)

# if there's overdispersion of ZIP - use a ZINB 

M4 <- zeroinfl(Counts ~ Canopy.std |
                 Canopy.std,
               dist = 'negbin',
               data = OBFL)

# check dispersion
E2 <- resid(M4, type = "pearson")
N  <- nrow(OBFL)
p  <- length(coef(M4)) + 1 # '+1' is due to theta
sum(E2^2) / (N - p)


# compare our 2 zero-inflated models 
lrtest(M3, M4)
summary(M4)



################################################################################
# deer counts model 

# check for zero inflation for the data set we're interested in 
100*sum(all_dat$Deer== 0)/nrow(all_dat)

# 56% 0s for the Deer - we'll try poisson and check for overdispersion 

# universal model with random effects 
mod1 <- glmer(Deer ~ urb_pca + rank_buff + med_inc  + (1|City) + 
                (1|season) + (1|Site), 
              data=all_dat, family = "poisson", 
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

# check for overdispersion
summary(mod1)
E2 <- resid(mod1, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(mod1))   
sum(E2^2) / (N - p)

# data is overdispersed 

# other option: make season/city fixed effects and do model selection
# mod <- glmer(Coyote ~ prop_veg + rank_buff + 
#   season + City + (1|Site), data = all_dat, family = poisson)
# 


# use negative binomial instead
# can't use lme4 for this distribution in glmms though


mod2 <- glmmTMB(Deer ~ urb_pca + rank_buff + med_inc  + (1|City) + 
                  (1|season) + (1|Site), 
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

# model didn't converge so lessened random effects 
zipmod <- glmmTMB(Deer ~ urb_pca + rank_buff + med_inc  + City + 
                    season + (1|Site), 
                  data = dat2,
                  ziformula=~1,
                  family="poisson")



summary(zipmod)
E2 <- resid(zipmod, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(zipmod))   
sum(E2^2) / (N - p)

# only  slightly overdispersed 

# compare AICs 
AIC(mod1, zipmod, mod2)

# if there's overdispersion of ZIP - use a ZINB 

M4 <- zeroinfl(Counts ~ Canopy.std |
                 Canopy.std,
               dist = 'negbin',
               data = OBFL)

# check dispersion
E2 <- resid(M4, type = "pearson")
N  <- nrow(OBFL)
p  <- length(coef(M4)) + 1 # '+1' is due to theta
sum(E2^2) / (N - p)


# compare our 2 zero-inflated models 
lrtest(M3, M4)
summary(M4)


