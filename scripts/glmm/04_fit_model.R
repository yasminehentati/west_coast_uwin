################################################################################
################# West Coast Env Health GLMM Analysis ##########################
#################      Step 4: Build Model       ###############################
################################################################################



library(readr)
library(here)
library(lme4)

install.packages("R2admb")
install.packages("glmmADMB", 
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                         getOption("repos")),
                 type="source")
library(glmmADMB)

# read in data 
LA_dat <- read_csv(here("data", "LA_counts_covs_10-25-22.csv"))
SF_dat <- read_csv(here("data", "SF_counts_covs_10-25-22.csv"))
WA_dat <- read_csv(here("data", "WA_counts_covs_10-25-22.csv"))

# combine data
cal_dat <- rbind(LA_dat, SF_dat)
all_dat <- rbind(cal_dat, WA_dat)


head(all_dat)

# make things factors
all_dat$City <- as.factor(all_dat$City)
all_dat$season <- as.factor(all_dat$season)
all_dat$Site <- as.factor(all_dat$Site)

################################################################################
# coyote counts model 

# check for zero inflation for the data set we're interested in 
100*sum(all_dat$Coyote == 0)/nrow(all_dat)

# 25% 0s for the coyotes - we'll try poisson and check for overdispersion 

# universal model with random effects 
mod1 <- glmer(Coyote ~ prop_veg + rank_buff
           + season + (1|City) + (1|Site), 
           data=all_dat, family = poisson)

# other option: make season/city fixed effects and do model selection
mod <- glmer(Coyote ~ prop_veg + rank_buff + 
                season + City + (1|Site), data = all_dat, family = poisson)


# check for overdispersion
summary(mod1)
E2 <- resid(mod1, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(mod1))   
sum(E2^2) / (N - p)

# if overdispersed - use negative binomial instead
# can't use lme4 for this distribution


mod2 <- glm.nb(Coyote ~ prop_veg + rank_buff + 
                 season + City + (1|Site),
             data = all_dat)

# check for underdispersion
E2 <- resid(M2, type = "pearson")
N  <- nrow(OBFL)
p  <- length(coef(M2)) + 1  # '+1' is for variance parameter in NB
sum(E2^2) / (N - p)

# if underdispersed - try zero-inflated poisson
# here we can choose different predictors for count and overdisp - 
# because we might expect different variables to drive presence/absence
# than the ones that drive total # of individuals (count) 

dat2 <- na.omit(all_dat)

zipmod <- glmmadmb(Coyote~ prop_veg + rank_buff + 
                             season + City + (1|Site),
                           data = all_dat,
                        zeroInflation=TRUE,
                        family="poisson")

summary(M3)

# if there's overdispersion - use a ZINB 

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