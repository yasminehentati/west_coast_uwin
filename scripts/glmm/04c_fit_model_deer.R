






################################################################################
# deer counts model 

# check for zero inflation for the data set we're interested in 
100*sum(all_dat$Deer== 0)/nrow(all_dat)

# 56% 0s for the Deer - we'll try poisson and check for overdispersion 

# universal model with random effects 
mod1 <- glmmTMB(Deer ~ urb_pca + rank_buff + (1|City) + season +
                  (1|Site), offset = log(J),
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


mod2 <- glmmTMB(Deer ~ urb_pca + rank_buff + med_inc + (1|City) + season +
                  (1|Site), offset = log(J),
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
zipmod <- glmmTMB(Deer ~ urb_pca + rank_buff +(1|City) + season +
                    (1|Site), 
                  data = all_dat,
                  ziformula=~1,
                  family=poisson)



summary(zipmod)
confint(zipmod)
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
mod1 <- glmmTMB(Deer ~ med_inc + urb_pca + rank_buff + (1|season) + (1|City) +
                  (1|Site), offset = log(J), 
                data=all_dat, family = "poisson")

# check for overdispersion
summary(mod1)
E2 <- resid(mod1, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(mod1))   
sum(E2^2) / (N - p)

# data is overdispersed 
# use negative binomial instead

mod2 <- glmmTMB(Deer ~ med_inc + urb_pca + rank_buff + (1|season) + (1|City) +
                  (1|Site), offset = log(J), 
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

zipmod <- glmmTMB(Deer ~ med_inc + urb_pca + rank_buff + (1|season) + (1|City) +
                    (1|Site), offset = log(J), data = all_dat, 
                  ziformula=~1,
                  family="poisson")



summary(zipmod)
E2 <- resid(zipmod, type = "pearson")
N  <- nrow(all_dat)
p  <- length(coef(zipmod))   
sum(E2^2) / (N - p)

# good 

# compare AICs 
AIC(mod1, zipmod, mod2)

#################################################################################


##################################################################################
# plot 

mydf <- ggpredict(zipmod, terms = "urb_pca")

ggplot(mydf, aes(x, predicted)) +
  geom_line() +
  labs(
    x ="urbanization_pca", y = "richness") + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) + 
  xlim(-1,1)
confint(zipmod)

mydf2 <- ggpredict(mod1, terms = "rank_buff")

ggplot(mydf2, aes(x, predicted)) +
  geom_line() +
  labs(
    x ="environmental health index", y = "richness") + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)


mydf3 <- ggpredict(mod1, terms = c("med_inc"))
confint(zipmod)
?ggpredict
ggplot(mydf3, aes(x, predicted)) +
  geom_line() +
  labs(
    x ="median income", y = "richness") + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)
summary(zipmod)
