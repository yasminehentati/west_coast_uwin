




all_dat_veg
################################################################################
# species richness


# universal model with random effects )
mod1 <- glmmTMB(richness ~ urb_pca + rank_buff + (1|City),
                data=all_dat_veg, family = "poisson")

summary(mod1)

# by city
all_dat_veg_tawa  <- all_dat_veg %>% dplyr::filter(City == "tawa")

class(all_dat_veg_tawa$urb_pca)

mod2 <- glm(richness ~ urb_pca + rank_buff + med_inc,
            data=all_dat_veg_tawa , family = "poisson")


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

################################################################################
# shannon diversity

# universal model with random effects 
mod1 <- glmmTMB(shannon.di ~  urb_pca + rank_buff + (1|City),
                data=all_dat_veg, family = "gaussian")

summary(mod1)

simulateResiduals(mod1,plot=T)



##################################################################################
# plot 

mydf <- ggpredict(mod1, terms = "urb_pca")

ggplot(mydf, aes(x, predicted)) +
  geom_line() +
  labs(
    x ="urbanization_pca", y = "shannon diversity") + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)


mydf2 <- ggpredict(mod1, terms = "rank_buff")

ggplot(mydf2, aes(x, predicted)) +
  geom_line() +
  labs(
    x ="environmental health index", y = "shannon diversity") + 
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
