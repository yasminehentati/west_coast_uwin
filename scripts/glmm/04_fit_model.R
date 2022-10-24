################################################################################
################# West Coast Env Health GLMM Analysis ##########################
#################      Step 4: Build Model       ###############################
################################################################################


mod <- glmer(counts ~ urbanization + income + env_health
           + (1|season) + (1|city) + (1|site), data=coyotecounts)

# other option: make season/city fixed effects and do model selection
mod <- glmer(counts ~ urbanization + income + env_health + 
               season + city + (1|site), data = coyotecounts)