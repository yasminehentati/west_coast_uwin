################################################################################
################# West Coast Env Health GLMM Analysis ##########################
#################             Maps              ###############################
################################################################################

# Code adapted from Timo Grossenbacher

# packages
library(mapview)
library(lme4)
library(tidyverse) # ggplot2, dplyr, tidyr, readr, purrr, tibble
library(magrittr) # pipes
library(lintr) # code linting
library(sf) # spatial data handling
library(raster) # raster handling (needed for relief)
library(viridis) # viridis color scale
library(cowplot) # stack ggplots
library(rmarkdown)
library(sf)
library(viridis)
library(here)


################################################################################
# read in data 

# income shapefiles 
sf_map <- st_read(here("data", "income_maps", "SF_income_envhealth.shp")) %>% 
  na.omit()
la_map <- st_read(here("data", "income_maps", "LA_income_envhealth.shp")) %>%
  na.omit()
wa_map <- st_read(here("data", "income_maps", "WA_income_envhealth.shp")) %>%
  na.omit() 
mapview(la_map)
# crop WA 
wa_map <- wa_map %>% st_transform(crs = "EPSG:4326") %>% 
  st_crop(c(xmin = -122.7, xmax = -122.32, ymin = 47.13, ymax = 47.3183))

campoints_WA <- st_read(here("data", "wa_camera_points.shp"))%>% 
  st_transform(st_crs(wa_map)) %>% 
  st_crop(wa_map)

# add water
water <- st_read(here("data", "water_bodies", 
                      "DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation.shp")) 
sf_use_s2(FALSE)
water <- water %>% st_transform(st_crs(wa_map)) %>%
  st_crop(wa_map)

# crop SF 

sf_map <- sf_map %>% st_transform(crs = "EPSG:4326") %>% 
  st_crop(c(xmin = -122.34, xmax = -122.07, ymin = 37.71, ymax = 37.823))

campoints_SF <- st_read(here("data", "sf_camera_points.shp"))%>% 
  st_transform(st_crs(sf_map)) %>% 
  st_crop(sf_map)


# add water
water_SF <- st_read(here("data", "water_bodies", 
                      "bayarea_allwater.shp")) 
sf_use_s2(FALSE)
water_SF <- water_SF %>% st_transform(st_crs(sf_map)) %>%
  st_crop(sf_map)



# crop LA 


# crop LA
la_map <- la_map %>% st_transform(crs = "EPSG:4326") %>% 
  st_crop(c(xmin = -118.81, xmax = -117.608, ymin = 33.709, ymax = 34.2844))
mapview(la_map)
campoints_LA <- st_read(here("data", "la_camera_points.shp"))%>% 
  st_transform(st_crs(la_map)) %>% 
  st_crop(la_map)

mapview(campoints_LA)
# add water
water_la <- st_read(here("data", "water_bodies", 
                      "Water.shp")) 
sf_use_s2(FALSE)


mapview(water_la)
################################################################################

##### create map color theme

theme_map <- function(...) {
  theme_minimal()  +
   theme(
 #    text = element_text(family = default_font_family,
 #                        color = default_font_color),
    # remove all axes
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    # add a subtle grid
   #  panel.grid.major = element_line(color = "#dbdbd9", size = 0.2),
panel.grid.minor = element_blank(),
# background colors
plot.background = element_rect(fill = NULL,
                                color = NA),
panel.background = element_rect(fill = NULL,
                                color = NA),
legend.background = element_rect(# fill = NULL,
                                  color = NA),
# borders and margins
plot.margin = unit(c(.5, .5, .2, .5), "cm"),
panel.border = element_blank(),
panel.spacing = unit(c(-.1, 0.2, .2, 0.2), "cm"),
# titles
legend.title = element_text(size = 11),
legend.text = element_text(size = 9, hjust = 0, # color = default_font_color
                           ),
plot.title = element_text(size = 15, hjust = 0.5, # color = default_font_color
                          ),
plot.subtitle = element_text(size = 10, hjust = 0.5, # color = default_font_color,
                             margin = margin(b = -0.1,
                                             t = -0.1,
                                             l = 2,
                                             unit = "cm"),
                             debug = F),
# captions
plot.caption = element_text(size = 7,
                            hjust = .5,
                            margin = margin(t = 0.2,
                                            b = 0,
                                            unit = "cm"),
                            color = "#939184"),
...
)
}

########################

# tacoma maps 
# code adapted from Timo Grossenbacher

# univariate chloropleth 

# define number of classes
no_classes <- 6



# extract quantiles
quantiles <- wa_map %>%
  pull(estimate) %>%
  quantile(probs = seq(0, 1, length.out = no_classes + 1)) %>%
  as.vector() # to remove names of quantiles, so idx below is numeric

# here we create custom labels
labels <- imap_chr(quantiles, function(., idx){
  return(paste0(round(quantiles[idx] / 1000, 0),
                "k",
                " – ",
                round(quantiles[idx + 1] / 1000, 0),
                "k"))
})

# we need to remove the last label 
# because that would be something like "478k - NA"
labels <- labels[1:length(labels) - 1]

# here we actually create a new 
# variable on the dataset with the quantiles
wa_map %<>%
  mutate(mean_quantiles = cut(estimate,
                              breaks = quantiles,
                              labels = labels,
                              include.lowest = T))

ggplot(
  # define main data source
  data = wa_map
) +
  # use the "alpha hack" (as the "fill" aesthetic is already taken)
  scale_alpha(name = "",
              range = c(0.6, 0),
              guide = F) + # suppress legend
  # add main fill aesthetic
  # use thin white stroke for municipality borders
  geom_sf(
    mapping = aes(
      fill = mean_quantiles
    ),
    color = "white",
    size = 0.1
  ) +
  # use the Viridis color scale
  scale_fill_viridis(
    option = "magma",
    name = "Average\nincome in CHF",
    alpha = 0.8, # make fill a bit brighter
    begin = 0.1, # this option seems to be new (compared to 2016):
    # with this we can truncate the
    # color scale, so that extreme colors (very dark and very bright) are not
    # used, which makes the map a bit more aesthetic
    end = 0.9,
    discrete = T, # discrete classes, thus guide_legend instead of _colorbar
    direction = 1, # dark is lowest, yellow is highest
    guide = guide_legend(
      keyheight = unit(5, units = "mm"),
      title.position = "top",
      reverse = T # display highest income on top
    )) +
  # draw lakes in light blue
  geom_sf(
    data = water,
    fill = "#D6F1FF",
    color = "transparent"
  ) +
  # add titles
  labs(x = NULL,
       y = NULL,
       title = "Switzerland's regional income",
       subtitle = "Average yearly income in Swiss municipalities, 2015") +
  # add theme
  theme_map()

wa_map$estimate
################################################################
# bivariate map 

##### get quantiles for our bivariate data 
# create 3 buckets for gini
quantiles_gini <- wa_map %>%
  pull(Rank) %>%
  quantile(probs = seq(0, 1, length.out = 4))

# create 3 buckets for mean income
quantiles_mean <- wa_map %>%
  pull(estimate) %>%
  quantile(probs = seq(0, 1, length.out = 4))

# create color scale that encodes two variables
# red for gini and blue for mean income
# the special notation with gather is due to readibility reasons
bivariate_color_scale <- tibble(
  "3 - 3" = "#3F2949", # high inequality, high income
  "2 - 3" = "#435786",
  "1 - 3" = "#4885C1", # low inequality, high income
  "3 - 2" = "#77324C",
  "2 - 2" = "#806A8A", # medium inequality, medium income
  "1 - 2" = "#89A1C8",
  "3 - 1" = "#AE3A4E", # high inequality, low income
  "2 - 1" = "#BC7C8F",
  "1 - 1" = "#CABED0" # low inequality, low income
) %>%
  gather("group", "fill")


# cut into groups defined above and join fill
wa_map %<>%
  mutate(
    gini_quantiles = cut(
      Rank,
      breaks = quantiles_gini,
      include.lowest = TRUE
    ),
    mean_quantiles = cut(
      estimate,
      breaks = quantiles_mean,
      include.lowest = TRUE
    ),
    # by pasting the factors together as numbers we match the groups defined
    # in the tibble bivariate_color_scale
    group = paste(
      as.numeric(gini_quantiles), "-",
      as.numeric(mean_quantiles)
    )
  ) %>%
  # we now join the actual hex values per "group"
  # so each municipality knows its hex value based on the his gini and avg
  # income value
  left_join(bivariate_color_scale, by = "group")



###### plot map 
wamap <- ggplot(
  # use the same dataset as before
  data = wa_map, 
) +
 
  # use the "alpha hack" (as the "fill" aesthetic is already taken)
  scale_alpha(name = "",
              range = c(0.6, 0),
              guide = F) + # suppress legend
  # color municipalities according to their gini / income combination
  geom_sf(
    aes(
      fill = fill
    ),
    # use thin white stroke for municipalities
    color = "white",
    size = 0.1
  )  + 
  # as the sf object  has a column with name "fill" that
  # contains the literal color as hex code for each municipality, we can use
  # scale_fill_identity here
  scale_fill_identity() +
  # draw lakes in light blue
  geom_sf(
    data = water,
    fill = "#D6F1FF",
    color = "transparent"
  ) +
  # add titles
#    labs(x = NULL,
  #      y = NULL,
  #      title = "Switzerland's regional income (in-)equality",
  #      subtitle = paste0("Average yearly income and income",
  #                        " (in-)equality in Swiss municipalities, 2015")) + #  + ,
       #caption = default_caption) +
  # add the theme
  theme_map() + 
  # add points
  geom_sf(data = campoints_WA, size = 2, fill = "black")

wamap

##### # draw legend

bivariate_color_scale %<>%
  separate(group, into = c("gini", "mean"), sep = " - ") %>%
  mutate(gini = as.integer(gini),
         mean = as.integer(mean))

legend <- ggplot() +
  geom_tile(
    data = bivariate_color_scale,
    mapping = aes(
      x = gini,
      y = mean,
      fill = fill)
  ) +
  scale_fill_identity() +
  labs(x = "Higher inequality ⟶️",
       y = "Higher income ⟶️") +
  theme_map() +
  # make font small enough
  theme(
    axis.title = element_text(size = 6)
  ) +
  # quadratic tiles
  coord_fixed()

ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, -0.05, 0.075, 0.2, 0.2)
  
# combine map and legend with cowplot 





################################################################################################
### sf maps 
# univariate chloropleth 

# define number of classes
no_classes <- 6

# extract quantiles
quantiles <- sf_map %>%
  pull(estimate) %>%
  quantile(probs = seq(0, 1, length.out = no_classes + 1)) %>%
  as.vector() # to remove names of quantiles, so idx below is numeric

# here we create custom labels
labels <- imap_chr(quantiles, function(., idx){
  return(paste0(round(quantiles[idx] / 1000, 0),
                "k",
                " – ",
                round(quantiles[idx + 1] / 1000, 0),
                "k"))
})

# we need to remove the last label 
# because that would be something like "478k - NA"
labels <- labels[1:length(labels) - 1]

# here we actually create a new 
# variable on the dataset with the quantiles
sf_map %<>%
  mutate(mean_quantiles = cut(estimate,
                              breaks = quantiles,
                              labels = labels,
                              include.lowest = T))

ggplot(
  # define main data source
  data = sf_map
) +
  # use the "alpha hack" (as the "fill" aesthetic is already taken)
  scale_alpha(name = "",
              range = c(0.6, 0),
              guide = F) + # suppress legend
  # add main fill aesthetic
  # use thin white stroke for municipality borders
  geom_sf(
    mapping = aes(
      fill = mean_quantiles
    ),
    color = "white",
    size = 0.1
  ) +
  # use the Viridis color scale
  scale_fill_viridis(
    option = "magma",
    name = "Average\nincome in CHF",
    alpha = 0.8, # make fill a bit brighter
    begin = 0.1, # this option seems to be new (compared to 2016):
    # with this we can truncate the
    # color scale, so that extreme colors (very dark and very bright) are not
    # used, which makes the map a bit more aesthetic
    end = 0.9,
    discrete = T, # discrete classes, thus guide_legend instead of _colorbar
    direction = 1, # dark is lowest, yellow is highest
    guide = guide_legend(
      keyheight = unit(5, units = "mm"),
      title.position = "top",
      reverse = T # display highest income on top
    )) +
  # draw lakes in light blue
  geom_sf(
    data = water_SF,
    fill = "#D6F1FF",
    color = "transparent"
  ) +
  # add titles
  labs(x = NULL,
       y = NULL,
       title = "Switzerland's regional income",
       subtitle = "Average yearly income in Swiss municipalities, 2015") +
  # add theme
  theme_map()


################################################################
# bivariate map 

##### get quantiles for our bivariate data 
# create 3 buckets for gini
quantiles_gini <- sf_map %>%
  pull(Rank) %>%
  quantile(probs = seq(0, 1, length.out = 4))

# create 3 buckets for mean income
quantiles_mean <- sf_map %>%
  pull(estimate) %>%
  quantile(probs = seq(0, 1, length.out = 4))

# create color scale that encodes two variables
# red for gini and blue for mean income
# the special notation with gather is due to readibility reasons
bivariate_color_scale <- tibble(
  "3 - 3" = "#3F2949", # high inequality, high income
  "2 - 3" = "#435786",
  "1 - 3" = "#4885C1", # low inequality, high income
  "3 - 2" = "#77324C",
  "2 - 2" = "#806A8A", # medium inequality, medium income
  "1 - 2" = "#89A1C8",
  "3 - 1" = "#AE3A4E", # high inequality, low income
  "2 - 1" = "#BC7C8F",
  "1 - 1" = "#CABED0" # low inequality, low income
) %>%
  gather("group", "fill")


# cut into groups defined above and join fill
sf_map %<>%
  mutate(
    gini_quantiles = cut(
      Rank,
      breaks = quantiles_gini,
      include.lowest = TRUE
    ),
    mean_quantiles = cut(
      estimate,
      breaks = quantiles_mean,
      include.lowest = TRUE
    ),
    # by pasting the factors together as numbers we match the groups defined
    # in the tibble bivariate_color_scale
    group = paste(
      as.numeric(gini_quantiles), "-",
      as.numeric(mean_quantiles)
    )
  ) %>%
  # we now join the actual hex values per "group"
  # so each municipality knows its hex value based on the his gini and avg
  # income value
  left_join(bivariate_color_scale, by = "group")





#### add annotations 

annotations <- tibble(
  label = c(
    "Grey areas mean\nlow income and\nlow inequality",
    "Blue areas mean\nhigh income and\nlow inequality",
    "Violet areas mean\nhigh income and\nhigh inequality",
    "Red areas mean\nlow income and\nhigh inequality"
  ),
  arrow_from = c(
    "548921,232972", # grey
    "771356,238335", # blue
    "781136,125067", # violet
    "616348,81869" # red
  ),
  arrow_to = c(
    "622435,206784", # grey
    "712671,261998", # blue
    "786229,149597", # violet
    "602334,122674" # red
  ),
  curvature = c(
    0.2, # grey
    0.1, # blue
    -0.1, # violet
    -0.2 # red
  ),
  nudge = c(
    "-3000,0", # grey
    "3000,5000", # blue
    "0,-5000", # violet
    "3000,0" # red
  ),
  just = c(
    "1,0", # grey
    "0,1", # blue
    "0.5,1", # violet
    "0,1" # red
  )
) %>%
  separate(arrow_from, into = c("x", "y")) %>%
  separate(arrow_to, into = c("xend", "yend")) %>%
  separate(nudge, into = c("nudge_x", "nudge_y"), sep = "\\,") %>%
  separate(just, into = c("hjust", "vjust"), sep = "\\,")




###### plot map 
sfmap <- ggplot(
  # use the same dataset as before
  data = sf_map, 
) +
  
  # use the "alpha hack" (as the "fill" aesthetic is already taken)
  scale_alpha(name = "",
              range = c(0.6, 0),
              guide = F) + # suppress legend
  # color municipalities according to their gini / income combination
  geom_sf(
    aes(
      fill = fill
    ),
    # use thin white stroke for municipalities
    color = "white",
    size = 0.1
  )  + 
  # as the sf object  has a column with name "fill" that
  # contains the literal color as hex code for each municipality, we can use
  # scale_fill_identity here
  scale_fill_identity() +
  # draw lakes in light blue
  geom_sf(
    data = water_SF,
    fill = "#D6F1FF",
    color = "transparent"
  ) +
  # add titles
 #  labs(x = NULL,
#        y = NULL,
 #       title = "Switzerland's regional income (in-)equality",
 #       subtitle = paste0("Average yearly income and income",
 #                         " (in-)equality in Swiss municipalities, 2015")) + #  + ,
  #caption = default_caption) +
  # add the theme
  theme_map() + 
  # add points
  geom_sf(data = campoints_SF, size = 2, fill = "black")



##### # draw legend

bivariate_color_scale %<>%
  separate(group, into = c("gini", "mean"), sep = " - ") %>%
  mutate(gini = as.integer(gini),
         mean = as.integer(mean))

legend <- ggplot() +
  geom_tile(
    data = bivariate_color_scale,
    mapping = aes(
      x = gini,
      y = mean,
      fill = fill)
  ) +
  scale_fill_identity() +
  labs(x = "Higher environmental health burden ⟶️",
       y = "Higher income ⟶️") +
  theme_map() +
  # make font small enough
  theme(
    axis.title = element_text(size = 20)
  ) +
  # quadratic tiles
  coord_fixed()
legend
ggdraw() + 
  draw_plot(sfmap, 0, 0, 1, 1) +  
  draw_plot(legend, 0.05, 0.2, 0.2, 0.2)

###################################################################################
# la maps 

# univariate chloropleth 

# define number of classes
no_classes <- 6

# extract quantiles
quantiles <- la_map %>%
  pull(estimate) %>%
  quantile(probs = seq(0, 1, length.out = no_classes + 1)) %>%
  as.vector() # to remove names of quantiles, so idx below is numeric

# here we create custom labels
labels <- imap_chr(quantiles, function(., idx){
  return(paste0(round(quantiles[idx] / 1000, 0),
                "k",
                " – ",
                round(quantiles[idx + 1] / 1000, 0),
                "k"))
})

# we need to remove the last label 
# because that would be something like "478k - NA"
labels <- labels[1:length(labels) - 1]

# here we actually create a new 
# variable on the dataset with the quantiles
la_map %<>%
  mutate(mean_quantiles = cut(estimate,
                              breaks = quantiles,
                              labels = labels,
                              include.lowest = T))

ggplot(
  # define main data source
  data = la_map
) +
  # use the "alpha hack" (as the "fill" aesthetic is already taken)
  scale_alpha(name = "",
              range = c(0.6, 0),
              guide = F) + # suppress legend
  # add main fill aesthetic
  # use thin white stroke for municipality borders
  geom_sf(
    mapping = aes(
      fill = mean_quantiles
    ),
    color = "white",
    size = 0.1
  ) +
  # use the Viridis color scale
  scale_fill_viridis(
    option = "magma",
    name = "Average\nincome in CHF",
    alpha = 0.8, # make fill a bit brighter
    begin = 0.1, # this option seems to be new (compared to 2016):
    # with this we can truncate the
    # color scale, so that extreme colors (very dark and very bright) are not
    # used, which makes the map a bit more aesthetic
    end = 0.9,
    discrete = T, # discrete classes, thus guide_legend instead of _colorbar
    direction = 1, # dark is lowest, yellow is highest
    guide = guide_legend(
      keyheight = unit(5, units = "mm"),
      title.position = "top",
      reverse = T # display highest income on top
    )) +
  # draw lakes in light blue
  # add titles
#   labs(x = NULL,
  #      y = NULL,
 #       title = "Switzerland's regional income",
 #       subtitle = "Average yearly income in Swiss municipalities, 2015") +
  # add theme
  theme_map()


################################################################
# bivariate map 

##### get quantiles for our bivariate data 
# create 3 buckets for gini
quantiles_gini <- la_map %>%
  pull(Rank) %>%
  quantile(probs = seq(0, 1, length.out = 4))

# create 3 buckets for mean income
quantiles_mean <- la_map %>%
  pull(estimate) %>%
  quantile(probs = seq(0, 1, length.out = 4))

# create color scale that encodes two variables
# red for gini and blue for mean income
# the special notation with gather is due to readibility reasons
bivariate_color_scale <- tibble(
  "3 - 3" = "#3F2949", # high inequality, high income
  "2 - 3" = "#435786",
  "1 - 3" = "#4885C1", # low inequality, high income
  "3 - 2" = "#77324C",
  "2 - 2" = "#806A8A", # medium inequality, medium income
  "1 - 2" = "#89A1C8",
  "3 - 1" = "#AE3A4E", # high inequality, low income
  "2 - 1" = "#BC7C8F",
  "1 - 1" = "#CABED0" # low inequality, low income
) %>%
  gather("group", "fill")


# cut into groups defined above and join fill
la_map %<>%
  mutate(
    gini_quantiles = cut(
      Rank,
      breaks = quantiles_gini,
      include.lowest = TRUE
    ),
    mean_quantiles = cut(
      estimate,
      breaks = quantiles_mean,
      include.lowest = TRUE
    ),
    # by pasting the factors together as numbers we match the groups defined
    # in the tibble bivariate_color_scale
    group = paste(
      as.numeric(gini_quantiles), "-",
      as.numeric(mean_quantiles)
    )
  ) %>%
  # we now join the actual hex values per "group"
  # so each municipality knows its hex value based on the his gini and avg
  # income value
  left_join(bivariate_color_scale, by = "group")




###### plot map 
lamap <- ggplot(
  # use the same dataset as before
  data = la_map, 
) +
  
  # use the "alpha hack" (as the "fill" aesthetic is already taken)
  scale_alpha(name = "",
              range = c(0.6, 0),
              guide = F) + # suppress legend
  # color municipalities according to their gini / income combination
  geom_sf(
    aes(
      fill = fill
    ),
    # use thin white stroke for municipalities
    color = "white",
    size = 0.1
  )  + 
  # as the sf object  has a column with name "fill" that
  # contains the literal color as hex code for each municipality, we can use
  # scale_fill_identity here
  scale_fill_identity() +
  
  # add titles
  #  labs(x = NULL,
  #        y = NULL,
  #       title = "Switzerland's regional income (in-)equality",
  #       subtitle = paste0("Average yearly income and income",
  #                         " (in-)equality in Swiss municipalities, 2015")) + #  + ,
  #caption = default_caption) +
  # add the theme
  theme_map() + 
  # add points
  geom_sf(data = campoints_LA, size = 2, fill = "black")



##### # draw legend

bivariate_color_scale %<>%
  separate(group, into = c("gini", "mean"), sep = " - ") %>%
  mutate(gini = as.integer(gini),
         mean = as.integer(mean))

legend <- ggplot() +
  geom_tile(
    data = bivariate_color_scale,
    mapping = aes(
      x = gini,
      y = mean,
      fill = fill)
  ) +
  scale_fill_identity() +
  labs(x = "Higher environmental health burden ⟶️",
       y = "Higher income ⟶️") +
  theme_map() +
  # make font small enough
  theme(
    axis.title = element_text(size = 20)
  ) +
  # quadratic tiles
  coord_fixed()
legend
ggdraw() + 
  draw_plot(sfmap, 0, 0, 1, 1) +  
  draw_plot(legend, 0.05, 0.2, 0.2, 0.2)

lamap
mapview(la_map)
########
# plot all maps and legend
allplot <- plot_grid(sfmap, wamap)
plot_grid(allplot, legend, ncol = 1, rel_heights = c(1, .1))
