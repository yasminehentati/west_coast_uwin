################################################################################
################# West Coast Env Health GLMM Analysis ##########################
#################             Plots              ###############################
################################################################################


# run "04_fit_model" before running this code 

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


################################################################################
# read in data 

# income shapefiles 
sf_map <- st_read(here("data", "income_maps", "sf_bay_med_income.shp"))
la_map <- st_read(here("data", "income_maps", "la_orange_county_med_income.shp"))
wa_map <- st_read(here("data", "income_maps", "seatac_med_income.shp"))

# environmental health data 
env_data <- read_csv(here("data", "env_health_ranks_all.csv"))

# subset out env data by area 
waenv <- env_data %>% dplyr::filter(substr(GEOID, 1, 5) 
                                 %in% c("53053", "53033"))

sfenv <- env_data %>% dplyr::filter(substr(GEOID, 1, 5) 
                                    %in% c("06055", "06041", 
                                           "06095", "06013",
                                           "06001", "06075",
                                           "06081", "06085"))

laenv <- env_data %>% dplyr::filter(substr(GEOID, 1, 5) 
                                    %in% c("06037", "06059"))
                  
# merge env health data to shapefile                  
sf_map <- sf_map %>% merge(sfenv, by = "GEOID", .keep_all_x=TRUE)
la_map <- la_map %>% merge(laenv, by = "GEOID", .keep_all_x=TRUE)
wa_map <- wa_map %>% merge(waenv, by = "GEOID", .keep_all_x=TRUE)


################################################################################

# maps 
# code adapted from Timo Grossenbacher

##### create map color theme

theme_map <- function(...) {
  theme_minimal() +
  theme(
    text = element_text(family = default_font_family,
                        color = default_font_color),
    # remove all axes
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    # add a subtle grid
    panel.grid.major = element_line(color = "#dbdbd9", size = 0.2),
panel.grid.minor = element_blank(),
# background colors
plot.background = element_rect(fill = default_background_color,
                               color = NA),
panel.background = element_rect(fill = default_background_color,
                                color = NA),
legend.background = element_rect(fill = default_background_color,
                                 color = NA),
# borders and margins
plot.margin = unit(c(.5, .5, .2, .5), "cm"),
panel.border = element_blank(),
panel.spacing = unit(c(-.1, 0.2, .2, 0.2), "cm"),
# titles
legend.title = element_text(size = 11),
legend.text = element_text(size = 9, hjust = 0,
                           color = default_font_color),
plot.title = element_text(size = 15, hjust = 0.5,
                          color = default_font_color),
plot.subtitle = element_text(size = 10, hjust = 0.5,
                             color = default_font_color,
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





##### get quantiles for our bivariate data 

# cut into groups defined above and join fill
municipality_prod_geo %<>%
  mutate(
    gini_quantiles = cut(
      gini,
      breaks = quantiles_gini,
      include.lowest = TRUE
    ),
    mean_quantiles = cut(
      mean,
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
map <- ggplot(
  # use the same dataset as before
  data = municipality_prod_geo
) +
  # first: draw the relief
  geom_raster(
    data = relief,
    aes(
      x = x,
      y = y,
      alpha = value
    )
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
  ) +
  # as the sf object municipality_prod_geo has a column with name "fill" that
  # contains the literal color as hex code for each municipality, we can use
  # scale_fill_identity here
  scale_fill_identity() +
  # use thicker white stroke for cantons
  geom_sf(
    data = canton_geo,
    fill = "transparent",
    color = "white",
    size = 0.5
  ) +
  # draw lakes in light blue
  geom_sf(
    data = lake_geo,
    fill = "#D6F1FF",
    color = "transparent"
  ) +
  # add titles
  labs(x = NULL,
       y = NULL,
       title = "Switzerland's regional income (in-)equality",
       subtitle = paste0("Average yearly income and income",
                         " (in-)equality in Swiss municipalities, 2015"),
       caption = default_caption) +
  # add the theme
  theme_map()

# add annotations one by one by walking over the annotations data frame
# this is necessary because we cannot define nudge_x, nudge_y and curvature
# in the aes in a data-driven way like as with x and y, for example
annotations %>%
  pwalk(function(...) {
    # collect all values in the row in a one-rowed data frame
    current <- tibble(...)
    
    # convert all columns from x to vjust to numeric
    # as pwalk apparently turns everything into a character (why???)
    current %<>%
      mutate_at(vars(x:vjust), as.numeric)
    
    # update the plot object with global assignment
    map <<- map +
      # for each annotation, add an arrow
      geom_curve(
        data = current,
        aes(
          x = x,
          xend = xend,
          y = y,
          yend = yend
        ),
        # that's the whole point of doing this loop:
        curvature = current %>% pull(curvature),
        size = 0.2,
        arrow = arrow(
          length = unit(0.005, "npc")
        )
      ) +
      # for each annotation, add a label
      geom_text(
        data = current,
        aes(
          x = x,
          y = y,
          label = label,
          hjust = hjust,
          vjust = vjust
        ),
        # that's the whole point of doing this loop:
        nudge_x = current %>% pull(nudge_x),
        nudge_y = current %>% pull(nudge_y),
        # other styles
        family = default_font_family,
        color = default_font_color,
        size = 3
      )
  })


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


# combine map and legend with cowplot 
ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.05, 0.075, 0.2, 0.2)