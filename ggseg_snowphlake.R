library(ggseg3d)
library(ggseg)
library(ggplot2)
library(dplyr)
library(tidyr)
# Enable this universe
options(repos = c(
  ggseg = 'https://ggseg.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))

# atlas data
dk_3d

filtered_dk =dk_3d %>% 
  filter(surf == "inflated" & hemi == "right") %>% 
  unnest(cols = ggseg_3d)

filtered_dk$region

aseg_3d

filtered_aseg = aseg_3d %>% 
  unnest(cols = ggseg_3d)

filtered_aseg

# ============= DK-atlas ========================================================================================

dk_data <- read.csv(file = 'dk_R_subtype0.csv')
dk_data

# DK for subtype 0

# make a change for better sequentiol HEX values
colors = c("indianred4",'indianred2','coral1','lightpink1','mistyrose1')

# LEFT HEMISPHERE
p = ggseg3d(.data = dk_data, 
        atlas = dk_3d,
        hemisphere = 'left',
        colour = "p", 
        palette = colors,
        text = "p") %>% 
  pan_camera("right medial")

p

# RIGHT HEMISPHERE
ggseg3d(.data = dk_data, 
        atlas = dk_3d,
        hemisphere = 'right',
        colour = "p", 
        palette = colors,
        text = "p") %>% 
  pan_camera("left medial")

# WHOLE BRAIN
ggseg3d(.data = dk_data, 
        atlas = dk_3d,
        hemisphere = c('left','right'),
        colour = "p", 
        # palette = colors,
        text = "p",
        options.legend = list(title=list(text=""))) %>% 
  pan_camera("left medial")


# ============= ASEG ========================================================================================

aseg_data <- read.csv(file = 'aseg_R_subtype0.csv')
aseg_data

# without point of reference
ggseg3d(.data = aseg_data, 
        atlas = aseg_3d, 
        colour = "p", 
        palette = colors,
        text = "p", 
        na.alpha= .5,
        options.legend = list(title=list(text=""))
  )

# added glassbrain
ggseg3d(.data = aseg_data, 
        atlas = aseg_3d, 
        colour = "p", 
        palette = colors,
        text = "p", 
        na.alpha= .5,
        options.legend = list(title=list(text=""))) %>% 
  add_glassbrain()


