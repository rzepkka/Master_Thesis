# print("Hello from R script!")

library(ggseg3d)
library(ggseg)
library(ggplot2)
library(dplyr)
library(tidyr)

# install.packages('caret')

dk_data <- read.csv(file = 'dk_R_subtype0.csv')
# 
# aseg_data <- read.csv(file = 'aseg_R_subtype0.csv')


# p = ggseg3d(.data = aseg_data, 
#         atlas = aseg_3d, 
#         colour = "p", 
#         palette = colors,
#         text = "p", 
#         na.alpha= .5,
#         options.legend = list(title=list(text=""))
# )

# p = ggseg3d(.data = dk_data, 
#             atlas = dk_3d,
#             hemisphere = 'left',
#             colour = "p", 
#             palette = colors,
#             text = "p")

# ggseg3d() %>%
#   remove_axes() %>%
#   pan_camera("right medial")

someData = dk_3d %>%
  filter(surf == "inflated" & hemi == "right") %>%
  unnest(ggseg_3d) %>%
  ungroup() %>%
  select(region) %>%
  na.omit() %>%
  mutate(p = sample(seq(0,.5, length.out = 100 ), nrow(.)) %>%
           round(2))


p = ggseg3d(.data = someData,
        atlas = dk_3d,
        colour = "p", text = "p") %>%
  pan_camera("right medial")

p
# 
# p = ggseg3d(.data = dk_data,
#         atlas = dk_3d,
#         hemisphere = 'left',
#         colour = "p",
#         text = "p")
# 
# typeof(p)
# 
# typeof(someData)
# 
# 
# aseg_3d



# ================================================================================================
# library(plotly)
# 
# fig <- plot_ly(
#   x = c(0, 1, 2, 0),
#   y = c(0, 0, 1, 2),
#   z = c(0, 2, 0, 1),
#   i = c(0, 0, 0, 1),
#   j = c(1, 2, 3, 2),
#   k = c(2, 3, 1, 3),
#   facecolor = toRGB(viridisLite::viridis(4))
# )
# 
# typeof(fig)
# 
# 
# 



