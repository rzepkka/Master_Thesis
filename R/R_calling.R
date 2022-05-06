# print("Hello from R script!")

library(ggseg3d)
library(ggseg)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)

library(htmlwidgets)

dk_data <- read.csv(file = 'data/dk_R_subtype0.csv')
# 
aseg_data <- read.csv(file = 'data/aseg_R_subtype0.csv')

# aseg = ggseg3d(.data = aseg_data,
#         atlas = aseg_3d,
#         colour = "p",
#         palette = colors,
#         text = "p",
#         na.alpha= .5,
#         options.legend = list(title=list(text=""))
# )
# 
# aseg

# p = ggseg3d(.data = dk_data,
#             atlas = dk_3d,
#             hemisphere = 'left',
#             colour = "p",
#             palette = colors,
#             text = "p")

# ggseg3d() %>%
#   remove_axes() %>%
#   pan_camera("right medial")



# someData = dk_3d %>%
#   filter(surf == "inflated" & hemi == "right") %>%
#   unnest(ggseg_3d) %>%
#   ungroup() %>%
#   select(region) %>%
#   na.omit() %>%
#   mutate(p = sample(seq(0,.5, length.out = 100 ), nrow(.)) %>%
#            round(2))


# p1 = ggseg3d(.data = someData,
#         atlas = dk_3d,
#         colour = "p", text = "p",
#         scene = 'scene1')
# 
# p1

colors = c("indianred4",'indianred2','coral1','lightpink1','mistyrose1')


dk_left = ggseg3d(.data = dk_data,
        atlas = dk_3d,
        hemisphere = 'left',
        colour = "p",
        palette = colors,
        text = "p",
        options.legend = list(title=list(text="Left")),
        scene = 'scene')
dk_left

dk_right = ggseg3d(.data = dk_data,
                  atlas = dk_3d,
                  hemisphere = 'right',
                  colour = "p",
                  palette = colors,
                  text = "p",,
                  options.legend = list(title=list(text="Right")),
                  scene = 'scene2')
dk_right

dk_whole = ggseg3d(.data = dk_data,
                   atlas = dk_3d,
                   hemisphere = c('left','right'),
                   colour = "p",
                   palette = colors,
                   text = "p",
                   options.legend = list(title=list(text="Cortical")),
                   scene = 'scene3')
dk_whole

aseg = ggseg3d(.data = aseg_data, 
        atlas = aseg_3d, 
        colour = "p", 
        palette = colors,
        text = "p", 
        options.legend = list(title=list(text="Subcortical")),
        scene = 'scene4')

aseg

# ========= 3D SUBPLOTS ===============

# subplot and define scene
fig1 <- subplot(dk_left, dk_right, dk_whole, aseg)

# 4-subplots
fig1 <- fig1 %>% layout(title = "3D Subplots",
                      scene = list(text ='Right hemishphere',domain=list(x=c(0,0.5),y=c(0.5,1)),
                                   aspectmode='auto', pan_camera ="left medial"),
                      scene2 = list(domain=list(x=c(0.5,1),y=c(0.5,1)),
                                    aspectmode='auto'),
                      scene3 = list(domain=list(x=c(0,0.5),y=c(0,0.5)),
                                    aspectmode='auto'),
                      scene4 = list(domain=list(x=c(0.5,1),y=c(0,0.5)),
                                    aspectmode='auto'))

fig1

saveWidget(fig1, "4_plots.html", selfcontained = F, libdir = "lib")


# 2 subplots

dk = ggseg3d(.data = dk_data,
                   atlas = dk_3d,
                   hemisphere = c('left','right'),
                   colour = "p",
                   palette = colors,
                   text = "p",
                   options.legend = list(title=list(text="Cortical")),
                   scene = 'scene2')
dk



aseg = ggseg3d(.data = aseg_data, 
               atlas = aseg_3d, 
               colour = "p", 
               palette = colors,
               text = "p", 
               # na.alpha= .5,
               options.legend = list(title=list(text="Subcortical")),
               scene = 'scene'
)

aseg


fig2 <- subplot(dk, aseg)
fig2 <- fig2 %>% layout(title = "3D Subplots",
                       scene = list(domain=list(x=c(0,1),y=c(0.5,1)),
                                     aspectmode='auto',
                                    xaxis=list(backgroundcolor='white')),
                       scene2 = list(domain=list(x=c(0,1),y=c(0,0.5)),
                                     aspectmode='auto'
                                     # ,
                                     # xaxis=list(backgroundcolor='white'),
                                     # yaxis=list(backgroundcolor='white'),
                                     # zaxis=list(backgroundcolor='white')
                                     ))

fig2

# p1.scene = 'scene1'
# p2.scene = 'scene2'

saveWidget(fig2, "2_plots.html", selfcontained = F, libdir = "lib")



