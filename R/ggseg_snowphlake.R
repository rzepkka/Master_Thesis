library(ggseg3d)
library(ggseg)
library(ggplot2)
library(dplyr)
library(tidyr)
library(htmltools)
library(htmlwidgets)

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

dk_data <- read.csv(file = 'data/dk_R_subtype0.csv')
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

# save_html(p, 'plot_htmltools.html', background = "white", libdir = "lib", lang = "en")
# saveWidget(p, "plot_htmlwidgets.html", selfcontained = F, libdir = "lib")

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
        palette = colors,
        text = "p",
        options.legend = list(title=list(text=""))) %>% 
  pan_camera("left medial")


# ============= ASEG ========================================================================================

aseg_data <- read.csv(file = 'data/aseg_R_subtype0.csv')
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

#=========== INPUTS TO 3D SUBPLOTS ========================================================================================

# LEFT HEMISPHERE
dk_left = ggseg3d(.data = dk_data,
                  atlas = dk_3d,
                  hemisphere = 'left',
                  colour = "p",
                  palette = colors,
                  text = "p",
                  options.legend = list(title=list(text="Left")),
                  scene = 'scene')
dk_left

# RIGHT HEMISPHERE
dk_right = ggseg3d(.data = dk_data,
                   atlas = dk_3d,
                   hemisphere = 'right',
                   colour = "p",
                   palette = colors,
                   text = "p",
                   options.legend = list(title=list(text="Right")),
                   scene = 'scene2')
dk_right

# WHOLE BRAIN
dk_whole = ggseg3d(.data = dk_data,
                   atlas = dk_3d,
                   hemisphere = c('left','right'),
                   colour = "p",
                   palette = colors,
                   text = "p",
                   options.legend = list(title=list(text="Cortical")),
                   scene = 'scene3')
dk_whole

# ASEG
aseg = ggseg3d(.data = aseg_data, 
               atlas = aseg_3d, 
               colour = "p", 
               palette = colors,
               text = "p", 
               options.legend = list(title=list(text="Subcortical")),
               scene = 'scene4')

aseg

# ========= 3D SUBPLOTS =============================================================================================================================================

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
             scene = 'scene')
dk

aseg = ggseg3d(.data = aseg_data, 
               atlas = aseg_3d, 
               colour = "p", 
               palette = colors,
               text = "p", 
               # na.alpha= .5,
               options.legend = list(title=list(text="Subcortical"),
                                     colorbar=list(limits=c(0,25))),
               scene = 'scene2',
)

aseg

colorbar(dk, limits=c(0,25))


fig2 <- subplot(dk, aseg)
fig2 <- fig2 %>% layout(title = "3D Visualization of the disease timeline",
                        scene = list(domain=list(x=c(0,1),y=c(0.5,1)),
                                     aspectmode='auto',
                                     xaxis=list(backgroundcolor='white')),
                        scene2 = list(domain=list(x=c(0,1),y=c(0,0.5)),
                                      aspectmode='auto'
                        ))

fig2

# p1.scene = 'scene1'
# p2.scene = 'scene2'

saveWidget(fig2, "2_plots.html", selfcontained = F, libdir = "lib")


ggseg3d(show.legend = TRUE, hemisphere = "left", limits=c(0,25))

