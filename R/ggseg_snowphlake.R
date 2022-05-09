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

dk_data <- read.csv(file = 'data/dk_R_Subtype 0.csv')
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
dk_data <- read.csv(file = 'data/dk_R_slider_Subtype 4.csv')
aseg_data <- read.csv(file = 'data/aseg_R_slider_Subtype 4.csv')

dk_data
aseg_data

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
               options.legend = list(title=list(text="Subcortical"),dtick=0.1,
                                     tickformatstops=list(dtickrange=c(0,1))),
               scene = 'scene2'
)


aseg
# ========= 3D SUBPLOTS =============================================================================================================================================

fig2 <- subplot(dk, aseg)
fig2 <- fig2 %>% layout(title = "Subtype 4 - 3D Visualization of the disease timeline",
                        scene = list(domain=list(x=c(0,1),y=c(0.5,1)),
                                     aspectmode='auto',
                                     xaxis=list(backgroundcolor='white')
                                     ),
                        scene2 = list(domain=list(x=c(0,1),y=c(0,0.5)),
                                      aspectmode='auto'
                        )) 
fig2 

# SAVE 3D PLOTS TO LOAD THEM INTO STREAMLIT 

saveWidget(fig2, "html_3D/slider/Subtype 4.html", selfcontained = F, libdir = "lib")


# =============================

p <- plot_ly(mtcars, x = ~wt, y = ~mpg, color = ~cyl)

# pass any colorbar attribute -- 
# https://plotly.com/r/reference/#scatter-marker-colorbar
colorbar(p, len = 0.5)

# Expand the limits of the colorbar
colorbar(p, limits = c(0, 20))
# values outside the colorbar limits are considered "missing"
colorbar(p, limits = c(5, 6))

# also works on colorbars generated via a z value
corr <- cor(diamonds[vapply(diamonds, is.numeric, logical(1))])
plot_ly(x = rownames(corr), y = colnames(corr), z = corr) %>%
  add_heatmap() %>%
  colorbar(limits = c(-1, 1))

