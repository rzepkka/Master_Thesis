import pandas as pd
import numpy as np
import pickle
import json
import warnings
warnings.filterwarnings("ignore")

import rpy2
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr, data
import rpy2.robjects as robjects
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri

from mapping_3D import dk_3D, aseg_3D, save_subtype_data


# =================== 1. LOAD PICKLE FILE ========================================================================================================================================================
read_input_file = open('data/EDADS_subtype_timelines_agecorrected_opt.pickle','rb')
load_inputs = pickle.load(read_input_file)
read_input_file.close()

T, S, Sboot = load_inputs

# =================== 2. GET SUBTYPE LABELS=====================================================================================================================================
def get_labels(S):
    unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
    subtype_labels = []
    for i in range(len(unique_subtypes)):
        subtype_labels.append('Subtype '+str(int(unique_subtypes[i])))        
    return subtype_labels

labels = get_labels(S=S)


# =================== 3. LOAD JSON FILES FOR BRAIN REGION MAPPINGS ===============================================================================================
f = open('data/DK_3D_combined.json')
map_dk = json.load(f)
f.close()

f = open('data/ASEG_3D_combined.json')
map_aseg = json.load(f)
f.close()


# =================== 4. RUN MAPPING FUNCTIONS - dk_3D() AND aseg_3D() FOR EACH Subtype AND SAVE DATA IN temp_folder
save_subtype_data(T, S, map_dk, map_aseg)


# =================== 5. SET UP R ========================================================================================================================

# import libraries
utils = importr('utils')
base = importr("base")
htmlwidgets = importr('htmlwidgets')
utils.chooseCRANmirror(ind=1)

# import libraries, set colors
robjects.r('''
	library(ggseg3d)
	library(ggseg)
	library(ggplot2)
	library(dplyr)
	library(tidyr)
	library(htmltools)
	library(htmlwidgets)
	library(plotly)

	options(warn=-1)

    colors = c("indianred4",'indianred2','coral1','lightpink1','mistyrose1')   
    ''')

# set colors
robjects.r('''
    colors = c("indianred4",'indianred2','coral1','lightpink1','mistyrose1')   
    ''')

# Custom ggseg3d() function and assign it to ggseg3d namespace
robjects.r('''
custom_ggseg <- function(.data=NULL, atlas="dk_3d",
                    surface = "LCBC", hemisphere = c("right","subcort"),
                    label = "region", text = NULL, colour = "colour",
                    palette = NULL, na.colour = "darkgrey", na.alpha = 1,
                    show.legend = TRUE, options.legend = NULL, ...) {
  
  # Grab the atlas, even if it has been provided as character string
  atlas3d = get_atlas(atlas,
                      surface = surface,
                      hemisphere = hemisphere)
  
  # If data has been supplied, merge it
  if(!is.null(.data)){
    atlas3d <- data_merge(.data, atlas3d)
  }
  
  pal.colours <- get_palette(palette)
  
  # If colour column is numeric, calculate the gradient
  if(is.numeric(unlist(atlas3d[,colour]))){
    
    if(is.null(names(palette))){
      pal.colours$values <- seq(0,1,length.out = nrow(pal.colours))
    }
    
    atlas3d$new_col = gradient_n_pal(pal.colours$orig, pal.colours$values,"Lab")(
      unlist(atlas3d[,colour]))
    fill = "new_col"
    
  }else{
    fill = colour
  }
  
  # initiate plot
  p = plotly::plot_ly()
  
  # add one trace per file inputed
  for(tt in 1:nrow(atlas3d)){
    
    col = rep(unlist(atlas3d[tt, fill]), nrow(atlas3d$mesh[[tt]]$faces))
    
    col = ifelse(is.na(col), na.colour, col)
    
    op = ifelse(is.na(unlist(atlas3d[tt, fill])), na.alpha, 1)
    
    txt = if(is.null(text)){
      text
    }else{
      paste0(text, ": ", unlist(atlas3d[tt, text]))
    }
    
    p = plotly::add_trace(p,
                          x = atlas3d$mesh[[tt]]$vertices$x,
                          y = atlas3d$mesh[[tt]]$vertices$y,
                          z = atlas3d$mesh[[tt]]$vertices$z,
                          
                          i = atlas3d$mesh[[tt]]$faces$i-1,
                          j = atlas3d$mesh[[tt]]$faces$j-1,
                          k = atlas3d$mesh[[tt]]$faces$k-1,
                          
                          facecolor = col,
                          type = "mesh3d",
                          text = txt,
                          showscale = FALSE,
                          opacity = op,
                          name = unlist(atlas3d[tt, label]),
                          ...
    )
  }
  
  # work around to get legend
  if(show.legend & is.numeric(unlist(atlas3d[,colour]))){
    
    dt_leg <- dplyr::mutate(pal.colours,
                            x = 0, y = 0, z = 0)
    p <- plotly::add_trace(p, data = dt_leg,
                           x = ~ x, y = ~ y, z = ~ z,
                           
                           intensity =  ~ values,
                           colorscale =  unname(dt_leg[,c("norm", "hex")]),
                           type = "mesh3d",
                           colorbar = options.legend
    )
  }
  
  return(p)
}

environment(custom_ggseg) <- asNamespace('ggseg3d')

''')

# Save html files for each subtype
for i in range(len(labels)):
    robjects.globalenv['i'] = i
        
    robjects.r('''
    
    dk_data <- read.csv(file = paste('/Users/macos/Documents/GitHub/Master_Thesis/temp_folder/dk_R_Subtype ',i,'.csv', sep=''))
    aseg_data <- read.csv(file = paste('/Users/macos/Documents/GitHub/Master_Thesis/temp_folder/aseg_R_Subtype ',i,'.csv', sep=''))   

    ''')
    
    
    output = robjects.r('''
      
    dk = custom_ggseg(.data = dk_data,
             atlas = dk_3d,
             hemisphere = c('left','right'),
             colour = "p",
             palette = colors,
             text = "p",
             options.legend = list(title=list(text="Cortical"),dtick=0.1,
                              tickformatstops=list(dtickrange=c(0,1))),
             scene = 'scene')

    aseg = custom_ggseg(.data = aseg_data, 
                   atlas = aseg_3d, 
                   colour = "p", 
                   palette = colors,
                   text = "p", 
                   options.legend = list(title=list(text="Subcortical"),dtick=0.1,
                                         tickformatstops=list(dtickrange=c(0,1))),
                   scene = 'scene2'
    )
    
    fig <- subplot(dk, aseg)
    fig <- fig %>% layout(title = paste('Subtype', i),
                            scene = list(domain=list(x=c(0,1),y=c(0.5,1)),
                                         aspectmode='auto',
                                         xaxis=list(backgroundcolor='white')
                            ),
                            scene2 = list(domain=list(x=c(0,1),y=c(0,0.5)),
                                          aspectmode='auto'
                            )) 
                              
    fig
   
   ''')
    
    
    htmlwidgets.saveWidget(output, f"temp_folder/Subtype {i}.html", selfcontained = False, libdir = f"sub {i}")


print('All HTML files successfully saved in: /temp_folder')
print('Run Streamlit APP by calling ">> streamlit run app.py" from command line')




















