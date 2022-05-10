install.packages('devtools')
# library(devtools)
# install_github("rzepkka/ggseg3D")

library(ggseg3d)

aseg_data <- read.csv(file = 'data/aseg_R_slider_Subtype 0.csv')

# ========== NEW GGSEG_3D =======================================

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
  
  pal.colours <- custom_palette(palette)
  
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
  
  p
}

# ========== COLORBAR =======================================
custom_palette <- function(palette){
  
  if(is.null(palette)){
    palette = c("skyblue", "dodgerblue")
  }
  
  if(!is.null(names(palette))){
    pal.colours <- names(palette)
    pal.values <- unname(palette)
    pal.norm <- range_norm(pal.values)
  }else{
    pal.colours <- palette
    pal.norm <- seq(0,1, length.out = length(pal.colours))
    pal.values <- seq(0,1, length.out = length(pal.colours))
  }
  
  # Might be a single colour
  pal.colours = if(length(palette) == 1){
    # If a single colour, dummy create a second
    # palette row for interpolation
    data.frame(values = c(pal.values,pal.values+1),
               norm = c(0, 1 ),
               orig = c(pal.colours,pal.colours),
               stringsAsFactors = F)
  }else{
    data.frame(values = pal.values,
               norm = pal.norm,
               orig = pal.colours,
               stringsAsFactors = F)
  }
  
  pal.colours$hex <- gradient_n_pal(
    colours = pal.colours$orig,
    values = pal.colours$values,
    space = "Lab")(pal.colours$values)
  
  pal.colours
  
}

# =========== ADD CUSTOM FUNCTIONS TO ENVIRONMENT ==============================

environment(custom_ggseg) <- asNamespace('ggseg3d')
environment(custom_palette) <- asNamespace('ggseg3d')



aseg = custom_ggseg(.data = aseg_data, 
               atlas = aseg_3d, 
               colour = "p", 
               palette = colors,
               text = "p", 
               # na.alpha= .5,
               # options.legend = list(title=list(text="Subcortical")),
               scene = 'scene'
)

aseg
