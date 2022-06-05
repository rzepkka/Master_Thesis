#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# my_function <- function(name){
#   print(paste("Hello ", name))
# }
 
# x <- args[1]
# # x<-'Julia'
# sprintf("%s is best", x)
# 
# my_function(x)  
  
# print(args[1])



# ==================== GGSEG =============================================================================================================================
suppressMessages(library(ggseg3d))
suppressMessages(library(ggseg))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(htmltools))
suppressMessages(library(htmlwidgets))
suppressMessages(library(plotly))

options(warn=-1)
colors = c("indianred4",'indianred2','coral1','lightpink1','mistyrose1')  

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


# ==================== SAVING =============================================================================================================================
i = strtoi(args[1])
dk_data <- read.csv(file = paste('temp_folder/dk_R_Subtype ',i,'.csv', sep=''))
aseg_data <- read.csv(file = paste('temp_folder/aseg_R_Subtype ',i,'.csv', sep=''))


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

saveWidget(fig, paste('temp_folder/Subtype ',i,'.html', sep=''), selfcontained = F, libdir = paste('subtype ',i,' files', sep=''))

print(paste('PROGRESS: Subtype ',i,' files SUCCESFULLY saved in: /temp_folder', sep=''))



