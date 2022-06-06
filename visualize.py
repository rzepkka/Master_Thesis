# Author: Vikram Venkatraghavan, Amsterdam UMC

from matplotlib import pyplot as plt 
import numpy as np 
import pandas as pd 
import plotly.express as px
import plotly.io as pio
import plotly.graph_objs as go
import plotly.offline as pyo
import ggseg
from collections import Counter
from plotly.subplots import make_subplots
import collections
import glob
from PIL import Image

from mapping_2D import dk_dict, aseg_dict
# from mapping_2D import dk_regions_2D, dk_dict_2D, aseg_dict_2D


# ============= PIE CHART =============================================================================================================================================================

def piechart_multiple(S, diagnosis, subtype_labels=None, chosen_subtypes = None):
    """
    Creates a pie chart of subtypes in the data   
    :param S: subtyping dictionary, subtypes for each patient individually
    :param diagnosis: np.array or list with diagnosis labels corresponding to records in S
    :param subtype_lables: a list with names of the subtype labels (optional)
    :param chosen_subtypes: a list with diagnosis labels to consider in the plot
    :return: plotly Pie Figure
    """
    
    # Get subtype labels
    unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
    if subtype_labels is None:
        subtype_labels = []
        for i in range(len(unique_subtypes)):
            subtype_labels.append('Subtype '+str(int(unique_subtypes[i])))
            
    # Set colors
    default_colors = ['#d2fafc','#62a2f5','#66cfff','#0e4287','#5894ed','#1217a3']
    color_map = {subtype_labels[i] : default_colors[i] for i in range(len(subtype_labels))}
    color_map['Outliers']='#8e9191'
            
    subtype_labels.append("Outliers")
    
            
    subtypes = list(S['subtypes'])
    subtypes = ['Outliers' if np.isnan(s) else s for s in subtypes]
    
    # Get diagnosis lables (exclude controls)
    diagnosis_labels = list(set(diagnosis))
    diagnosis_labels.remove('CN')
    diagnosis_labels.sort()
    
    if chosen_subtypes is None:
        chosen_subtypes=diagnosis_labels

    df = {'Diagnosis':diagnosis, 'Subtype':subtypes}
    df = pd.DataFrame(df)
    df = df[df['Diagnosis']!='CN']
        
    # Create DataFrame for quering to plotting
    df_subtypes = subtype_labels*len(diagnosis_labels)
    df_diagnosis = [d for d in diagnosis_labels for i in range(len(subtype_labels))]
    counts = list(df.groupby(['Diagnosis','Subtype']).size())
    df = pd.DataFrame({'Diagnosis':df_diagnosis, 'Subtype':df_subtypes, 'Counts': counts})

    # Query diagnosis for plotting
    df_plot = df[df['Diagnosis'].isin(chosen_subtypes)]

    fig = px.pie(df_plot, values='Counts', names='Subtype', title='Subtypes',color='Subtype',
                color_discrete_map=color_map)

    # style the plot
    fig.update_traces(textposition='inside', textinfo='value+percent')   
    
    fig.update_layout(title='',
                  title_x=0.5,
                  title_font_size=24,
                  legend_font_size=18)

    return fig


# ============= EVENT CENTERS =============================================================================================================================================================

def event_centers(T, S, color_list=['#000000'], chosen_subtypes = None,
        subtype_labels = None, orderBy = None, width=1050, height=900, slider = None, fontsize=[34,18,14,22]):
    
    """
    Creates event centers box plot for multiple subtypes  
    :param T: Timeline object
    :param S: subtyping dictionary, subtypes for each patient individually
    :param color_list: a list with color names corresponding to each subtype, len(color_list) = len(subtypes); preferably hex values
    :param chosen_subtypes: a list with names of the subtypes to visualize
    :param subtype_lables: a list with names of the subtype labels (optional)
    :param orderBy: string, name of the subtype to order the boxplots by (optional)
    :param width: chosen width of the returned plot (optional)
    :param height: chosen height of the returned plot (optional)
    :param slider: int, value of the slider from 2D visualizations (optional)
    :return: plotly box figure
    """

    unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
    if subtype_labels is None:
        subtype_labels = []
        for i in range(len(unique_subtypes)):
            subtype_labels.append('Subtype '+str(int(unique_subtypes[i])))
                    
    if orderBy is None:
        orderBy = subtype_labels[0]
                
    if chosen_subtypes is None:
        chosen_subtypes = subtype_labels
        
    num_subtypes = len(subtype_labels)
    
    labels = T.biomarker_labels
    labels_cleaned = map(lambda x: x.replace("-"," "), labels)
    labels_cleaned = map(lambda x: x.replace("_"," "), labels_cleaned)
    labels_cleaned = list(map(lambda x: x.lower(), labels_cleaned))
    
    # key: value --> ordering: region_name
    labels_dict = {num: label for num, label in enumerate(labels_cleaned)}
    
    color_map = {chosen_subtypes[i]: color_list[i] for i in range(len(color_list))}

    # EVENT-CENTERS
    evn = []
    reg = []
    subs = []

    for b in range(T.bootstrap_repetitions):
        for i, s in enumerate(subtype_labels):
            for r in range(len(labels)):
                
                # SUBTYPES 
                subs.append(s)
                
                # EVENT-CENTERS
                evn.append(T.bootstrap_sequence_model[b]['event_centers'][i][r])
                
                # CORRESPONDING REGIONS
                label_number = T.bootstrap_sequence_model[b]['ordering'][i][r]
                reg.append(labels_dict[label_number])
                
                    
    dic = {'Region':reg, 'Subtype':subs, 'Score':evn}
    df = pd.DataFrame(dic)
        
    fig = px.box(df[df['Subtype'].isin(chosen_subtypes)], 
                 x="Score", 
                 y="Region", 
                 color = 'Subtype',
                color_discrete_map=color_map,
                 title=f"Event Centers", width=width, height=height, 
                 labels={"Score": "Disease Stage",  "Region": "Region Names"})
    
    df_sortBy = df[df['Subtype']==orderBy].drop(columns=['Subtype'])

    # GROUP BY MEDIAN
    df_sorted = df_sortBy.groupby('Region').quantile(q=0.5).sort_values(by='Score', ascending = True)

    # GROUP BY MEAN
    # df_sorted = df_sortBy.groupby('Region').aggregate('mean').sort_values(by='Score', ascending = True)


    labels_sorted = list(df_sorted.index)
    labels_sorted.reverse()

    font_title, font_axes, font_ticks, font_legend = fontsize

    fig.update_yaxes(categoryarray=labels_sorted, 
                    categoryorder="array", 
                    title_font_size = font_axes, 
                    tickfont_size=font_ticks)

    fig.update_xaxes(title_font_size = font_axes, 
                    tickfont_size = font_ticks)
    
    fig.update_layout(xaxis = dict(tickmode = 'linear', 
                                   tick0 = 0.0, 
                                   dtick = 0.1),
                      title_font_size=font_title,
                      title_x=0.5,
                      # hovermode=False,
                      legend_font_size=font_legend)

    fig.add_vline(x=slider, line_width=2, line_dash="dash", line_color="red",
                  annotation_text=f"Slider value = {slider}",
                  annotation_position="top left",
                  annotation_font_color="red"
                  )

    return fig


# ============= PATIENT STAGING =============================================================================================================================================================

def patient_staging(S, diagnosis, color_list=['#000000'], num_bins=10, bin_width=0.02, width=1200, height=900, 
                        fontsize=[34,18,14,22], opacity=0.8, barmode = 'group'):
    """
    Creates a barplot for diagnosis
    :param S: subtyping dictionary, subtypes for each patient individually
    :param diagnosis: np.array or list with diagnosis labels corresponding to records in S
    :param color_list: list with color hex values
    :param num_bins: int, how many bins should be displayed (optional)
    :param bin_width: int (optional)
    :param width: int (optional)
    :param height: int (optional)
    :param fontsize: a list of 4 ints, corresponding to [font_title, font_axes, font_ticks, font_legend] respectively (optional)
    :param opacity: float (optional)
    :param barmode: string, 'group' or 'stack'
    :return: plotly go Bar figure
    """  
    
    # Convert NaNs to 0.0
    staging = np.array([np.float64(0.0) if np.isnan(stage) else stage for stage in S['staging']])
   
    # Count number of each subtype occurences
    counter = dict(Counter(diagnosis))
        
    # Get labels
    labels = list(set(diagnosis[diagnosis!='CN']))
    labels.sort()
    
    # Get indexes
    diagnosis = np.array(diagnosis)
    staging = np.array(staging)
    
    # Get indexes for each diagnostic label
    idx_list = []
    for l in labels:
        if l!='CN':
            idx = np.where(diagnosis==l)
            idx = idx[0]
            idx_list.append(idx)             
      
    # Bar settings
    freq,binc=np.histogram(staging,bins=np.linspace(0,1,num_bins+1)) # Center the bins
    bin_width = np.repeat(bin_width, num_bins)
                  
    bar_width = np.repeat(bin_width, num_bins)
    counter = dict(Counter(diagnosis))

    fig = go.Figure()

    for count, idx in enumerate(idx_list):
        freq,binc=np.histogram(staging[idx],bins=binc) # , range = (0.,1.)
        freq = (1.*freq)/len(staging)
        label = labels[count]

        fig.add_trace(go.Bar(
                    x=binc,
                    y=freq,
                    name=f'{label} (n = {counter[label]})', 
                    width=bin_width,
                    marker_color=color_list[count],
                    opacity=opacity
        )) 
                
    font_title, font_axes, font_ticks, font_legend = fontsize

    fig.update_layout(
        title="Patient Staging",
        title_font_size=font_title,
        title_x=0.5,
        xaxis_title="Disease Stage",
        yaxis_title="Frequency of occurences",
        xaxis = dict(
            tickmode = 'linear',
            tick0 = 0.0,
            dtick = 0.1
        ),
        barmode=barmode,
        legend_font_size=font_legend,
        legend=dict(
            yanchor="top",
            y=0.95,
            xanchor="right",
            x=0.95),
        autosize = False,
        width=width,
        height=height
    )
    
    fig.update_xaxes(range=[-0.05, 1.0])
    
    fig.update_yaxes(title_font_size = font_axes, 
                    tickfont_size=font_ticks)
    
    fig.update_xaxes(title_font_size = font_axes, 
                    tickfont_size = font_ticks)

    return fig



def staging_boxes(S, diagnosis, color_list='#000000', width=950, height=400, fontsize=[34,18,14,22]):
    """
    Creates a boxplot
    :param S: subtyping dictionary, subtypes for each patient individually
    :param diagnosis: np.array or list; with diagnosis labels corresponding to records in S
    :param color_list: list with color hex values (optional)
    :param width: int (optional)
    :param height: int (optional)
    :param fontsize: a list of 4 ints, corresponding to [font_title, font_axes, font_ticks, font_legend] respectively (optional)
    :return: plotly go Box figure
    """
    
    # Convert NaNs to 0.0
    staging = np.array([np.float64(0.0) if np.isnan(stage) else stage for stage in S['staging']])
   
    # Count number of each subtype occurences
    counter = dict(Counter(diagnosis))
        
    # Get labels
    labels = list(set(diagnosis[diagnosis!='CN']))
    labels.sort()
    
    # Get indexes
    diagnosis = np.array(diagnosis)
    staging = np.array(staging)
    
    # Get indexes for each diagnostic label
    idx_list = []
    for l in labels:
        if l!='CN':
            idx = np.where(diagnosis==l)
            idx = idx[0]
            idx_list.append(idx)
           
    fig = go.Figure()

    for count, idx in enumerate(idx_list):
        fig.add_trace(go.Box(x=staging[idx], name=labels[count],
                             fillcolor=color_list[count],
                            line_color='#000000',
                            opacity=0.8))

    fig.update_xaxes(range=[-0.05, 1.0])

    font_title, font_axes, font_ticks, font_legend = fontsize

    fig.update_layout(
            # title="Staging - Boxplots",
            title_font_size=font_title,
            title_x=0.5,
            xaxis_title="Disease Stage",
            yaxis_title="Diagnosis",
            xaxis = dict(
                tickmode = 'linear',
                tick0 = 0.0,
                dtick = 0.1
            ),
            legend_font_size=font_legend,
            showlegend=False,
            autosize = False,
            width=width,
            height=height
        )
    
    fig.update_yaxes(title_font_size = font_axes, 
                    tickfont_size=font_ticks)
    
    fig.update_xaxes(title_font_size = font_axes, 
                    tickfont_size = font_ticks)

    return fig



# ============= ATYPICALITY =============================================================================================================================================================

def atypicality(S, diagnosis, color_list=['#000000'], num_bins=10, bin_width=0.02, width=1200, height=900, 
                        fontsize=[34,18,14,22], opacity=0.8):
    """
    Creates a barplot
    :param S: subtyping dictionary, subtypes for each patient individually
    :param diagnosis: np.array or list; with diagnosis labels corresponding to records in S
    :param color_list: list with color hex values (optional)
    :param num_bins: int, how many bins should be displayed
    :param bin_width: int (optional)
    :param width: int (optional)
    :param height: int (optional)
    :param fontsize: a list of 4 ints, corresponding to [font_title, font_axes, font_ticks, font_legend] respectively (optional)
    :param opacity: float (optional)

    :return: plotly go Bar figure
    """  
    
    # Convert NaNs to 0.0
    atypical = np.array([np.float64(0.0) if np.isnan(x) else x for x in S['atypicality']])
   
    # Count number of each subtype occurences
    counter = dict(Counter(diagnosis))
        
    # Get labels
    labels = list(set(diagnosis[diagnosis!='CN']))
    labels.sort()

    # Get indexes
    diagnosis = np.array(diagnosis)
    atypical = np.array(atypical)
    
    # Get indexes for each diagnostic label
    idx_list = []
    for l in labels:
        if l!='CN':
            idx = np.where(diagnosis==l)
            idx = idx[0]
            idx_list.append(idx)

    # Bar settings
    num_bins = num_bins
    bin_width = np.repeat(bin_width, num_bins)
          
    color_list = color_list
        
    num_bins = num_bins
    bar_width = np.repeat(0.02, num_bins)
    # counter = dict(Counter(diagnosis))
    counter = dict(Counter(diagnosis[diagnosis!='CN']))

    fig = go.Figure()
    
    for count, idx in enumerate(idx_list):
        freq,binc=np.histogram(atypical[idx],bins=num_bins)
        freq = (1.*freq)/len(atypical)
        
        label = labels[count]

        fig.add_trace(go.Bar(
                    x=binc[:-1],
                    y=freq,
                    name=f'{label} (n = {counter[label]})',
                    width=bin_width,
                    marker_color=color_list[count],
                    opacity=opacity
        )) 

    font_title, font_axes, font_ticks, font_legend = fontsize
                
    fig.update_layout(
        title="Atypicality",
        title_font_size=font_title,
        title_x=0.5,
        xaxis_title="Value",
        yaxis_title="Frequency of occurences",
        xaxis = dict(
            tickmode = 'linear',
            tick0 = 0.0,
            dtick = 2
        ),
        barmode='group',
        legend_font_size=font_legend,
        autosize = False,
        width=width,
        height=height
    )
    
    fig.update_xaxes(range=[np.min(atypical)-1.5, np.max(atypical)])
    
    fig.update_yaxes(title_font_size = font_axes, 
                    tickfont_size=font_ticks)
    
    fig.update_xaxes(title_font_size = font_axes, 
                    tickfont_size = font_ticks)

    return fig


def atypicality_boxes(S, diagnosis, color_list='#000000', width=950, height=400, fontsize=[34,18,14,22]):
    """
    Creates a boxplot
    :param S: subtyping dictionary, subtypes for each patient individually
    :param diagnosis: np.array or list; with diagnosis labels corresponding to records in S
    :param color_list: list with color hex values
    :param width: int (optional)
    :param height: int (optional)
    :param fontsize: a list of 4 ints, corresponding to [font_title, font_axes, font_ticks, font_legend] respectively (optional)
    :return: plotly go Box figure
    """
    
    # Convert NaNs to 0.0
    atypical = np.array([np.float64(0.0) if np.isnan(x) else x for x in S['atypicality']])
   
    # Count number of each subtype occurences
    counter = dict(Counter(diagnosis))
        
    # Get labels
    labels = list(set(diagnosis[diagnosis!='CN']))
    labels.sort()
    
    # Get indexes
    diagnosis = np.array(diagnosis)
    atypical = np.array(atypical)
    
    # Get indexes for each diagnostic label
    idx_list = []
    for l in labels:
        if l!='CN':
            idx = np.where(diagnosis==l)
            idx = idx[0]
            idx_list.append(idx)
        
    
    fig = go.Figure()

    for count, idx in enumerate(idx_list):
        fig.add_trace(go.Box(x=atypical[idx], name=labels[count],
                             fillcolor=color_list[count],
                            line_color='#000000'))

    font_title, font_axes, font_ticks, font_legend = fontsize
    
    fig.update_layout(
            # title="Atypicality - Boxplots",
            title_font_size=font_title,
            title_x=0.5,
            xaxis_title="Value",
            yaxis_title="Diagnosis",
            xaxis = dict(
                tickmode = 'linear',
                tick0 = 0.0,
                dtick = 5
            ),
            legend_font_size=font_legend,
            showlegend=False,
            autosize = False,
            width=width,
            height=height
        )

    fig.update_xaxes(range=[np.min(atypical)-1.5, np.max(atypical)])

    
    fig.update_yaxes(title_font_size = font_axes, 
                    tickfont_size=font_ticks)
    
    fig.update_xaxes(title_font_size = font_axes, 
                    tickfont_size = font_ticks)

    return fig


# ============= SCATTEEPLOT =============================================================================================================================================================

def staging_scatterplot(S, diagnosis, subtype_labels = None, chosen_subtypes = None, color_list = ['#000000'], width=1100, height=800, fontsize=[34,18,14,22]):
    """
    Creates a scatterplot of staging ~ atypicality
    :param S: subtyping dictionary, subtypes for each patient individually
    :param subtype_labels: list with label name for the subtypes (optional)
    :param chosen_subtypes: a list with subtype labels to consider in the plot
    :param color_list: list with color hex values (optional)
    :param width: int (optional)
    :param height: int (optional)
    :param fontsize: a list of 4 ints, corresponding to [font_title, font_axes, font_ticks, font_legend] respectively (optional)
    :return: plotly scatterplot
    """     
    
    # Get subtype labels
    unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
    if subtype_labels is None:
        subtype_labels = []
        for i in range(len(unique_subtypes)):
            subtype_labels.append('Subtype '+str(int(unique_subtypes[i])))
            
    subtype_map = {unique_subtypes[i]: subtype_labels[i] for i in range(len(subtype_labels))}
    
    # Get diagnosis lables (exclude controls)
    diagnosis_labels = list(set(diagnosis))
    diagnosis_labels.remove('CN')
    diagnosis_labels.sort()

    if chosen_subtypes is None:
        chosen_subtypes=diagnosis_labels
       
    # Create DataFrame
    staging = list(S['staging'])
    atypical = list(S['atypicality'])
    diagnosis = list(diagnosis)
    subtype = list(S['subtypes'])
    
    df = pd.DataFrame(list(zip(staging, atypical,subtype, diagnosis)),
               columns =['Stage', 'Atypicality','Subtype','Diagnosis'])
    df = df[df['Diagnosis'] != 'CN']
    df['Subtype'] = df['Subtype'].apply(lambda x: x if np.isnan(x) else subtype_map[x])
    df = df.dropna(axis=0, subset=['Subtype'])
    
    color_map = {diagnosis_labels[i]: color_list[i] for i in range(len(diagnosis_labels))}

    font_title, font_axes, font_ticks, font_legend = fontsize    

    df_plot = df[df['Subtype'].isin(chosen_subtypes)] # df['Diagnosis']
    
    fig = px.scatter(df_plot, x='Stage', y='Atypicality', color='Diagnosis', color_discrete_map=color_map) # color='Subtype'

    # style the plot
    fig.update_layout(
        title="Staging ~ Atypicality",
        title_font_size=font_title,
        title_x=0.5,
        xaxis_title="Stage",
        yaxis_title="Atypicality",
        xaxis = dict(
            tickmode = 'linear',
            tick0 = 0.0,
            dtick = 2
        ),
        barmode='group',
        legend_font_size=font_legend,
        autosize = False,
        width=width,
        height=height
    )
    
    fig.update_xaxes(range=[np.min(atypical)-1.5, 
                            np.max(atypical)],
                     title_font_size = font_axes, 
                    tickfont_size = font_ticks
                    )
    
    fig.update_yaxes(title_font_size = font_axes, 
                    tickfont_size=font_ticks)

    return fig  


# ============= 2D PLOTTING =============================================================================================================================================================

# def plot_dk_atlas(T,S, map_dk, subtype_labels = None, subtype = None, slider = None):     

#     """
#     Creates a dictionary, which can be used as input to ggseg.plot_dk() and plots it
#     :param T: Timeline object
#     :param S: subtyping dictionary, subtypes for each patient individually
#     :param subtype_labels: a list with names of the subtypes (optional)
#     :param subtype: name or index of the subtype to visualise (optional)  
#     :param slider: int (optional)
#     :returns a figures by plt.show() -> ggseg.plot_dk() 
#     """   
    
#     if slider is None:
#         dk = dk_dict(T, S, mapped_dict = map_dk, subtype = subtype)  
#     else:
#         dk_ = dk_dict(T, S, mapped_dict = map_dk, subtype = subtype)
#         dk = {k: v for k, v in dk_.items() if v <= slider}
        
    
#     if subtype is None:
#         pass
#     else:
#         return ggseg.plot_dk(dk, cmap='Reds_r', figsize=(6,6),
#                   vminmax = [0,1],
#                   background='black', edgecolor='white', bordercolor='gray', title=f'{subtype}',fontsize = 24)

def plot_dk_atlas(T,S, map_dk, subtype_labels = None, subtype = None, slider = None, save = False, filename='file'):     

    """
    Creates a dictionary, which can be used as input to ggseg.plot_dk() and plots it
    :param T: Timeline object
    :param S: subtyping dictionary, subtypes for each patient individually
    :param subtype_labels: a list with names of the subtypes (optional)
    :param subtype: name or index of the subtype to visualise (optional)  
    :param slider: int (optional)
    :returns a figures by plt.show() -> ggseg.plot_dk() 
    """   
    
    if slider is None:
        dk = dk_dict(T, S, mapped_dict = map_dk, subtype = subtype)  
    else:
        dk_ = dk_dict(T, S, mapped_dict = map_dk, subtype = subtype)
        dk = {k: v for k, v in dk_.items() if v <= slider}
    

    
    if subtype is None:
        # subtype = 'default = 0'
        pass
    
    # save the images for animation
    elif save is True:
                
        custom_dk(dk, cmap='Reds_r', figsize=(6,6),
                  vminmax = [0,1],
                  background='black', edgecolor='white', bordercolor='gray', title=f'Subtype 0',fontsize = 24,
                 filename=filename)          
    
    else:
         return ggseg.plot_dk(dk, cmap='Reds_r', figsize=(6,6),
              vminmax = [0,1],
              background='black', edgecolor='white', bordercolor='gray', 
                title=f'{subtype}',fontsize = 24)


# def plot_aseg_atlas(T,S, map_aseg, subtype_labels = None, subtype = None, slider = None):     

#     """
#     Creates a dictionary, which can be used as input to ggseg.plot_aseg() function
#     :param T: Timeline object
#     :param S: subtyping dictionary, subtypes for each patient individually
#     :param subtype_labels: a list with names of the subtypes (optional)
#     :param subtype: name or index of the subtype to visualise (optional)  
#     :param slider: int (optional)
#     :returns a figures by plt.show() -> ggseg.plot_aseg()
#     """
#     if slider is None:  
#         aseg = aseg_dict(T,S, map_aseg,subtype = subtype)
#     else:
#         aseg_ = aseg_dict(T,S, map_aseg, subtype = subtype)
#         aseg = {k: v for k, v in aseg_.items() if v <= slider}

#     if subtype is None:
#         # subtype = 'Subtype 0'
#         pass 
#     else:
#         return ggseg.plot_aseg(aseg, cmap='Reds_r', figsize=(6,6),
#                 vminmax = [0,1],
#                 background='black', edgecolor='white', bordercolor='gray', title=f'{subtype}', fontsize = 18)

def plot_aseg_atlas(T,S, map_aseg, subtype_labels = None, subtype = None, slider = None, save = False, filename='file'):     

    """
    Creates a dictionary, which can be used as input to ggseg.plot_aseg() function
    :param T: Timeline object
    :param S: subtyping dictionary, subtypes for each patient individually
    :param subtype_labels: a list with names of the subtypes (optional)
    :param subtype: name or index of the subtype to visualise (optional)  
    :param slider: int (optional)
    :returns a figures by plt.show() -> ggseg.plot_aseg()
    """
    if slider is None:  
        aseg = aseg_dict(T,S, map_aseg,subtype = subtype)
    else:
        aseg_ = aseg_dict(T,S, map_aseg, subtype = subtype)
        aseg = {k: v for k, v in aseg_.items() if v <= slider}

    if subtype is None:
        # subtype = 'Subtype 0'
        pass 
    
    elif save is True:
        
        custom_aseg(aseg, cmap='Reds_r', figsize=(6,6),
                  vminmax = [0,1],
                  background='black', edgecolor='white', bordercolor='gray', title=f'Subtype 0',fontsize = 20,
                 filename=filename)      
        
    else:
        ggseg.plot_aseg(aseg, cmap='Reds_r', figsize=(6,6),
                vminmax = [0,1],
                background='black', edgecolor='white', bordercolor='gray', title=f'{subtype}', fontsize = 18)

# Customized function to save plt

def custom_dk(data, cmap='Spectral', background='k', edgecolor='w', ylabel='',
             figsize=(15, 15), bordercolor='w', vminmax=[], title='',
             fontsize=15, filename="file"):
    
    import ggseg
    import matplotlib.pyplot as plt
    import os.path as op
    from glob import glob
    from ggseg import _create_figure_, _render_regions_, _get_cmap_, _render_data_, _add_colorbar_

    wd = op.join(op.dirname(ggseg.__file__), 'data', 'dk')

    # A figure is created by the joint dimensions of the whole-brain outlines
    whole_reg = ['lateral_left', 'medial_left', 'lateral_right',
                 'medial_right']
    files = [open(op.join(wd, e)).read() for e in whole_reg]
    ax = _create_figure_(files, figsize, background, title, fontsize, edgecolor)

    # Each region is outlined
    reg = glob(op.join(wd, '*'))
    files = [open(e).read() for e in reg]
    _render_regions_(files, ax, bordercolor, edgecolor)

    # For every region with a provided value, we draw a patch with the color
    # matching the normalized scale
    cmap, norm = _get_cmap_(cmap, data.values(), vminmax=vminmax)
    _render_data_(data, wd, cmap, norm, ax, edgecolor)

    # DKT regions with no provided values are rendered in gray
    data_regions = list(data.keys())
    dkt_regions = [op.splitext(op.basename(e))[0] for e in reg]
    NA = set(dkt_regions).difference(data_regions).difference(whole_reg)
    files = [open(op.join(wd, e)).read() for e in NA]
    _render_regions_(files, ax, 'gray', edgecolor)

    # A colorbar is added
    _add_colorbar_(ax, cmap, norm, edgecolor, fontsize*0.75, ylabel)

    plt.savefig(f'video/{filename}.png', bbox_inches='tight', pad_inches=0.2)
    plt.close()

    print('PROGRESS: Animations for cortinal regions created.')


def custom_aseg(data, cmap='Spectral', background='k', edgecolor='w', ylabel='',
              figsize=(15, 5), bordercolor='w', vminmax=[],
              title='', fontsize=15, filename=""):
    import matplotlib.pyplot as plt
    import os.path as op
    from glob import glob
    import ggseg
    from ggseg import _create_figure_, _render_regions_, _get_cmap_, _render_data_, _add_colorbar_

    wd = op.join(op.dirname(ggseg.__file__), 'data', 'aseg')
    reg = [op.basename(e) for e in glob(op.join(wd, '*'))]

    # Select data from known regions (prevents colorbar from going wild)
    known_values = []
    for k, v in data.items():
        if k in reg:
            known_values.append(v)

    whole_reg = ['Coronal', 'Sagittal']
    files = [open(op.join(wd, e)).read() for e in whole_reg]

    # A figure is created by the joint dimensions of the whole-brain outlines
    ax = _create_figure_(files, figsize, background,  title, fontsize, edgecolor)

    # Each region is outlined
    reg = glob(op.join(wd, '*'))
    files = [open(e).read() for e in reg]
    _render_regions_(files, ax, bordercolor, edgecolor)

    # For every region with a provided value, we draw a patch with the color
    # matching the normalized scale
    cmap, norm = _get_cmap_(cmap, known_values, vminmax=vminmax)
    _render_data_(data, wd, cmap, norm, ax, edgecolor)

    # The following regions are ignored/displayed in gray
    NA = ['Cerebellum-Cortex', 'Cerebellum-White-Matter', 'Brain-Stem']
    files = [open(op.join(wd, e)).read() for e in NA]
    _render_regions_(files, ax, '#111111', edgecolor)

    # A colorbar is added
    _add_colorbar_(ax, cmap, norm, edgecolor, fontsize*0.75, ylabel)

    plt.savefig(f'video/{filename}.png', bbox_inches='tight', pad_inches=0.1)
    plt.close()

    print('PROGRESS: Animations for subcortinal regions created.')




