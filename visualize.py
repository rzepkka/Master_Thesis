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

from mapping_2D import dk_regions_2D, dk_dict_2D, aseg_dict_2D


def subtype_metrics(T):

    fig, ax = plt.subplots(3,1,figsize=(7, 7))
    ax[0].set_title('Metric: Residual sum of squares (RSS)')
    ax[0].plot(range(2,1+T.n_maxsubtypes),-np.diff(T.subtyping_model.rss_data),color='k')
    ax[0].plot(range(2,1+T.n_maxsubtypes),-np.diff(T.subtyping_model.rss_random))
    ax[0].set_xlabel('Number of subtypes')
    ax[0].set_ylabel('Change in RSS')
    ax[0].grid(visible=True,linestyle='--')
    ax[0].legend(['Data','Random'])
    m = np.max(-np.diff(T.subtyping_model.rss_random))
    ax[0].plot([2,T.n_maxsubtypes],[m,m],linestyle='dashed')

    ax[1].set_title('Silhouette score (SS)')
    ax[1].plot(range(1,1+T.n_maxsubtypes),T.subtyping_model.silhouette_score)
    ax[1].set_xlabel('Number of subtypes')
    ax[1].set_ylabel('SS')
    ax[1].grid(visible=True,linestyle='--')

    ax[2].set_title('Cophenetic correlation (CC)')
    ax[2].plot(range(1,1+T.n_maxsubtypes),T.subtyping_model.cophenetic_correlation)
    ax[2].set_xlabel('Number of subtypes')
    ax[2].set_ylabel('CC')
    ax[2].grid(visible=True,linestyle='--')
    fig.tight_layout()

    return fig, ax

def subtypes_piechart(S,diagnosis,diagnostic_labels_for_plotting,title = None,subtype_labels = None):

    fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
    
    unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
    if subtype_labels is None:
        subtype_labels = []
        for i in range(len(unique_subtypes)):
            subtype_labels.append('Subtype '+str(int(unique_subtypes[i])))
    
    n_data = np.zeros(len(unique_subtypes))
    for d in diagnostic_labels_for_plotting:
        idx_d = diagnosis == d
        for u in range(len(unique_subtypes)):
            n_data[u] = n_data[u] + np.sum(S['subtypes'][idx_d] == unique_subtypes[u])

    wedges, _ = ax.pie(n_data, wedgeprops=dict(width=0.5), startangle=-40)

    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"),
            bbox=bbox_props, zorder=0, va="center")

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(subtype_labels[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
                    horizontalalignment=horizontalalignment, **kw)
    
    if title is not None:
        ax.set_title(title)

    plt.show()

    return fig, ax
    

# ============= PIE CHART =============================================================================================================================================================

def subtype_piechart(S, subtype_labels=None):
    
    # Get subtype labels
    unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
    if subtype_labels is None:
        subtype_labels = []
        for i in range(len(unique_subtypes)):
            subtype_labels.append('Subtype '+str(int(unique_subtypes[i])))
    
    # Count number of each subtype
    counter = dict(Counter(S['subtypes'][~np.isnan(S['subtypes'])]))
    counter = collections.OrderedDict(sorted(counter.items()))
    
    # data for plotting
    values = list(counter.values())
    labels = list(subtype_labels)
    
    fig = go.Figure(data=[go.Pie(labels=labels, values=values, textinfo='value+percent')]) # , pull=[0, 0, 0.2, 0]
        
    fig.update_layout(title="Subtypes Pie Chart",
                      title_x=0.5,
                      title_font_size=24,
                      legend_font_size=18)
    
    return fig

# ============= EVENT CENTERS =============================================================================================================================================================

def event_centers(T, S, color_list=['#000000'], chosen_subtypes = None,
        subtype_labels = None, orderBy = None, width=1050, height=900, slider = None, fontsize=[34,18,14,22]):
    
    """
    Creates event centers box plots for multiple subtypes
    
    :param T: Timeline object
    :param S:
    :param color_list: a list with color names corresponding to each subtype, len(color_list) = len(subtypes). Preferably hex values
    :param chosen_subtypes: a list with names of the subtypes to visualize
    :param subtype_lables: a list with names of the subtype labels 
    :param orderBy: string, name of the subtype to order the boxplots by; default None
    :param width: chosen width of the returned plot
    :param height: chosen height of the returned plot
    :param slider: int, value of the slider from 2D visualizations
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
                      hovermode=False)

    fig.update_layout(legend_font_size=font_legend,
                    legend=dict(
                        yanchor="top",
                        y=0.90,
                        xanchor="right",
                        x=0.95))

    fig.add_vline(x=slider, line_width=2, line_dash="dash", line_color="red",
                  annotation_text=f"Slider value = {slider}",
                  annotation_position="top left",
                  annotation_font_color="red"
                  )

    return fig


# ============= PATIENT STAGING =============================================================================================================================================================

def patient_staging(S, diagnosis, color_list=['#000000'], num_bins=10, bin_width=0.02, width=1200, height=900, 
                        fontsize=[34,18,14,22], opacity=0.8):
    """
    Creates a barplot
    :param S: dictionary, Snowphlake output
    :param diagnosis: np.array or list; with diagnosis labels corresponding to records in S
    :param color_list: list with color hex values
    :param num_bins: int, how many bins should be displayed
    :param bin_width: int
    :return: plotly go Bar figure
    """  
    
    # Convert NaNs to 0.0
    staging = np.array([np.float64(0.0) if np.isnan(stage) else stage for stage in S['staging']])
   
    # Count number of each subtype occurences
    counter = dict(Counter(diagnosis))
        
    # Get labels
    labels = list(set(diagnosis))[::-1]
    
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
    num_bins = num_bins
    bin_width = np.repeat(bin_width, num_bins)
          
    color_list = color_list
        
    count=-1    
    num_bins = num_bins
    bar_width = np.repeat(0.02, num_bins)
    counter = dict(Counter(diagnosis))

    fig = go.Figure()

    for count, idx in enumerate(idx_list):
        freq,binc=np.histogram(staging[idx],bins=num_bins)
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
        barmode='group',
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
    :param S: dictionary, Snowphlake output
    :param diagnosis: np.array or list; with diagnosis labels corresponding to records in S
    :param color_list: list with color hex values
    :return: plotly go Box figure
    """
    
    # Convert NaNs to 0.0
    staging = np.array([np.float64(0.0) if np.isnan(stage) else stage for stage in S['staging']])
   
    # Count number of each subtype occurences
    counter = dict(Counter(diagnosis))
        
    # Get labels
    labels = list(set(diagnosis))[::-1]
    
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
            # legend=dict(
            #     yanchor="top",
            #     y=0.97,
            #     xanchor="right",
            #     x=0.97),
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
    :param S: dictionary, Snowphlake output
    :param diagnosis: np.array or list; with diagnosis labels corresponding to records in S
    :param color_list: list with color hex values
    :param num_bins: int, how many bins should be displayed
    :param bin_width: int
    :return: plotly go Bar figure
    """  
    
    # Convert NaNs to 0.0
    atypical = np.array([np.float64(0.0) if np.isnan(x) else x for x in S['atypicality']])
   
    # Count number of each subtype occurences
    counter = dict(Counter(diagnosis))
        
    # Get labels
    # labels = list(set(diagnosis))
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
        
    count=-1    
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
        legend=dict(
            yanchor="top",
            y=0.95,
            xanchor="right",
            x=0.95),
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
    :param S: dictionary, Snowphlake output
    :param diagnosis: np.array or list; with diagnosis labels corresponding to records in S
    :param color_list: list with color hex values
    :return: plotly go Box figure
    """
    
    # Convert NaNs to 0.0
    atypical = np.array([np.float64(0.0) if np.isnan(x) else x for x in S['atypicality']])
   
    # Count number of each subtype occurences
    counter = dict(Counter(diagnosis))
        
    # Get labels
    labels = list(set(diagnosis))
    
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
            legend=dict(
                yanchor="top",
                y=0.97,
                xanchor="right",
                x=0.97),
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


def staging_scatterplot(S, diagnosis, color_list = ['#000000']):
    """
    Creates a scatterplot of staging vs. atypicality
    :param S: dictionary, Snowphlake output
    :return: plotly scatterplot
    """     
    staging = list(S['staging'])
    atypical = list(S['atypicality'])
    diagnosis = list(diagnosis)
    
    df = pd.DataFrame(list(zip(staging, atypical,diagnosis)),
               columns =['Stage', 'Atypicality','Diagnosis'])
    df = df[df['Diagnosis'] != 'CN']
    
    fig = px.scatter(df, x='Stage', y='Atypicality', color='Diagnosis')

    return fig  

# ============= 2D PLOTTING =============================================================================================================================================================

def plot_dk_atlas(T,S, subtype_labels = None, subtype = None, slider = None):     
    """
    Creates a dictionary, which can be used as input to ggseg.plot_dk() function
    :param T: timeline object from snowphlake
    :param S: dictionary from snowphlake
    :param subtype_labels: a list with names of the subtypes (optional)
    :param subtype: name or index of the subtype to visualise (optional)  
    :param slider: int
    :returns a figures by plt.show() -> ggseg.plot_dk() 
    """
    mapped_dict = dk_regions_2D(T)    
    
    if slider is None:
        dk = dk_dict_2D(T, S, mapped_dict = mapped_dict, subtype = subtype)  
    else:
        dk_ = dk_dict_2D(T, S, mapped_dict = mapped_dict, subtype = subtype)
        dk = {k: v for k, v in dk_.items() if v <= slider}
           
    if subtype is None:
        # subtype = 'default = 0'
        pass
    else:
        return ggseg.plot_dk(dk, cmap='Reds_r', 
                        figsize=(6,6),
                        vminmax = [0,1],
                        background='black', 
                        edgecolor='white', 
                        bordercolor='gray', 
                        title=f'Cortical regions',
                        fontsize = 24)


def plot_aseg_atlas(T,S, subtype_labels = None, subtype = None, slider = None):     
    """
    Creates a dictionary, which can be used as input to ggseg.plot_aseg() function
    :param T: timeline object from snowphlake
    :param S: dictionary from snowphlake
    :param subtype_labels: a list with names of the subtypes (optional)
    :param subtype: name or index of the subtype to visualise (optional)  
    :param slider: int
    :returns a figures by plt.show() -> ggseg.plot_aseg()
    """
    if slider is None:  
        aseg = aseg_dict_2D(T,S, subtype = subtype)
    else:
        aseg_ = aseg_dict_2D(T,S, subtype = subtype)
        aseg = {k: v for k, v in aseg_.items() if v <= slider}

    if subtype is None:
        # subtype = 'Subtype 0'
        pass 
    else:
        return ggseg.plot_aseg(aseg, 
                            cmap='Reds_r', 
                            figsize=(6,6),
                            vminmax = [0,1],
                            background='black', 
                            edgecolor='white', 
                            bordercolor='gray', 
                            title=f'Subcortical regions', 
                            fontsize = 18)







