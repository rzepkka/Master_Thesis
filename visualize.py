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
import os
from random import gauss
from scipy.stats import norm
import scipy.stats 

from mapping_2D import dk_dict, aseg_dict

def get_labels(S):
    unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
    subtype_labels = []
    for i in range(len(unique_subtypes)):
      subtype_labels.append('Subtype '+str(int(unique_subtypes[i])))        
    return subtype_labels

# ============= PATIENT-SPECIFIC PLOTS =============================================================================================================================================================

def subtype_probabilities(info, S, patient_id=0, subtype_labels = None, color=['#000000'],fontlist = [24, 18, 14, 22], width=900, height=600):
    """
    Creates a barplot for subtype probabilities
    :param info: pandas DataFrame with patients' data
    :param S: subtyping dictionary, subtypes for each patient individually
    :param patient_id: ID of a patient to visualize (optional)
    :param subtype_labels: list with label name for the subtypes (optional)
    :param colort: hex color value (optional)
    :param fontsize: a list of 4 ints, corresponding to [font_title, font_axes, font_ticks, font_legend] respectively (optional)
    :param width: int (optional)
    :param height: int (optional)
    :return: plotly express bar figure
    """  
    
    if patient_id not in list(info['PTID']) or patient_id is None:
        return 'Wrong patient ID', 'No prediction'
    else:
    
        # Get subtype labels
        unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
        if subtype_labels is None:
            subtype_labels = []
            for i in range(len(unique_subtypes)):
                subtype_labels.append('Subtype '+str(int(unique_subtypes[i])))

        subtype_map = {unique_subtypes[i]: subtype_labels[i] for i in range(len(subtype_labels))}

        subtypes = S['subtypes']
        subtypes = ['Outlier' if np.isnan(s) else subtype_map[s] for s in subtypes]    

        weights = S['subtypes_weights']

        # Create weight DataFrame
        df_weights = pd.DataFrame(weights, columns=subtype_labels)
        df_weights['Sum']=df_weights[subtype_labels].sum(axis = 1, skipna = True)
        df_weights['Prediction'] = subtypes
        
        # Count probabilities
        df_prob = pd.DataFrame()

        # TO CHANHE WHEN I GET DATA
        df_prob['Patient ID'] = info['PTID']
        for s in subtype_labels:
            df_prob[s]=round(df_weights[s]/df_weights['Sum']*100,2)
            
        data = df_prob[subtype_labels][df_prob['Patient ID']==patient_id]        
        prediction = np.array(df_weights['Prediction'][df_prob['Patient ID']==patient_id])[0]

        df = pd.DataFrame(data.values[0], data.columns)
        df = df.rename(columns={0: "Probability"})

        fig = px.bar(df, x=df.index, y="Probability",
                text=data.values[0],
                text_auto=True,
                width=width,
                height=height)

        # Styling 
        font_title, font_axes, font_ticks, font_bars = fontlist

        fig.update_layout(
            title_text='Subtype probabilities', # title of plot
            title_x=0.5,
            title_font_size=font_title,
            xaxis_title_text='Subtype', # xaxis label
            yaxis_title_text='Probability (%)', # yaxis label
            bargap=0.2, # gap between bars of adjacent location coordinates
        )

        fig.update_traces(marker_color=color,
                        textfont_size=font_bars,
                        texttemplate='%{text} %')

        fig.update_yaxes(title_font_size = font_axes, 
                        tickfont_size=font_ticks)

        fig.update_xaxes(title_font_size = font_axes, 
                        tickfont_size = font_ticks)


        return fig, prediction


def individual_staging(data, S, Sboot, diagnosis, patient_id,  color_list='#000000', width=950, height=400, fontsize=[24,20,20,22]):
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
                    
    fig = make_subplots(rows=2, cols = 1, 
                        subplot_titles=['Disease stage of the patient (representing mean and uncertainty of estimation)','Reference disease stage for different diagnostic groups'],
                       row_heights=[200,400])
 
    # ADD PATIENT-SPECIFIC
    boot=[]
    for b in range(len(Sboot)):
        boot.append(Sboot[b]['staging'][data['PTID']==patient_id][0])
        
        
        
    fig.add_trace(go.Box(x=boot, 
                         name='',
                         fillcolor=color_list[-1],
                         line_color='#000000',
                         opacity=0.8), row=1, col=1)

    for count, idx in enumerate(idx_list):
        fig.add_trace(go.Box(x=staging[idx], name=labels[count],
                             fillcolor=color_list[count],
                            line_color='#000000',
                            opacity=0.4),row=2,col=1)

    fig.update_xaxes(range=[0.0, 1.0])

    font_title, font_axes, font_ticks, font_legend = fontsize

    fig.update_layout(
            title_font_size=font_title,
            title_x=0.5,
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

    
    fig.update_xaxes(title_text="Disease Stage", tick0 = 0.0, dtick = 0.1,  row=2, col=1)
    
    for i in fig['layout']['annotations']:
        i['font'] = dict(size=20)

    return fig


def biomarker_distribution(data, T, subtype, patient_id=None):
    """
    Creates bimodal distribution plots for each biomarker, for chosen subtype; 
    vertical line corresponds to the patient-specific information
    :param data: pandas DataFrame with patients' data
    :param T: Snowphlake timeline object
    :param subtype: chosen disease subtype
    :param patient_id: ID of the chosen patient (optional)
    :return: plotly subplots with go Scatter figures
    """
    
    titles = [label.lower().replace("_"," ") for label in list(T.biomarker_labels)]
    num_rows = int(np.ceil(len(titles)/3))
    num_cols = 3
    fig = make_subplots(rows=num_rows, cols = num_cols, subplot_titles=titles)
    labels = T.biomarker_labels

    biomarker=0
    for row in range(1,num_rows+1):
        for col in range(1, num_cols+1):

            if biomarker >= len(labels):
                break
            
            if biomarker == (len(labels)-1):
                showlegend=True
            else: showlegend=False

            Dallis = data[labels[biomarker]]
            x_grid = np.linspace(np.min(data[labels[biomarker]]), np.max(data[labels[biomarker]]),1000)

            # CASES
            mu_cases = T.mixture_model.cases[biomarker]['mu'][0][subtype]
            sigma_cases = T.mixture_model.cases[biomarker]['std'][0][subtype]
            norm_pre = scipy.stats.norm(loc=mu_cases, scale=sigma_cases)
            likeli_pre=norm_pre.pdf(x_grid)

            likeli_pre=likeli_pre*(T.mixture_model.mixing[biomarker][subtype])

            # CONTROLS
            mu_controls = T.mixture_model.controls[biomarker]['mu'][0]
            sigma_controls = T.mixture_model.controls[biomarker]['std'][0]

            norm_post = scipy.stats.norm(loc=mu_controls, scale=sigma_controls)
            likeli_post=norm_post.pdf(x_grid)
            likeli_post=likeli_post*(1-T.mixture_model.mixing[biomarker][subtype])

            fig.add_trace(go.Scatter(x=x_grid, y=likeli_pre,
                                mode='lines',
                                name='Cases',
                                line=dict(color='red', width=2),
                                showlegend=showlegend,
                                    ), row=row, col=col)

            fig.add_trace(go.Scatter(x=x_grid, y=likeli_post,
                                mode='lines',
                                name='Controls',
                                line=dict(color='green', width=2),
                                showlegend=showlegend,
                                    ), row=row, col=col)


            if patient_id not in list(data['PTID']) or patient_id is None:
                continue
            else:
                patient = np.array(data[labels[biomarker]][data['PTID']==patient_id])[0]
                fig.add_vline(x=patient, line_width=2, line_dash="dash", line_color="red", row=row, col=col)            

            biomarker+=1
        
    # STYLING
    fig.update_layout(title = 'Biomarker Distribution',
                    title_x=0.5,
                    title_font_size=24,
                     height=1600,
                     width=1000,
                     legend_font_size=22
                     )
            
    return fig


# ============= PIE CHART =============================================================================================================================================================

def piechart_multiple(S, diagnosis, subtype_labels=None, chosen_subtypes = None):
    """
    Creates a pie chart of subtypes in the data   
    :param S: subtyping dictionary, subtypes for each patient individually
    :param diagnosis: np.array or list with diagnosis labels corresponding to records in S
    :param subtype_lables: a list with names of the subtype labels (optional)
    :param chosen_subtypes: a list with diagnosis labels to consider in the plot
    :return: plotly express Pie Figure
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

    # STYLING
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
    :param fontsize: a list of 4 ints, corresponding to [font_title, font_axes, font_ticks, font_legend] respectively (optional)
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

    # STYLING
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

    # STYLING
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



def staging_boxes(S, diagnosis, color_list=['#000000'], width=950, height=400, fontsize=[34,18,14,22]):
    """
    Creates a boxplot
    :param S: subtyping dictionary, subtypes for each patient individually
    :param diagnosis: np.array or list; with diagnosis labels corresponding to records in S
    :param color_list: a list with color hex values (optional)
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

    # STYLING
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


# ============= 2D PLOTTING =============================================================================================================================================================

def plot_dk_atlas(T,S, map_dk, subtype_labels = None, subtype = None, slider = None, save = False, filename='file'):     

    """
    Creates a dictionary, which can be used as input to ggseg.plot_dk() and plots it
    :param T: Timeline object
    :param S: subtyping dictionary, subtypes for each patient individually
    :param map_dk: dictionary with loaded JSON file for mapping cortical regions
    :param subtype_labels: a list with names of the subtypes (optional)
    :param subtype: name or index of the subtype to visualise (optional)  
    :param slider: int (optional)
    :param save: boolean, indicates whether plot should be displayed or saved in /temp_folder (optional)
    :param filename: string, how the saved filed should be names (if save=True, optional)
    :returns a figures by plt.show() from ggseg.plot_dk(), or saves the file using custom_dk() funtion
    """   
    
    if slider is None:
        dk = dk_dict(T, S, mapped_dict = map_dk, subtype = subtype,subtype_labels=subtype_labels)  
    else:
        dk_ = dk_dict(T, S, mapped_dict = map_dk, subtype = subtype,subtype_labels=subtype_labels)
        dk = {k: v for k, v in dk_.items() if v <= slider}
    
    if subtype is None:
        pass
    
    # save the images for animation
    if save is True:
                
        custom_dk(dk, cmap='Reds_r', figsize=(6,6),
                  vminmax = [0,1],
                  background='black', edgecolor='white', bordercolor='gray', title=f'{subtype}',fontsize = 24,
                 filename=filename)          
    
    else:
         ggseg.plot_dk(dk, cmap='Reds_r', figsize=(6,6),
              vminmax = [0,1],
              background='black', edgecolor='white', bordercolor='gray', 
                title=f'{subtype}',fontsize = 24)



def plot_aseg_atlas(T,S, map_aseg, subtype_labels = None, subtype = None, slider = None, save = False, filename='file'):     

    """
    Creates a dictionary, which can be used as input to ggseg.plot_aseg() function
    :param T: Timeline object
    :param S: subtyping dictionary, subtypes for each patient individually
    :param map_aseg: dictionary with loaded JSON file for mapping subcortical regions
    :param subtype_labels: a list with names of the subtypes (optional)
    :param subtype: name or index of the subtype to visualise (optional)  
    :param slider: int (optional)
    :param save: boolean, indicates whether plot should be displayed or saved in /temp_folder (optional)
    :param filename: string, how the saved filed should be names (if save=True, optional)
    :returns a figures by plt.show() from ggseg.plot_aseg(), or saves the file using custom_aseg() funtion
    """
    if slider is None:  
        aseg = aseg_dict(T,S, map_aseg,subtype = subtype,subtype_labels=subtype_labels)
    else:
        aseg_ = aseg_dict(T,S, map_aseg, subtype = subtype,subtype_labels=subtype_labels)
        aseg = {k: v for k, v in aseg_.items() if v <= slider}

    if subtype is None:
        pass 
    
    if save is True:
        
        custom_aseg(aseg, cmap='Reds_r', figsize=(6,6),
                  vminmax = [0,1],
                  background='black', edgecolor='white', bordercolor='gray', title=f'{subtype}',fontsize = 20,
                 filename=filename)      
        
    else:
        ggseg.plot_aseg(aseg, cmap='Reds_r', figsize=(6,6),
                vminmax = [0,1],
                background='black', edgecolor='white', bordercolor='gray', title=f'{subtype}', fontsize = 18)

def custom_dk(data, cmap='Spectral', background='k', edgecolor='w', ylabel='',
             figsize=(15, 15), bordercolor='w', vminmax=[], title='',
             fontsize=15, filename="file"):
    
    """
    plot_dk() function from ggseg library, customized so that output plots can be saven in /temp_folder and used for 2D animations
    """

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

    plt.savefig(f'temp_folder/snapshots/{filename}.png', bbox_inches='tight', pad_inches=0.20)
    # plt.savefig(f'video/{filename}.png', bbox_inches='tight', pad_inches=0.20)

    plt.close()
    

def custom_aseg(data, cmap='Spectral', background='k', edgecolor='w', ylabel='',
              figsize=(15, 5), bordercolor='w', vminmax=[],
              title='', fontsize=15, filename=""):

    """
    plot_aseg() function from ggseg library, customized so that output plots can be saven in /temp_folder and used for 2D animations
    """

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

    plt.savefig(f'temp_folder/snapshots/{filename}.png', bbox_inches='tight', pad_inches=0.1)
    # plt.savefig(f'video/{filename}.png', bbox_inches='tight', pad_inches=0.1)

    plt.close()

    # print('PROGRESS: Animations for subcortinal regions created.')



# ============= ADDITIONAL =============================================================================================================================================================

# Display patient's prediction on patient-specifin info page
def get_prediction(data, S, patient_id, subtype_labels=None):

    if patient_id not in list(data['PTID']) or patient_id is None:
            return 'Wrong patient ID', 'No prediction'
    else:
        unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
        if subtype_labels is None:
            subtype_labels = []
            for i in range(len(unique_subtypes)):
                subtype_labels.append('Subtype '+str(int(unique_subtypes[i])))

    subtype_map = {unique_subtypes[i]: subtype_labels[i] for i in range(len(subtype_labels))}
    subtypes = S['subtypes']
    subtypes = ['Outlier' if np.isnan(s) else subtype_map[s] for s in subtypes]    


    # subtypes
    patients = data['PTID']
    df = pd.DataFrame({'ID':patients, 'Prediction':subtypes})

    prediction = np.array(df['Prediction'][df['ID']==patient_id])[0]

    return prediction



# Create 2D gifs
def make_gif(frame_folder, subtype, atlas):
    file_list = glob.glob(f'{frame_folder}/*.png')
    file_list.sort()
    frames = [Image.open(image) for image in file_list]
    frame_one = frames[0]
    frame_one.save(f"temp_folder/2D_animations/{atlas}-{subtype}.gif", format="GIF", append_images=frames,
               save_all=True, duration=200, loop=0) 


# Embed PDF presentation in the App
def displayPDF(file):
    # Opening file from file path
    with open(file, "rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode('utf-8')

    # Embedding PDF in HTML
    pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="1300" height="1100" type="application/pdf"></iframe>'
    # Displaying File
    st.markdown(pdf_display, unsafe_allow_html=True)

