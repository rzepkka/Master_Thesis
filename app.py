import streamlit as st
import streamlit.components as stc
import streamlit.components.v1 as components

import scipy.stats 
import pandas as pd
import numpy as np
import pickle
import numpy as np
import os
import webbrowser
from collections import Counter
import collections

import matplotlib.pyplot as plt
from matplotlib import rc

#
import plotly.graph_objs as go
from ipywidgets import Output

from visualize import event_centers, plot_dk_atlas, plot_aseg_atlas, patient_staging, staging_boxes, subtype_piechart

from visualize import atypicality, atypicality_boxes, staging_scatterplot, piechart_multiple

st.set_option('deprecation.showPyplotGlobalUse', False)


# LOAD PICKLE FILE
# read_input_file = open('data/ADC_FTLD_subtypes_agecorrected_zscore_final.pickle','rb')
read_input_file = open('data/EDADS_subtype_timelines_agecorrected_opt.pickle','rb')
load_inputs = pickle.load(read_input_file)
read_input_file.close()

T, S, Sboot = load_inputs

diagnosis = np.load('data/diagnosis.npy', allow_pickle=True)

# Get labels for options in select boxes
def get_labels(S):
    unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
    subtype_labels = []
    for i in range(len(unique_subtypes)):
        subtype_labels.append('Subtype '+str(int(unique_subtypes[i])))        
    return subtype_labels

labels = get_labels(S=S)

def main():

    st.set_page_config(layout="wide")
    st.sidebar.title("Adjusting Plots")

    # Connect .css file for styling the app components
    def local_css(file_name):
        with open(file_name) as f:
            st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
    local_css("style.css")

    st.header('Disease progression timeline')

    # SELECT PLOT
    plot_type_list = ['Disease timeline']
    chosen_plot_type = st.sidebar.selectbox("Select Plot", plot_type_list)

    # CHOOSE WIDTH AND HEIGHT
    chosen_width = st.sidebar.number_input('Select width of the plot in px',value=1000, max_value=1130)
    chosen_height = st.sidebar.number_input('Select height of the plot in px',value=800)

    if chosen_plot_type == 'Disease timeline':

        col_piechart, col_piechart_select = st.columns([5.2,3])

        with col_piechart_select:
            diagnosis_labels = list(set(diagnosis))
            diagnosis_labels.remove('CN')
            choose_pieplot = st.multiselect('Select diagnosis', diagnosis_labels, default=diagnosis_labels)


        with col_piechart:
            plot_piechart = piechart_multiple(S=S,
                                            diagnosis=diagnosis,
                                            chosen_subtypes=choose_pieplot)
            st.plotly_chart(plot_piechart)

        # SELECTION
        # options = ['Subtype 0','Subtype 1', 'Subtype 2', 'Subtype 3', 'Subtype 4']
        num_subtypes = len(labels)
 
        subtype_list = []

        col_select, col_select_blank = st.columns([5.2,3])
        with col_select:
            subtype_visualize = st.selectbox('Select a subtype to visualize:',labels)       
        subtype_list.append(subtype_visualize)

        # BUTTONS
        html_file = subtype_visualize.replace(" ","_")

        if st.sidebar.button('Open 3D visualization in a new tab'):
            filename = "file://"+os.getcwd()+ f"/html_3D/{html_file}.html"
            webbrowser.open(filename, new = 2)

        if st.sidebar.button('Download 3D visualisation as HTML file'):
            st.sidebar.warning('Feature not added yet')

        color_list = []
        default_color_list = ['#0000ff', '#880000', '#ffa07a', '#04977d', '#fd8ef3']  

        # ======================= 2D ===============================================================================================================
        
        col_cortical, col_subcortical, col_blank = st.columns([2,3.2,3])

        if subtype_visualize != None:

            col_slider, col_slider_blank = st.columns([5.2,3])
            with col_slider:
                slider = st.slider(label = 'Choose regions to display', 
                                    help = 'Only display regions with ordering <= to the chosen value',
                                    min_value = 0.0, 
                                    max_value = 1.0, 
                                    value=1.0, step=0.01)
        
            with col_cortical:
                ggseg_dk = plot_dk_atlas(T = T, 
                                        S = S, 
                                        subtype = subtype_visualize, 
                                        slider = slider)
                st.pyplot(ggseg_dk)

            with col_subcortical:     
                ggseg_aseg = plot_aseg_atlas(T = T, 
                                        S = S, 
                                        subtype = subtype_visualize, 
                                        slider = slider)       
                st.pyplot(ggseg_aseg)

        # ======================= EVENT CENTERS ===============================================================================================================

        # list for additional subtypes to compare
        options_compare = []
        for label in labels:
            if label == subtype_visualize:
                pass
            else: 
                options_compare.append(label)

        col_event_centers, col_event_centers_options = st.columns([3,1])

        with col_event_centers_options:
            st.subheader('Style event centers graph')
            chosen_subtypes = st.multiselect('Select additional subtypes to compare:', options_compare)

            title_font = st.number_input('Title size:',value=34)
            title_axes = st.number_input('Axis labels size:',value=18)
            title_ticks = st.number_input('Ticks size:',value=14)
            title_legend = st.number_input('Legend size:',value=22)

            font_list = [title_font, title_axes, title_ticks, title_legend]

            for subtype in chosen_subtypes:
                subtype_list.append(subtype)

            for idx, subtype in enumerate(subtype_list):
                subtype_color = st.text_input(f'Select color for {subtype}', value = f'{default_color_list[idx]}',placeholder='e.g. #000000')
                color_list.append(subtype_color)

        eventCenters = event_centers(T = T,
                                    S = S, 
                                    color_list = color_list,
                                    chosen_subtypes = subtype_list,
                                    subtype_labels = labels, 
                                    orderBy = subtype_visualize,
                                    width = chosen_width,
                                    height = chosen_height,
                                    slider = slider,
                                    fontsize = font_list)

    
        with col_event_centers:
            st.plotly_chart(eventCenters)

        # ADD DIVIDER
        st.markdown('---')

        # ======================= PATIENT STAGING ===============================================================================================================

        col_staging, col_staging_options = st.columns([3,1])

        diagnosis_labels = list(set(diagnosis))
        diagnosis_labels.sort()

        color_list = ['#00ff00','#377eb8','#4daf4a','#e41a1c']
        color_diagnosis =[]

        with col_staging:
            chosen_plot = st.selectbox('Select plot:',['Patient Staging', 'Atypicality', 'Scatterplot'])

        #  ================== PATIENT STAGING =================================================================================================================

        with col_staging_options:
            st.subheader('Style the graphs') 
            opacity = st.number_input('Select opacity value:', value=0.8)
            title_font = st.number_input('Title font:',value=34)
            title_axes = st.number_input('Axis labels font:',value=18)
            title_ticks = st.number_input('Ticks size font:',value=14)
            title_legend = st.number_input('Legend font:',value=22)
            font_list = [title_font, title_axes, title_ticks, title_legend]


        if chosen_plot =='Patient Staging':

            with col_staging_options:

                num_bins = st.number_input('Select the number of bins', value = 10)
                bin_width = st.number_input('Select bin width:', value = 0.04)

                for idx, label in enumerate(diagnosis_labels):
                        if label != 'CN':
                            color = st.text_input(f'Select color for {label}', value = f'{color_list[idx]}',placeholder='e.g. #000000')
                            color_diagnosis.append(color)

                # BARPLOT
                plot_staging = patient_staging(S=S,
                                        diagnosis=diagnosis, 
                                        color_list = color_diagnosis,
                                        num_bins=num_bins, 
                                        bin_width=bin_width,
                                        width = chosen_width,
                                        height = chosen_height,
                                        fontsize=font_list,
                                        opacity=opacity)

                # BOX
                box_staging = staging_boxes(S=S,
                                        diagnosis=diagnosis,
                                        color_list=color_diagnosis,
                                        width=chosen_width,
                                        fontsize=font_list)

            with col_staging:
                st.plotly_chart(plot_staging)
                st.plotly_chart(box_staging)  

        elif chosen_plot =='Atypicality':   

            with col_staging_options:

                num_bins = st.number_input('Select the number of bins', value = 15)
                bin_width = st.number_input('Select bin width:', value = 1.2)

                for idx, label in enumerate(diagnosis_labels):
                    if label != 'CN':
                        color = st.text_input(f'Select color for {label}', value = f'{color_list[idx]}',placeholder='e.g. #000000')
                        color_diagnosis.append(color)

                plot_atypicality = atypicality(S=S,
                                diagnosis=diagnosis, 
                                color_list = color_diagnosis,
                                num_bins=num_bins, 
                                bin_width=bin_width,
                                width = chosen_width,
                                height = chosen_height,
                                fontsize=font_list,
                                opacity=opacity)


                box_atypicality = atypicality_boxes(S=S,
                                        diagnosis=diagnosis,
                                        color_list=color_diagnosis,
                                        width=chosen_width,
                                        fontsize=font_list)

            with col_staging:
                    st.plotly_chart(plot_atypicality)
                    st.plotly_chart(box_atypicality)

        else:

            with col_staging_options:

                for idx, label in enumerate(diagnosis_labels):
                    if label != 'CN':
                        color = st.text_input(f'Select color for {label}', value = f'{color_list[idx]}',placeholder='e.g. #000000')
                        color_diagnosis.append(color)

                plot_scatter = staging_scatterplot(S=S,
                                                diagnosis=diagnosis,
                                                color_list=color_diagnosis,
                                                width = chosen_width,
                                                height = chosen_height,)

            with col_staging:
                st.plotly_chart(plot_scatter)

        # ADD DIVIDER
        st.markdown('---')

main()









